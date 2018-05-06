// accelerators/lbvh.cpp*
#include "accelerators/lbvh.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>
#include "core/light.h"

namespace pbrt {

	STAT_MEMORY_COUNTER("Memory/LBVH tree", treeBytes);
	STAT_RATIO("LBVH/Lights per leaf node", totalLights, totalLeafNodes);
	STAT_COUNTER("LBVH/Interior nodes", interiorNodes);
	STAT_COUNTER("LBVH/Leaf nodes", leafNodes);

	struct LinearLBVHNode {
		Bounds3f bounds;
		union {
			int primitivesOffset;   // leaf
			int secondChildOffset;  // interior
		};
		uint16_t nPrimitives;  // 0 -> interior node
		uint8_t axis;          // interior node: xyz
		uint8_t pad[1];        // ensure 32 byte total size
	};

	struct Bounds_o {
		Vector3f axis;
		float theta_o;
		float theta_e;

		Bounds_o() {}

		Bounds_o(Vector3f axis, float theta_o, float theta_e)
			: axis(axis),
			theta_o(theta_o),
			theta_e(theta_e) {}
	};

	// BVHAccel Local Declarations
	struct LBVHLightInfo {

		size_t lightNumber;
		Bounds3f bounds_w;
		Bounds_o bounds_o;
		Point3f centroid;

		LBVHLightInfo() {}

		LBVHLightInfo(size_t lightNumber, Light &light)
			: lightNumber(lightNumber),
			centroid(.5f * bounds_w.pMin + .5f * bounds_w.pMax) {

			bounds_o = Bounds_o(light.Axis(), light.Theta_o(), light.Theta_e());
			bounds_w = light.Bounds();
		}
	};

	struct LBVHBuildNode {
		// LBVHBuildNode Public Methods
		void InitLeaf(int first, int n, const Bounds3f &b_w, const Bounds_o &b_o) {
			firstLightOffset = first;
			nLights = n;
			bounds_w = b_w;
			bounds_o = b_o;
			children[0] = children[1] = nullptr;
			++leafNodes;
			++totalLeafNodes;
			totalLights += n;
		}
		void InitInterior(int split, LBVHBuildNode *c0, LBVHBuildNode *c1, const Bounds_o &b_o) {
			children[0] = c0;
			children[1] = c1;
			bounds_w = Union(c0->bounds_w, c1->bounds_w);
			bounds_o = b_o;
			splitAxis = split;
			nLights = 0;
			++interiorNodes;
		}
		Bounds3f bounds_w;
		Bounds_o bounds_o;
		LBVHBuildNode *children[2];
		int splitAxis, firstLightOffset, nLights;
		// TODO: initializing energy
		float energy;
	};

	Bounds_o Union_o(const Bounds_o &b1, const Bounds_o &b2) {
		Vector3f axis;
		float theta_o;
		float theta_e;
		return Bounds_o(axis, theta_o, theta_e);
	}

	// LBVHAccel Method Definitions
	LBVHAccel::LBVHAccel(std::vector<std::shared_ptr<Light>> l,
		int maxPrimsInNode)
		: maxLightsInNode(std::min(255, maxPrimsInNode)),
		lights(std::move(l)) {
		ProfilePhase _(Prof::AccelConstruction);
		if (lights.empty()) return;
		// Build LBVH from _lights_

		// Initialize _lightInfo_ array for lights
		std::vector<LBVHLightInfo> lightInfo(lights.size());
		for (size_t i = 0; i < lights.size(); ++i)
			lightInfo[i] = { i, *lights[i] };

		// Build LBVH tree for lights using _lightInfo_
		MemoryArena arena(1024 * 1024);
		int totalNodes = 0;
		std::vector<std::shared_ptr<Light>> orderedLights;
		orderedLights.reserve(lights.size());
		LBVHBuildNode *root;
		root = recursiveBuild(arena, lightInfo, 0, lights.size(), &totalNodes, orderedLights);
		lights.swap(orderedLights);
		lightInfo.resize(0);
		LOG(INFO) << StringPrintf("BVH created with %d nodes for %d "
			"primitives (%.2f MB), arena allocated %.2f MB",
			totalNodes, (int)lights.size(),
			float(totalNodes * sizeof(LinearLBVHNode)) /
			(1024.f * 1024.f),
			float(arena.TotalAllocated()) /
			(1024.f * 1024.f));

		// Compute representation of depth-first traversal of LBVH tree
		treeBytes += totalNodes * sizeof(LinearLBVHNode) + sizeof(*this) +
			lights.size() * sizeof(lights[0]);
		nodes = AllocAligned<LinearLBVHNode>(totalNodes);
		int offset = 0;
		flattenLBVHTree(root, &offset);
		CHECK_EQ(totalNodes, offset);
	}

	struct BucketInfo {
		int count = 0;
		Bounds3f bounds_w;
		Bounds_o bounds_o;
	};

	LBVHBuildNode *LBVHAccel::recursiveBuild(
		MemoryArena &arena, std::vector<LBVHLightInfo> &lightInfo, int start,
		int end, int *totalNodes,
		std::vector<std::shared_ptr<Light>> &orderedLights) {
		CHECK_NE(start, end);
		LBVHBuildNode *node = arena.Alloc<LBVHBuildNode>();
		(*totalNodes)++;
		// Compute bounds of all lights in BVH node
		Bounds3f bounds_w;
		Vector3f axis = lightInfo[(start + end) / 2].bounds_o.axis;
		float maxE = 0.f;
		float maxO = 0.f;
		for (int i = start; i < end; ++i) {
			bounds_w = Union(bounds_w, lightInfo[i].bounds_w);
			float e = 0.f;
			float o = 0.f;
			if ((e = AbsDot(axis, lightInfo[i].bounds_o.axis) + lightInfo[i].bounds_o.theta_e) > maxE) {
				maxE = e;
			}
			if ((o = e + lightInfo[i].bounds_o.theta_o) > maxO) {
				maxO = o;
			}
		}
		Bounds_o bounds_o = Bounds_o(axis, maxE, maxO);
		int nLights = end - start;
		if (nLights == 1) {
			// Create leaf _LBVHBuildNode_
			int firstPrimOffset = orderedLights.size();
			for (int i = start; i < end; ++i) {
				int lightNum = lightInfo[i].lightNumber;
				orderedLights.push_back(lights[lightNum]);
			}
			node->InitLeaf(firstPrimOffset, nLights, bounds_w, bounds_o);
			return node;
		}
		else {
			// Compute bound of light centroids, choose split dimension _dim_
			Bounds3f centroidBounds;
			for (int i = start; i < end; ++i)
				centroidBounds = Union(centroidBounds, lightInfo[i].centroid);
			int dim = centroidBounds.MaximumExtent();

			// Partition lights into two sets and build children
			int mid = (start + end) / 2;
			if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
				// Create leaf _LBVHBuildNode_
				int firstPrimOffset = orderedLights.size();
				for (int i = start; i < end; ++i) {
					int lightNum = lightInfo[i].lightNumber;
					orderedLights.push_back(lights[lightNum]);
				}
				node->InitLeaf(firstPrimOffset, nLights, bounds_w, bounds_o);
				return node;
			}
			else {
				// Partition primitives using approximate SAOH
				if (nLights <= 2) {
					// Partition primitives into equally-sized subsets
					mid = (start + end) / 2;
					std::nth_element(&lightInfo[start], &lightInfo[mid],
						&lightInfo[end - 1] + 1,
						[dim](const LBVHLightInfo &a,
							const LBVHLightInfo &b) {
						return a.centroid[dim] <
							b.centroid[dim];
					});
				}
				else {
					// Allocate _BucketInfo_ for SAH partition buckets
					PBRT_CONSTEXPR int nBuckets = 12;
					BucketInfo buckets[nBuckets];

					// Initialize _BucketInfo_ for SAH partition buckets
					for (int i = start; i < end; ++i) {
						int b = nBuckets *
							centroidBounds.Offset(
								lightInfo[i].centroid)[dim];
						if (b == nBuckets) b = nBuckets - 1;
						CHECK_GE(b, 0);
						CHECK_LT(b, nBuckets);
						buckets[b].count++;
						buckets[b].bounds_w =
							Union(buckets[b].bounds_w, lightInfo[i].bounds_w);
					}

					// Compute costs for splitting after each bucket
					Float cost[nBuckets - 1];
					for (int i = 0; i < nBuckets - 1; ++i) {
						Bounds3f b0, b1;
						int count0 = 0, count1 = 0;
						float leftmaxE = 0, rightmaxE = 0;
						float leftmaxO = 0, rightmaxO = 0;
						for (int j = 0; j <= i; ++j) {
							b0 = Union(b0, buckets[j].bounds_w);
							count0 += buckets[j].count;
							float e = 0.f;
							float o = 0.f;
							if ((e = AbsDot(axis, lightInfo[j].bounds_o.axis) + lightInfo[j].bounds_o.theta_e) > maxE) {
								leftmaxE = e;
							}
							if ((o = e + lightInfo[j].bounds_o.theta_o) > maxO) {
								leftmaxO = o;
							}
						}
						for (int j = i + 1; j < nBuckets; ++j) {
							b1 = Union(b1, buckets[j].bounds_w);
							count1 += buckets[j].count;
							float e = 0.f;
							float o = 0.f;
							if ((e = AbsDot(axis, lightInfo[j].bounds_o.axis) + lightInfo[j].bounds_o.theta_e) > maxE) {
								leftmaxE = e;
							}
							if ((o = e + lightInfo[j].bounds_o.theta_o) > maxO) {
								leftmaxO = o;
							}
						}
						// not sure if the integral is correct
						float leftAngle = 2 * Pi * (1 - cos(leftmaxO) + (2 * sin(leftmaxO) * leftmaxE + cos(leftmaxO - cos(2 * leftmaxE + leftmaxO))));
						float rightAngle = 2 * Pi * (1 - cos(rightmaxO) + (2 * sin(rightmaxO) * rightmaxE + cos(rightmaxO - cos(2 * rightmaxE + rightmaxO))));
						// TODO: missing energy
						cost[i] = 1 +
							(count0 * b0.SurfaceArea() * leftAngle +
								count1 * b1.SurfaceArea() * rightAngle) /
							bounds_w.SurfaceArea();
					}

					// Find bucket to split at that minimizes SAH metric
					Float minCost = cost[0];
					int minCostSplitBucket = 0;
					for (int i = 1; i < nBuckets - 1; ++i) {
						if (cost[i] < minCost) {
							minCost = cost[i];
							minCostSplitBucket = i;
						}
					}

					// Either create leaf or split lights at selected SAH
					// bucket
					Float leafCost = nLights;
					if (nLights > maxLightsInNode || minCost < leafCost) {
						LBVHLightInfo *pmid = std::partition(
							&lightInfo[start], &lightInfo[end - 1] + 1,
							[=](const LBVHLightInfo &pi) {
							int b = nBuckets *
								centroidBounds.Offset(pi.centroid)[dim];
							if (b == nBuckets) b = nBuckets - 1;
							CHECK_GE(b, 0);
							CHECK_LT(b, nBuckets);
							return b <= minCostSplitBucket;
						});
						mid = pmid - &lightInfo[0];
					}
					else {
						// Create leaf _BVHBuildNode_
						int firstPrimOffset = orderedLights.size();
						for (int i = start; i < end; ++i) {
							int primNum = lightInfo[i].lightNumber;
							orderedLights.push_back(lights[primNum]);
						}
						node->InitLeaf(firstPrimOffset, nLights, bounds_w, bounds_o);
						return node;
					}
				}
				node->InitInterior(dim,
					recursiveBuild(arena, lightInfo, start, mid,
						totalNodes, orderedLights),
					recursiveBuild(arena, lightInfo, mid, end,
						totalNodes, orderedLights), bounds_o);
			}
		}
		return node;
	}
}