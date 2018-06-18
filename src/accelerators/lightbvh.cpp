// accelerators/lightbvh.cpp*
#include "accelerators/lightbvh.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>
#include "core/light.h"
#include <iostream>
#include "sampler.h"

namespace pbrt {

	STAT_MEMORY_COUNTER("Memory/LightBVH tree", treeBytes);
	STAT_RATIO("LightBVH/Lights per leaf node", totalLights, totalLeafNodes);
	STAT_COUNTER("LightBVH/Interior nodes", interiorNodes);
	STAT_COUNTER("LightBVH/Leaf nodes", leafNodes);

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
	struct LightBVHLightInfo {

		size_t lightNumber;
		Bounds3f bounds_w;
		Bounds_o bounds_o;
		Point3f centroid;
		float energy = 0.f;

		LightBVHLightInfo() {}

		LightBVHLightInfo(size_t lightNumber, Light &light)
			: lightNumber(lightNumber) {
			bounds_o = Bounds_o(light.Axis(), light.Theta_o(), light.Theta_e());
			bounds_w = light.Bounds();
			centroid = .5f * bounds_w.pMin + .5f * bounds_w.pMax;
			Spectrum s = light.Power();
			int n = s.nSamples;
			for (int i = 0; i < n; i++) {
				energy += s[i];
			}
		}
	};

	struct LightBVHNode {
		// LightBVHNode Public Methods
		void InitLeaf(int first, int n, const Bounds3f &b_w, const Bounds_o &b_o, float e) {
			firstLightOffset = first;
			nLights = n;
			bounds_w = b_w;
			bounds_o = b_o;
			children[0] = children[1] = nullptr;
			++leafNodes;
			++totalLeafNodes;
			totalLights += n;
			energy = e;
			centroid = .5f * bounds_w.pMin + .5f * bounds_w.pMax;
		}
		void InitInterior(int split, LightBVHNode *c0, LightBVHNode *c1, const Bounds_o &b_o, float e) {
			children[0] = c0;
			children[1] = c1;
			bounds_w = Union(c0->bounds_w, c1->bounds_w);
			bounds_o = b_o;
			splitAxis = split;
			nLights = 0;
			++interiorNodes;
			energy = e;
			centroid = .5f * bounds_w.pMin + .5f * bounds_w.pMax;
		}
		void PrintNode(int depth, std::vector<LightBVHLightInfo> &lightInfo) {
			std::string s;
			for (int i = 0; i < depth; i++)  s += "-";
			if (nLights != 0) {
				for (int i = 0; i < nLights; i++) {
					Point3f p = centroid;
					std::cout << s << "n = " << nLights << ", " << "x = " << p.x << ", y = " << p.y << ", z = " << p.z << std::endl;
				}
			}
			else {
				std::cout << s << std::endl;
				children[0]->PrintNode(depth + 1, lightInfo);
				children[1]->PrintNode(depth + 1, lightInfo);
			}
		}
		Bounds3f bounds_w;
		Bounds_o bounds_o;
		Point3f centroid;
		LightBVHNode *children[2];
		int splitAxis, firstLightOffset, nLights;
		float energy;
		// TODO: number of emitters under this node required?
	};

	// LightBVHAccel Method Definitions
	LightBVHAccel::LightBVHAccel(std::vector<std::shared_ptr<Light>> l,
		int maxPrimsInNode)
		: maxLightsInNode(std::min(255, maxPrimsInNode)),
		lights(std::move(l)) {
		ProfilePhase _(Prof::AccelConstruction);
		if (lights.empty()) return;
		// Build LBVH from _lights_

		// Initialize _lightInfo_ array for lights
		std::vector<LightBVHLightInfo> lightInfo(lights.size());
		for (size_t i = 0; i < lights.size(); ++i)
			lightInfo[i] = { i, *lights[i] };

		// Build LightBVH tree for lights using _lightInfo_
		int totalNodes = 0;
		std::vector<std::shared_ptr<Light>> orderedLights;
		orderedLights.reserve(lights.size());
		MemoryArena arena(1024 * 1024);
		root = recursiveBuild(arena, lightInfo, 0, lights.size(), &totalNodes, orderedLights);
		root->PrintNode(0, lightInfo);
		lights.swap(orderedLights);
		lightInfo.resize(0);
	}

	LightBVHNode* LightBVHAccel::recursiveBuild(
		MemoryArena &arena, std::vector<LightBVHLightInfo> &lightInfo, int start,
		int end, int *totalNodes,
		std::vector<std::shared_ptr<Light>> &orderedLights) {
		CHECK_NE(start, end);
		// TODO: Change malloc back to arena allocation when I know what the problem is
		LightBVHNode* node = (LightBVHNode*)malloc(sizeof LightBVHNode);
		(*totalNodes)++;
		// Compute bound of light centroids, choose split dimension _dim_
		Bounds3f centroidBounds;
		for (int i = start; i < end; ++i)
			centroidBounds = Union(centroidBounds, lightInfo[i].centroid);
		int dim = centroidBounds.MaximumExtent();
		// Compute bounds of all lights in BVH node
		Bounds3f totalBounds_w;
		int mid = (start + end) / 2;
		std::nth_element(&lightInfo[start], &lightInfo[mid],
			&lightInfo[end - 1] + 1,
			[dim](const LightBVHLightInfo &a,
				const LightBVHLightInfo &b) {
			return a.bounds_o.axis[dim] <
				b.bounds_o.axis[dim];
		});
		Vector3f axis = lightInfo[mid].bounds_o.axis;
		float maxE = 0.f;
		float maxO = 0.f;
		for (int i = start; i < end; ++i) {
			totalBounds_w = Union(totalBounds_w, lightInfo[i].bounds_w);
			float e = 0.f;
			float o = 0.f;
			if ((e = acos(AbsDot(axis, lightInfo[i].bounds_o.axis)) + lightInfo[i].bounds_o.theta_e) > maxE) {
				maxE = e;
			}
			if ((o = e + lightInfo[i].bounds_o.theta_o) > maxO) {
				maxO = o;
			}
		}
		Bounds_o bounds_o = Bounds_o(axis, maxE, maxO);
		float totalAngle = 2 * Pi * (1 - cos(maxO) + (2 * sin(maxO) * maxE + cos(maxO) - cos(2 * maxE + maxO)) / 4);
		int nLights = end - start;
		// total energy for the node
		float totalEnergy = 0;
		if (nLights == 1) {
			// Create leaf _LBVHBuildNode_
			int firstPrimOffset = orderedLights.size();
			for (int i = start; i < end; ++i) {
				int lightNum = lightInfo[i].lightNumber;
				orderedLights.push_back(lights[lightNum]);
			}
			node->InitLeaf(firstPrimOffset, nLights, totalBounds_w, bounds_o, lightInfo[0].energy);
			return node;
		}
		else {
			// Partition lights into two sets and build children
			// Partition primitives using approximate SAOH
			if (nLights <= 2) {
				// Partition primitives into equally-sized subsets
				std::nth_element(&lightInfo[start], &lightInfo[mid],
					&lightInfo[end - 1] + 1,
					[dim](const LightBVHLightInfo &a,
						const LightBVHLightInfo &b) {
					return a.centroid[dim] <
						b.centroid[dim];
				});
				// TODO: calculate totalEnergy?
			}
			else {
				// split the lights in 12 intervals with the same size and compute the cost for each split
				// 11 splitindices = 12 splits
				int buckets = std::min(12, nLights);
				std::vector<float> cost(buckets - 1);
				float lightsPerInterval = (float)nLights / buckets;
				for (int i = 0; i < buckets - 1; i++) {

					std::nth_element(&lightInfo[(int)(lightsPerInterval * i)], &lightInfo[(int)(lightsPerInterval * (i + 1))],
						&lightInfo[end - 1] + 1,
						[dim](const LightBVHLightInfo &a,
							const LightBVHLightInfo &b) {
						return a.centroid[dim] <
							b.centroid[dim];
					});
				}
				for (int i = 0; i < buckets - 1; i++) {
					Bounds3f b0, b1;
					float leftmaxE = 0, rightmaxE = 0;
					float leftmaxO = 0, rightmaxO = 0;
					float leftEnergy = 0, rightEnergy = 0;
					int intervalEnd = (int)(lightsPerInterval * (i + 1));
					std::nth_element(&lightInfo[start], &lightInfo[start + intervalEnd / 2],
						&lightInfo[start + intervalEnd] + 1,
						[dim](const LightBVHLightInfo &a,
							const LightBVHLightInfo &b) {
						return a.bounds_o.axis[dim] <
							b.bounds_o.axis[dim];
					});
					Vector3f leftAxis = lightInfo[mid].bounds_o.axis;
					std::nth_element(&lightInfo[intervalEnd + 1], &lightInfo[intervalEnd + 1 + (end - 1 - (intervalEnd + 1)) / 2],
						&lightInfo[end - 1] + 1,
						[dim](const LightBVHLightInfo &a,
							const LightBVHLightInfo &b) {
						return a.bounds_o.axis[dim] <
							b.bounds_o.axis[dim];
					});
					Vector3f rightAxis = lightInfo[intervalEnd + 1 + (end - 1 - (intervalEnd + 1)) / 2].bounds_o.axis;
					for (int j = 0; j <= (i + 1) * lightsPerInterval; j++) {
						LightBVHLightInfo l = lightInfo[start + j];
						b0 = Union(b0, l.bounds_w);
						float e = 0.f;
						float o = 0.f;
						if ((e = acos(AbsDot(leftAxis, l.bounds_o.axis)) + l.bounds_o.theta_e) > leftmaxE) {
							leftmaxE = e;
						}
						if ((o = e + l.bounds_o.theta_o) > leftmaxO) {
							leftmaxO = o;
						}
						leftEnergy += l.energy;
					}
					for (int j = (i + 1) * lightsPerInterval + 1; j < nLights; j++) {
						LightBVHLightInfo l = lightInfo[start + j];
						b1 = Union(b1, l.bounds_w);
						float e = 0.f;
						float o = 0.f;
						if ((e = acos(AbsDot(rightAxis, l.bounds_o.axis)) + l.bounds_o.theta_e) > rightmaxE) {
							rightmaxE = e;
						}
						if ((o = e + l.bounds_o.theta_o) > rightmaxO) {
							rightmaxO = o;
						}
						rightEnergy += l.energy;
					}
					totalEnergy = leftEnergy + rightEnergy;
					float leftAngle = 2 * Pi * (1 - cos(leftmaxO) + (2 * sin(leftmaxO) * leftmaxE + cos(leftmaxO) - cos(2 * leftmaxE + leftmaxO)) / 4);
					float rightAngle = 2 * Pi * (1 - cos(rightmaxO) + (2 * sin(rightmaxO) * rightmaxE + cos(rightmaxO) - cos(2 * rightmaxE + rightmaxO)) / 4);
					cost[i] = (b0.SurfaceArea() * leftAngle * leftEnergy +
						b1.SurfaceArea() * rightAngle * rightEnergy) /
						totalBounds_w.SurfaceArea() * totalAngle * totalEnergy;
				}


				// Find bucket to split at that minimizes SAOH metric
				Float minCost = cost[0];
				int minCostSplitBucket = 0;
				for (int i = 1; i < buckets - 1; i++) {
					if (cost[i] < minCost) {
						minCost = cost[i];
						minCostSplitBucket = i;
					}
				}

				// split lights at selected SAOH bucket
				mid = start + (int)((minCostSplitBucket + 1) * lightsPerInterval) + 1;
			}
			node->InitInterior(dim,
				recursiveBuild(arena, lightInfo, start, mid,
					totalNodes, orderedLights),
				recursiveBuild(arena, lightInfo, mid, end,
					totalNodes, orderedLights), bounds_o, totalEnergy);
		}
		return node;
	}

	std::shared_ptr<LightBVHAccel> CreateLightBVHAccelerator(
		std::vector<std::shared_ptr<Light>> lights, const ParamSet &ps) {
		// TODO: max lights in node?
		int maxLightsInNode = 1;
		return std::make_shared<LightBVHAccel>(std::move(lights), maxLightsInNode);
	}

	int LightBVHAccel::Sample(const Interaction &it, const Scene &scene, Sampler &sampler,
		bool handleMedia, const Distribution1D *lightDistrib) {
		float sample1D = sampler.Get1D();
		return TraverseNode(root, sample1D, it, scene, sampler, handleMedia, lightDistrib);
	}

	int LightBVHAccel::TraverseNode(LightBVHNode *node, float sample1D, const Interaction &it, const Scene &scene, Sampler &sampler, bool handleMedia, const Distribution1D *lightDistrib) {
		// im already at a leaf of the tree
		std::cout << node->energy << "stfu";
		if (node->nLights != 0) {
			return floor(sample1D * node->nLights);
		}
		// finding out if I have to take the left or the right path
		Point3f o = it.p;
		float firstImportance = calculateImportance(o, node->children[0]);
		float secondImportance = calculateImportance(o, node->children[1]);
		// normalize the importance
		firstImportance = firstImportance / (firstImportance + secondImportance);
		secondImportance = secondImportance / (firstImportance + secondImportance);
		// left side traversal
		if (sample1D < firstImportance) {
			return TraverseNode(node, sample1D / firstImportance, it, scene, sampler, handleMedia, lightDistrib);
		}
		// right side traversal
		return TraverseNode(node, (sample1D - firstImportance) / secondImportance, it, scene, sampler, handleMedia, lightDistrib);
	}

	float LightBVHAccel::calculateImportance(Point3f o, LightBVHNode* node) {
		float theta_e = node->bounds_o.theta_e;
		float theta_o = node->bounds_o.theta_o;
		Vector3f d = node->centroid - o;
		float distance = d.Length();
		float theta = acos(Dot(node->bounds_o.axis, -d));
		float theta_u;
		// numeric inaccuracies?
		float ep = 0.0001;

		float t0 = -std::numeric_limits<float>::infinity(), t1 = std::numeric_limits<float>::infinity();
		for (int i = 0; i < 3; ++i) {
			// Update interval for _i_th bounding box slab
			float invRayDir = 1 / d[i];
			float tNear = (node->bounds_w.pMin[i] - o[i]) * invRayDir;
			float tFar = (node->bounds_w.pMax[i] - o[i]) * invRayDir;

			// Update parametric interval from slab intersection $t$ values
			if (tNear > tFar) std::swap(tNear, tFar);

			// Update _tFar_ to ensure robust ray--bounds intersection
			tFar *= 1 + 2 * gamma(3);
			t0 = tNear > t0 ? tNear : t0;
			t1 = tFar < t1 ? tFar : t1;
		}

		// shading point is already in the box -> can always find a theta_u with 0
		if (t0 < 0 || t1 < 0) {
			theta_u = 0;
		}
		else {
			Point3f is = o + d * t0;
			int side;
			Point3f pMin = node->bounds_w.pMin;
			Point3f pMax = node->bounds_w.pMax;
			for (int i = 0; i < 3; i++) {
				if (abs(is[i] - node->bounds_w.pMin[i]) < ep) {
					side = i * 2;
				}
				else if (abs(is[i] - node->bounds_w.pMax[i]) < ep) {
					side = i * 2 + 1;
				}
			}
			Point3f c[4];
			// ugly code defining the corners to check
			switch (side) {
			case 0:
				c[0] = Point3f(pMin.x, pMin.y, pMin.z);
				c[1] = Point3f(pMin.x, pMin.y, pMax.z);
				c[2] = Point3f(pMin.x, pMax.y, pMin.z);
				c[3] = Point3f(pMin.x, pMax.y, pMax.z);
				break;
			case 1:
				c[0] = Point3f(pMax.x, pMin.y, pMin.z);
				c[1] = Point3f(pMax.x, pMin.y, pMax.z);
				c[2] = Point3f(pMax.x, pMax.y, pMin.z);
				c[3] = Point3f(pMax.x, pMax.y, pMax.z);
				break;
			case 2:
				c[0] = Point3f(pMin.x, pMin.y, pMin.z);
				c[1] = Point3f(pMax.x, pMin.y, pMin.z);
				c[2] = Point3f(pMin.x, pMin.y, pMax.z);
				c[3] = Point3f(pMax.x, pMin.y, pMax.z);
				break;
			case 3:
				c[0] = Point3f(pMin.x, pMax.y, pMin.z);
				c[1] = Point3f(pMax.x, pMax.y, pMin.z);
				c[2] = Point3f(pMin.x, pMax.y, pMax.z);
				c[3] = Point3f(pMax.x, pMax.y, pMax.z);
				break;
			case 4:
				c[0] = Point3f(pMin.x, pMin.y, pMin.z);
				c[1] = Point3f(pMax.x, pMin.y, pMin.z);
				c[2] = Point3f(pMin.x, pMax.y, pMin.z);
				c[3] = Point3f(pMax.x, pMax.y, pMin.z);
				break;
			case 5:
				c[0] = Point3f(pMin.x, pMin.y, pMax.z);
				c[1] = Point3f(pMax.x, pMin.y, pMax.z);
				c[2] = Point3f(pMin.x, pMax.y, pMax.z);
				c[3] = Point3f(pMax.x, pMax.y, pMax.z);
				break;
			default: std::cout << "There has been an issue finding the side of the bounds that was hit";
			}
			// setting theta_u to the maximum angle
			for (int i = 0; i < 4; i++) {
				theta_u = std::max(theta_u, acos(Dot(d, c[i] - o)));
			}
		}

		return node->energy * cos(std::max(0.f, std::min(theta - theta_o - theta_u, theta_e))) / (distance * distance);
	}
}