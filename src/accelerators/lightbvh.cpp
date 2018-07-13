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
		void InitLeaf(int first, int n, const Bounds3f &b_w, const Bounds_o &b_o, float e, int num) {
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
			lightNum = num;
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
		int splitAxis, firstLightOffset, nLights, lightNum;
		float energy;
		// TODO: number of emitters under this node required?
	};

	// LightBVHAccel Method Definitions
	LightBVHAccel::LightBVHAccel(std::vector<std::shared_ptr<Light>> l, float splitThreshold)
		: lights(std::move(l)), splitThreshold(splitThreshold) {
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
		// root->PrintNode(0, lightInfo);
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
		int nLights = end - start;

		// index of the last element of the first part and partial sort for median axis
		int mid = (end + start - 1) / 2;
		std::nth_element(&lightInfo[start], &lightInfo[mid],
			&lightInfo[end - 1] + 1,
			[dim](const LightBVHLightInfo &a,
				const LightBVHLightInfo &b) {
			return a.bounds_o.axis[dim] <
				b.bounds_o.axis[dim];
		});

		// calculate theta_o and theta_e for the node before the possible split
		Vector3f axis = lightInfo[mid].bounds_o.axis;
		float maxE = 0.f;
		float maxO = 0.f;
		calculateThetas(lightInfo, start, end - 1, axis, &maxO, &maxE);
		float totalAngle = 2 * Pi * (1 - cos(maxO) + (2 * sin(maxO) * maxE + cos(maxO) - cos(2 * maxE + maxO)) / 4);
		Bounds_o bounds_o = Bounds_o(axis, maxO, maxE);

		// total energy for the node
		float totalEnergy = 0;

		// only one light --> leaf creation
		if (nLights == 1) {
			// Create leaf _LightBVHBuildNode_
			int firstPrimOffset = orderedLights.size();
			int lightNum = lightInfo[start].lightNumber;
			orderedLights.push_back(lights[lightNum]);
			node->InitLeaf(firstPrimOffset, nLights, lightInfo[start].bounds_w, bounds_o, lightInfo[start].energy, lightNum);
			return node;
		}
		else {
			// two lights --> split in the middle
			if (nLights <= 2) {
				// Partition primitives into equally-sized subsets
				std::nth_element(&lightInfo[start], &lightInfo[mid],
					&lightInfo[end - 1] + 1,
					[dim](const LightBVHLightInfo &a,
						const LightBVHLightInfo &b) {
					return a.centroid[dim] <
						b.centroid[dim];
				});
				totalEnergy = lightInfo[start].energy + lightInfo[start + 1].energy;
			}
			// more than two lights --> calculating best split
			else {
				// split the lights in a certain number of intervals with the same size and compute the cost for each split
				// at the moment we are using the amount of lights as the number of buckets but I do not know how well it scales
				int buckets = std::min(12, nLights);
				std::vector<float> cost(buckets - 1);
				// a float defining the amount of lights that goes in every bucket
				float lightsPerBucket = (float)nLights / buckets;
				float totalCost;
				for (int i = 0; i < buckets - 1; i++) {
					// defines the end of the interval we are regarding right now, the int states the index of the last element of the first part
					int intervalEnd = start + lightsPerBucket * i;
					// first partial sort to find out which elements goes to which part, intervalEnd is the median for the axis with maximum extend
					std::nth_element(&lightInfo[(int)(start)], &lightInfo[intervalEnd],
						&lightInfo[end - 1] + 1,
						[dim](const LightBVHLightInfo &a,
							const LightBVHLightInfo &b) {
						return a.centroid[dim] <
							b.centroid[dim];
					});
					Bounds3f b0, b1;
					float leftE = 0, rightE = 0;
					float leftO = 0, rightO = 0;
					float leftEnergy = 0, rightEnergy = 0;
					// two integers for the median indices of left and right side
					int leftMid = start + (intervalEnd - start) / 2;
					int rightMid = intervalEnd + 1 + (end - 1 - intervalEnd) / 2;
					// partial sorting for both sides
					std::nth_element(&lightInfo[start], &lightInfo[leftMid],
						&lightInfo[intervalEnd] + 1,
						[dim](const LightBVHLightInfo &a,
							const LightBVHLightInfo &b) {
						return a.bounds_o.axis[dim] <
							b.bounds_o.axis[dim];
					});
					std::nth_element(&lightInfo[intervalEnd + 1], &lightInfo[rightMid],
						&lightInfo[end - 1] + 1,
						[dim](const LightBVHLightInfo &a,
							const LightBVHLightInfo &b) {
						return a.bounds_o.axis[dim] <
							b.bounds_o.axis[dim];
					});

					// setting the median axis
					Vector3f leftAxis = lightInfo[leftMid].bounds_o.axis;
					Vector3f rightAxis = lightInfo[rightMid].bounds_o.axis;

					// left side Theta_e and Theta_o calculations
					calculateThetas(lightInfo, start, intervalEnd, leftAxis, &leftO, &leftE);

					// left side bounds and energy calculations
					for (int j = start; j < intervalEnd + 1; j++) {
						LightBVHLightInfo l = lightInfo[j];
						b0 = Union(b0, l.bounds_w);
						leftEnergy += l.energy;
					}

					// right side Theta_e and Theta_o calculations
					calculateThetas(lightInfo, intervalEnd + 1, end - 1, rightAxis, &rightO, &rightE);

					// right side bounds and energy calculations
					for (int j = intervalEnd + 1; j < end; j++) {
						LightBVHLightInfo l = lightInfo[j];
						b1 = Union(b1, l.bounds_w);
						rightEnergy += l.energy;
					}

					// calculating the cost for the current split

					float leftAngle = 2 * Pi * (1 - cos(leftO) + (2 * sin(leftO) * leftE + cos(leftO) - cos(2 * leftE + leftO)) / 4);
					float rightAngle = 2 * Pi * (1 - cos(rightO) + (2 * sin(rightO) * rightE + cos(rightO) - cos(2 * rightE + rightO)) / 4);
					float leftCost = b0.SurfaceArea() * leftAngle * leftEnergy;
					float rightCost = b1.SurfaceArea() * rightAngle * rightEnergy;
					if (i == 0) {
						totalEnergy = leftEnergy + rightEnergy;
						Bounds3f totalBounds = Union(b0, b1);
						totalCost = totalBounds.SurfaceArea() * totalAngle * totalEnergy;
					}
					cost[i] = (leftCost + rightCost) / totalCost;
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

				// split lights at selected SAOH bucket, mid is the last element of the first part
				mid = start + minCostSplitBucket * lightsPerBucket;
				std::nth_element(&lightInfo[(int)(start)], &lightInfo[mid],
					&lightInfo[end - 1] + 1,
					[dim](const LightBVHLightInfo &a,
						const LightBVHLightInfo &b) {
					return a.centroid[dim] <
						b.centroid[dim];
				});
			}
			// recursively build children
			node->InitInterior(dim,
				recursiveBuild(arena, lightInfo, start, mid + 1,
					totalNodes, orderedLights),
				recursiveBuild(arena, lightInfo, mid + 1, end,
					totalNodes, orderedLights), bounds_o, totalEnergy);
		}
		return node;
	}

	void LightBVHAccel::calculateThetas(std::vector<LightBVHLightInfo> &lightInfo, int startIndex, int endIndex, Vector3f axis, float *theta_o, float *theta_e) {
		for (int j = startIndex; j < endIndex + 1; j++) {
			LightBVHLightInfo l = lightInfo[j];
			float o = 0.f;
			float e = 0.f;
			float angle = acos(Dot(axis, l.bounds_o.axis));
			if ((o = angle + l.bounds_o.theta_o) > *theta_o) {
				*theta_o = o;
			}
			if ((e = angle + l.bounds_o.theta_e + l.bounds_o.theta_o) > *theta_e) {
				*theta_e = e;
			}
		}
		*theta_e -= *theta_o;
	}

	std::shared_ptr<LightBVHAccel> CreateLightBVHAccelerator(
		std::vector<std::shared_ptr<Light>> lights, float splitThreshold) {
		return std::make_shared<LightBVHAccel>(std::move(lights), splitThreshold);
	}

	int LightBVHAccel::SampleOneLight(const Interaction &it, Sampler &sampler, float *pdf) {
		float sample1D = sampler.Get1D();
		*pdf = 1;
		return TraverseNodeForOneLight(root, sample1D, it, pdf);
	}

	int LightBVHAccel::TraverseNodeForOneLight(LightBVHNode *node, float sample1D, const Interaction &it, float *pdf) {
		// im already at a leaf of the tree
		if (node->nLights == 1) {
			return node->lightNum;
		}
		// finding out if I have to take the left or the right path
		float firstImportance = calculateImportance(it, node->children[0]);
		float secondImportance = calculateImportance(it, node->children[1]);
		// normalize the importance
		float totalImportance = firstImportance + secondImportance;
		// return -1 when the contribution of the sampled light will be zero (because of orientation)
		if (totalImportance == 0) {
			return -1;
		}
		firstImportance /= totalImportance;
		secondImportance /= totalImportance;
		// left child traversal
		if (sample1D < firstImportance) {
			*pdf *= firstImportance;
			return TraverseNodeForOneLight(node->children[0], sample1D / firstImportance, it, pdf);
		}
		// right child traversal
		*pdf *= secondImportance;
		return TraverseNodeForOneLight(node->children[1], (sample1D - firstImportance) / secondImportance, it, pdf);
	}

	std::vector<std::pair<int, float>> LightBVHAccel::SampleMultipleLights(const Interaction &it, Sampler &sampler) {
		float sample1D = sampler.Get1D();
		std::vector<std::pair<int, float>> lightVector;
		TraverseNodeForMultipleLights(root, sample1D, it, &lightVector);
		return lightVector;
	}

	void LightBVHAccel::TraverseNodeForMultipleLights(LightBVHNode *node, float sample1D, const Interaction &it, std::vector<std::pair<int, float>> *lightVector) {
		// im already at a leaf of the tree
		if (node->nLights == 1) {
			lightVector->push_back(std::pair<int, float>(node->lightNum, 1.f));
			return;
		}
		bool split = false;
		Bounds3f b = node->bounds_w;
		Point3f pMax = b.pMax;
		Point3f pMin = b.pMin;
		Point3f o = it.p;
		if (o.x >= pMin.x && o.x <= pMax.x && o.y >= pMin.y && o.y <= pMax.y && o.z >= pMin.z && o.z <= pMax.z) {
			split = true;
		}
		else {
			float maxAngle = 0;
			Vector3f c[8];
			for (int i = 0; i < 8; i++) {
				c[i] = Normalize(b.Corner(i) - o);
			}
			for (int i = 0; i < 7; i++) {
				for (int j = i + 1; j < 8; j++) {
					maxAngle = std::max(maxAngle, std::acos(Dot(c[i], c[j])));
				}
			}
			if (maxAngle / Pi > splitThreshold) {
				split = true;
			}
		}
		if (split) {
			TraverseNodeForMultipleLights(node->children[0], sample1D, it, lightVector);
			TraverseNodeForMultipleLights(node->children[1], sample1D, it, lightVector);
		}
		else {
			float pdf = 1.f;
			int lightNum = TraverseNodeForOneLight(node, sample1D, it, &pdf);
			lightVector->push_back(std::pair<int, float>(lightNum, pdf));
		}
	}

	float LightBVHAccel::calculateImportance(const Interaction &it, LightBVHNode* node) {
		Point3f o = it.p;
		// turn normal if pointed to wrong side
		Normal3f n = it.n;
		if (Dot(it.wo, n) < 0) {
			n *= -1;
		}
		// check if the bounding box of node is behind the shaded point
		bool wrongDirection = true;
		if (it.IsSurfaceInteraction()) {
			for (int i = 0; i < 8; i++) {
				// no need to normalize since just comparing to 0
				if (Dot(n, node->bounds_w.Corner(i) - o) >= 0 ) {
					wrongDirection = false;
					break;
				}
			}
			if (wrongDirection) {
				return 0;
			}
		}
		float theta_e = node->bounds_o.theta_e;
		float theta_o = node->bounds_o.theta_o;
		Vector3f d = node->centroid - o;
		float distance = d.Length();
		d = Normalize(d);
		float theta = acos(Dot(node->bounds_o.axis, -d));
		float theta_u = 0;
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

		// shading point is already in the box -> can always find a theta_u with 1 --> angleImportance is always 1
		float angleImportance = 1;
		// shading point is not in the box
		if (t0 > 0 && t1 > 0) {
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
				theta_u = std::max(theta_u, acos(Dot(d, Normalize(c[i] - o))));
			}
			// experimental: test if point is in the cone of the sampled light
			//if (theta - theta_o - theta_u - theta_e > 0) {
			angleImportance = cos(std::max(0.f, std::min(theta - theta_o - theta_u, theta_e)));
			//}
		}
		return node->energy * angleImportance / (distance * distance);
	}
}