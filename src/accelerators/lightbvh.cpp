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
#include "integrator.h"

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
		void InitLeaf(int n, const Bounds3f &b_w, const Bounds_o &b_o, float e, int num) {
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
			nLights = c0->nLights + c1->nLights;
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
		float energy;
		LightBVHNode *children[2];
		int splitAxis, nLights, lightNum;
	};

	struct LinearLightBVHNode {
		Bounds3f bounds_w;
		Bounds_o bounds_o;
		float energy;
		Point3f centroid;
		union {
			int lightNum;   // leaf
			int secondChildOffset;  // interior
		};
		uint16_t nLights;  
		uint8_t splitAxis;          
	};


	// LightBVHAccel Method Definitions
	LightBVHAccel::LightBVHAccel(std::vector<std::shared_ptr<Light>> l, float splitThreshold)
		: lights(std::move(l)), splitThreshold(splitThreshold), rng(0) {
		ProfilePhase _(Prof::AccelConstruction);
		if (lights.empty()) return;
		// Build LBVH from _lights_

		// Initialize _lightInfo_ array for lights
		std::vector<LightBVHLightInfo> lightInfo(lights.size());
		for (size_t i = 0; i < lights.size(); ++i)
			lightInfo[i] = { i, *lights[i] };

		// Build LightBVH tree for lights using _lightInfo_
		int totalNodes = 0;
		MemoryArena arena(1024 * 1024);
		root = recursiveBuild(arena, lightInfo, 0, lights.size(), &totalNodes);
		// root->PrintNode(0, lightInfo);
		lightInfo.resize(0);
		nodes = AllocAligned<LinearLightBVHNode>(totalNodes);
		int offset = 0;
		flattenLightBVHTree(root, &offset);
		CHECK_EQ(totalNodes, offset);
	}

	LightBVHNode* LightBVHAccel::recursiveBuild(MemoryArena &arena, std::vector<LightBVHLightInfo> &lightInfo, int start, int end, int *totalNodes) {
		CHECK_NE(start, end);
		LightBVHNode *node = (LightBVHNode*)malloc(sizeof LightBVHNode);
		(*totalNodes)++;
		// number of lights under this node
		int nLights = end - start;
		// current middle element to split in two parts later
		int mid = (end + start - 1) / 2;
		// split the lights in a certain number of intervals with the same size and compute the cost for each split
		// at the moment we are using the amount of lights as the number of buckets but I do not know how well it scales
		int buckets = std::min(nLights, nLights);
		// cost vector for each split + dim
		std::vector<float> cost = std::vector<float>((buckets - 1) * 3);
		// a float defining the amount of lights that goes in every bucket
		float lightsPerBucket = (float)nLights / buckets;
		// the three orientation bounds of this node for the three dimensions
		Bounds_o bounds_o[3];
		// total energy for the node
		float totalEnergy = 0;
		// for each dimension calculate the splits
		for (int dim = 0; dim < 3; dim++) {
			//// index of the last element of the first part and partial sort for median axis
			//std::nth_element(&lightInfo[start], &lightInfo[mid],
			//	&lightInfo[end - 1] + 1,
			//	[dim](const LightBVHLightInfo &a,
			//		const LightBVHLightInfo &b) {
			//	return a.bounds_o.axis[dim] <
			//		b.bounds_o.axis[dim];
			//});

			Vector3f axis(0.f, 0.f, 0.f);
			for (int i = start; i < end; i++) {
				axis += lightInfo[i].bounds_o.axis;
			}
			if (axis == Vector3f(0.f, 0.f, 0.f)) {
				axis = lightInfo[start].bounds_o.axis;
			}
			else {
				axis = Normalize(axis);
			}

			// calculate theta_o and theta_e for the node before the possible split
			float maxE = 0.f;
			float maxO = 0.f;
			calculateThetas(lightInfo, start, end - 1, axis, &maxO, &maxE);
			float totalAngle = 2 * Pi * ((1 - cos(maxO)) + ((2 * sin(maxO) * maxE + cos(maxO) - cos(2 * maxE + maxO)) / 4));
			bounds_o[dim] = Bounds_o(axis, maxO, maxE);

			// only one light --> leaf creation
			if (nLights == 1) {
				// Create leaf _LightBVHBuildNode_
				int lightNum = lightInfo[start].lightNumber;
				node->InitLeaf(nLights, lightInfo[start].bounds_w, bounds_o[dim], lightInfo[start].energy, lightNum);
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
					float totalCost;
					for (int i = 0; i < buckets - 1; i++) {
						// defines the end of the interval we are regarding right now, the int states the index of the last element of the first part
						int intervalEnd = start + lightsPerBucket * (i + 1) - 1;
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
						//// two integers for the median indices of left and right side
						//int leftMid = start + (intervalEnd - start) / 2;
						//int rightMid = intervalEnd + 1 + (end - 1 - intervalEnd) / 2;
						//// partial sorting for both sides for median axis
						//std::nth_element(&lightInfo[start], &lightInfo[leftMid],
						//	&lightInfo[intervalEnd] + 1,
						//	[dim](const LightBVHLightInfo &a,
						//		const LightBVHLightInfo &b) {
						//	return a.bounds_o.axis[dim] <
						//		b.bounds_o.axis[dim];
						//});
						//std::nth_element(&lightInfo[intervalEnd + 1], &lightInfo[rightMid],
						//	&lightInfo[end - 1] + 1,
						//	[dim](const LightBVHLightInfo &a,
						//		const LightBVHLightInfo &b) {
						//	return a.bounds_o.axis[dim] <
						//		b.bounds_o.axis[dim];
						//});

						//// setting the median axis
						//Vector3f leftAxis = lightInfo[leftMid].bounds_o.axis;
						//Vector3f rightAxis = lightInfo[rightMid].bounds_o.axis;

						Vector3f leftAxis(0.f, 0.f, 0.f);
						Vector3f rightAxis(0.f, 0.f, 0.f);

						for (int j = start; j < intervalEnd + 1; j++) {
							leftAxis += lightInfo[j].bounds_o.axis;
						}

						for (int j = intervalEnd + 1; j < end; j++) {
							rightAxis += lightInfo[j].bounds_o.axis;
						}

						if (leftAxis == Vector3f(0.f, 0.f, 0.f)) {
							leftAxis = lightInfo[start].bounds_o.axis;
						}
						else {
							leftAxis = Normalize(leftAxis);
						}

						if (rightAxis == Vector3f(0.f, 0.f, 0.f)) {
							rightAxis = lightInfo[intervalEnd + 1].bounds_o.axis;
						}
						else {
							rightAxis = Normalize(rightAxis);
						}

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
						float leftAngle = 2 * Pi * ((1 - cos(leftO)) + ((2 * sin(leftO) * leftE + cos(leftO) - cos(2 * leftE + leftO)) / 4));
						float rightAngle = 2 * Pi * ((1 - cos(rightO)) + ((2 * sin(rightO) * rightE + cos(rightO) - cos(2 * rightE + rightO)) / 4));
						float leftCost = b0.SurfaceArea() * leftAngle * leftEnergy;
						float rightCost = b1.SurfaceArea() * rightAngle * rightEnergy;
						if (i == 0) {
							totalEnergy = leftEnergy + rightEnergy;
							Bounds3f totalBounds = Union(b0, b1);
							totalCost = totalBounds.SurfaceArea() * totalAngle * totalEnergy;
						}
						cost[i + (dim * (buckets - 1))] = (leftCost + rightCost) / totalCost;
					}
				}
			}
		}
		// Find bucket to split at that minimizes SAOH metric
		Float minCost = cost[0];
		int minCostSplitBucket = 0;
		for (int i = 1; i < (buckets - 1) * 3; i++) {
			if (cost[i] < minCost) {
				minCost = cost[i];
				minCostSplitBucket = i;
			}
		}
		// find out dimension and buckets accordingly
		int dim = minCostSplitBucket / (buckets - 1);
		minCostSplitBucket %= (buckets - 1);

		// split lights at selected SAOH bucket, mid is the last element of the first part
		mid = start + lightsPerBucket * (minCostSplitBucket + 1) - 1;
		std::nth_element(&lightInfo[(int)(start)], &lightInfo[mid],
			&lightInfo[end - 1] + 1,
			[dim](const LightBVHLightInfo &a,
				const LightBVHLightInfo &b) {
			return a.centroid[dim] <
				b.centroid[dim];
		});
		// recursively build children
		node->InitInterior(dim,
			recursiveBuild(arena, lightInfo, start, mid + 1,
				totalNodes),
			recursiveBuild(arena, lightInfo, mid + 1, end,
				totalNodes), bounds_o[dim], totalEnergy);
		return node;
	}

	int LightBVHAccel::flattenLightBVHTree(LightBVHNode *node, int *offset) {
		LinearLightBVHNode *linearNode = &nodes[*offset];
		linearNode->bounds_w = node->bounds_w;
		linearNode->bounds_o = node->bounds_o;
		linearNode->energy = node->energy;
		linearNode->centroid = node->centroid;
		linearNode->nLights = node->nLights;
		int myOffset = (*offset)++;
		if (node->nLights == 1) {
			linearNode->lightNum = node->lightNum;
		}
		else {
			// Create interior flattened BVH node
			linearNode->splitAxis = node->splitAxis;
			flattenLightBVHTree(node->children[0], offset);
			linearNode->secondChildOffset =
				flattenLightBVHTree(node->children[1], offset);
		}
		return myOffset;
	}

	void LightBVHAccel::calculateThetas(std::vector<LightBVHLightInfo> &lightInfo, int startIndex, int endIndex, Vector3f axis, float *theta_o, float *theta_e) {
		for (int j = startIndex; j < endIndex + 1; j++) {
			LightBVHLightInfo l = lightInfo[j];
			*theta_o = std::max(*theta_o, acos(Dot(axis, l.bounds_o.axis)) + l.bounds_o.theta_o);
			*theta_e = std::max(*theta_e, *theta_o + l.bounds_o.theta_e);
		}
		*theta_e -= *theta_o;
		*theta_e = std::max(0.f, std::min(*theta_e, Pi - *theta_o));
		*theta_o = std::min(Pi, *theta_o);
	}

	std::shared_ptr<LightBVHAccel> CreateLightBVHAccelerator(
		std::vector<std::shared_ptr<Light>> lights, float splitThreshold) {
		return std::make_shared<LightBVHAccel>(std::move(lights), splitThreshold);
	}

	int LightBVHAccel::SampleOneLight(const Interaction &it, float *pdf) {
		float sample1D = rng.UniformFloat();
		*pdf = 1;
		return TraverseNodeForOneLight(nodes, sample1D, it, pdf);
	}

	int LightBVHAccel::TraverseNodeForOneLight(LinearLightBVHNode *node, float sample1D, const Interaction &it, float *pdf) {
		// im already at a leaf of the tree
		if (node->nLights == 1) {
			return node->lightNum;
		}
		// finding out if I have to take the left or the right path
		float firstImportance = calculateImportance(it, &node[1]);
		float secondImportance = calculateImportance(it, &nodes[node->secondChildOffset]);
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
			return TraverseNodeForOneLight(&node[1], sample1D / firstImportance, it, pdf);
		}
		// right child traversal
		*pdf *= secondImportance;
		return TraverseNodeForOneLight(&nodes[node->secondChildOffset], (sample1D - firstImportance) / secondImportance, it, pdf);
	}

	std::vector<std::pair<int, float>> LightBVHAccel::SampleMultipleLights(const Interaction &it) {
		std::vector<std::pair<int, float>> lightVector;
		TraverseNodeForMultipleLights(nodes, it, &lightVector);
		return lightVector;
	}

	void LightBVHAccel::TraverseNodeForMultipleLights(LinearLightBVHNode *node, const Interaction &it, std::vector<std::pair<int, float>> *lightVector) {
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
			float r1 = rng.UniformFloat();
			float r2 = rng.UniformFloat();
			Point2f uScattering(r1, r2);
			BxDFType bsdfFlags = BSDF_ALL;
			const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
			Vector3f wi;
			float scatteringPdf;
			BxDFType sampledType;
			isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
				bsdfFlags, &sampledType);
			scatteringPdf /= Pi;
			//Normal3f n = it.n;
			//if (Dot(it.wo, n) < 0) {
			//	n *= -1;
			//}
			if (scatteringPdf != 0) {
				wi = Normalize(wi);
				Vector3f d = Normalize(node->centroid - o);
				float maxCos = acos(Dot(wi, d)) / Pi;
				float maxAngle = 0;
				Vector3f c[8];
				for (int i = 0; i < 8; i++) {
					c[i] = Normalize(b.Corner(i) - o);
				}
				for (int i = 0; i < 7; i++) {
					for (int j = i + 1; j < 8; j++) {
						maxAngle = std::max(maxAngle, acos(Dot(c[i], c[j])));
					}
				}
				maxAngle /= Pi;
				if (maxAngle * maxCos * scatteringPdf > splitThreshold) {
					split = true;
				}
				//}
			}
		}
		if (split) {
			TraverseNodeForMultipleLights(&node[1], it, lightVector);
			TraverseNodeForMultipleLights(&nodes[node->secondChildOffset], it, lightVector);
		}
		else {
			float pdf = 1.f;
			float sample1D = rng.UniformFloat();
			int lightNum = TraverseNodeForOneLight(node, sample1D, it, &pdf);
			lightVector->push_back(std::pair<int, float>(lightNum, pdf));
		}
	}

	float LightBVHAccel::calculateImportance(const Interaction &it, LinearLightBVHNode* node) {
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
				if (Dot(n, node->bounds_w.Corner(i) - o) >= 0) {
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

		float t0 = -std::numeric_limits<float>::infinity();
		for (int i = 0; i < 3; ++i) {
			float invRayDir = 1 / d[i];
			float tNear = std::min((node->bounds_w.pMin[i] - o[i]) * invRayDir, (node->bounds_w.pMax[i] - o[i]) * invRayDir);
			t0 = std::max(t0, tNear);
		}

		// shading point is already in the box -> can always find a theta_u with 1 --> angleImportance is always 1
		float angleImportance = 1;
		// shading point is not in the box
		if (t0 > 0) {
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
			if (theta - theta_o - theta_u - theta_e < 0) {
				angleImportance = cos(std::max(0.f, theta - theta_o - theta_u));
			}
			else {
				return 0;
			}
			//angleImportance = cos(std::max(0.f, std::min(theta - theta_o - theta_u, theta_e)));
		}
		return node->energy * angleImportance / (distance * distance);
	}
}