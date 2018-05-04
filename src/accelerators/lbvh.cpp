// accelerators/lbvh.cpp*
#include "accelerators/lbvh.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>
#include "core/light.h"

namespace pbrt {

	STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
	STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
	STAT_COUNTER("BVH/Interior nodes", interiorNodes);
	STAT_COUNTER("BVH/Leaf nodes", leafNodes);

	struct Bounds_o {
		Vector3f axis;
		float theta_o;
		float theta_e;

		Bounds_o(Vector3f axis, float theta_o, float theta_e)
			: axis(axis),
			theta_o(theta_o),
			theta_e(theta_e) {}
	};

	// BVHAccel Local Declarations
	struct LBVHPrimitiveInfo {

		float energy;
		size_t primitiveNumber;
		Bounds3f bounds_w;
		Bounds_o bounds_o;
		Point3f centroid;

		LBVHPrimitiveInfo(float energy, size_t primitiveNumber, const Bounds3f &bounds_w, Light &light)
			: energy(energy),
			primitiveNumber(primitiveNumber),
			bounds_w(bounds_w),
			centroid(.5f * bounds_w.pMin + .5f * bounds_w.pMax),
			bounds_o(Bounds_o(light.Axis, light.Theta_o, light.Theta_e)) {}
	};

	struct LBVHBuildNode {
		// LBVHBuildNode Public Methods
		void InitLeaf(int first, int n, const Bounds3f &b) {
			firstPrimOffset = first;
			nPrimitives = n;
			bounds = b;
			children[0] = children[1] = nullptr;
			++leafNodes;
			++totalLeafNodes;
			totalPrimitives += n;
		}
		void InitInterior(int axis, LBVHBuildNode *c0, LBVHBuildNode *c1) {
			children[0] = c0;
			children[1] = c1;
			bounds = Union(c0->bounds, c1->bounds);
			splitAxis = axis;
			nPrimitives = 0;
			++interiorNodes;
		}
		Bounds3f bounds;
		LBVHBuildNode *children[2];
		int splitAxis, firstPrimOffset, nPrimitives;
	};
}