#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_LBVH_H
#define PBRT_ACCELERATORS_LBVH_H

// accelerators/lbvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <atomic>

namespace pbrt {
	struct Bounds_o;
	struct LBVHBuildNode;

	// LBVHAccel Forward Declarations
	struct LBVHLightInfo;
	struct MortonPrimitive;
	struct LinearLBVHNode;

	class LBVHAccel : public Aggregate {
	public:

		// BVHAccel Public Methods
		LBVHAccel(std::vector<std::shared_ptr<Light>> l,
			int maxPrimsInNode = 1);
		Bounds3f WorldBound() const;
		~LBVHAccel();

	private:
		// BVHAccel Private Methods
		LBVHBuildNode * recursiveBuild(
			MemoryArena &arena, std::vector<LBVHLightInfo> &LightInfo,
			int start, int end, int *totalNodes,
			std::vector<std::shared_ptr<Light>> &orderedLights);
		LBVHBuildNode *HLBVHBuild(
			MemoryArena &arena, const std::vector<LBVHLightInfo> &LightInfo,
			int *totalNodes,
			std::vector<std::shared_ptr<Light>> &orderedLights) const;
		LBVHBuildNode *emitLBVH(
			LBVHBuildNode *&buildNodes,
			const std::vector<LBVHLightInfo> &primitiveInfo,
			MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
			std::vector<std::shared_ptr<Light>> &orderedLights,
			std::atomic<int> *orderedPrimsOffset, int bitIndex) const;
		LBVHBuildNode *buildUpperSAH(MemoryArena &arena,
			std::vector<LBVHBuildNode *> &treeletRoots,
			int start, int end, int *totalNodes) const;
		int flattenLBVHTree(LBVHBuildNode *node, int *offset);

		// BVHAccel Private Data
		const int maxLightsInNode;
		std::vector<std::shared_ptr<Light>> lights;
		LinearLBVHNode *nodes = nullptr;
	};
}

#endif  // PBRT_ACCELERATORS_LBVH_H