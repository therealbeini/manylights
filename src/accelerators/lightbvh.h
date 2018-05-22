#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_LBVH_H
#define PBRT_ACCELERATORS_LBVH_H

// accelerators/lightbvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <atomic>

namespace pbrt {
	struct Bounds_o;
	struct LightBVHBuildNode;

	// LightBVHAccel Forward Declarations
	struct LightBVHLightInfo;
	struct MortonPrimitive;
	struct LinearLBVHNode;

	class LightBVHAccel {
	public:

		// BVHAccel Public Methods
		LightBVHAccel(std::vector<std::shared_ptr<Light>> l,
			int maxPrimsInNode = 1);

	private:
		// BVHAccel Private Methods
		LightBVHBuildNode * recursiveBuild(
			MemoryArena &arena, std::vector<LightBVHLightInfo> &LightInfo,
			int start, int end, int *totalNodes,
			std::vector<std::shared_ptr<Light>> &orderedLights);

		// BVHAccel Private Data
		const int maxLightsInNode;
		std::vector<std::shared_ptr<Light>> lights;
		LinearLBVHNode *nodes = nullptr;
	};

std::shared_ptr<LightBVHAccel> CreateLBVHAccelerator(
	std::vector<std::shared_ptr<Light>> lights, const ParamSet &ps);

} // namespace pbrt

#endif  // PBRT_ACCELERATORS_LBVH_H