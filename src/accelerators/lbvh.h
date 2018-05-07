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
		bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
		bool IntersectP(const Ray &ray) const;

	private:
		// BVHAccel Private Methods
		LBVHBuildNode * recursiveBuild(
			MemoryArena &arena, std::vector<LBVHLightInfo> &LightInfo,
			int start, int end, int *totalNodes,
			std::vector<std::shared_ptr<Light>> &orderedLights);

		// BVHAccel Private Data
		const int maxLightsInNode;
		std::vector<std::shared_ptr<Light>> lights;
		LinearLBVHNode *nodes = nullptr;
	};

std::shared_ptr<LBVHAccel> CreateLBVHAccelerator(
	std::vector<std::shared_ptr<Light>> lights, const ParamSet &ps);

} // namespace pbrt

#endif  // PBRT_ACCELERATORS_LBVH_H