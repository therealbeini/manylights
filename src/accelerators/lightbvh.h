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
	struct LightBVHNode;

	// LightBVHAccel Forward Declarations
	struct LightBVHLightInfo;

	class LightBVHAccel {
	public:

		LightBVHNode * root = nullptr;
		std::vector<std::shared_ptr<Light>> lights;
		float splitThreshold;

		// BVHAccel Public Methods
		LightBVHAccel(std::vector<std::shared_ptr<Light>> l, float splitThreshold);
		int SampleOneLight(const Interaction &it, Sampler &sampler, float *pdf);
		std::vector<std::pair<int, float>> SampleMultipleLights(const Interaction &it, Sampler &sampler);

	private:
		// BVHAccel Private Methods
		LightBVHNode* recursiveBuild(
			MemoryArena &arena, std::vector<LightBVHLightInfo> &LightInfo,
			int start, int end, int *totalNodes,
			std::vector<std::shared_ptr<Light>> &orderedLights);
		void calculateThetas(std::vector<LightBVHLightInfo> &lightInfo, int startIndex, int endIndex, Vector3f axis, float *theta_o, float *theta_e);
		int TraverseNodeForOneLight(LightBVHNode *node, float sample1D, const Interaction &it, float *pdf);
		void TraverseNodeForMultipleLights(LightBVHNode *node, float sample1D, const Interaction &it, std::vector<std::pair<int, float>> *lightVector);
		float calculateImportance(const Interaction & it, LightBVHNode * node);

		// BVHAccel Private Data
	};

	std::shared_ptr<LightBVHAccel> CreateLightBVHAccelerator(
		std::vector<std::shared_ptr<Light>> lights, float splitThreshold);

} // namespace pbrt

#endif  // PBRT_ACCELERATORS_LIGHTBVH_H