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

	// BVHAccel Forward Declarations
	struct LBVHPrimitiveInfo;
}

#endif  // PBRT_ACCELERATORS_LBVH_H