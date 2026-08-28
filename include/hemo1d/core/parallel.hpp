#pragma once
 
#include "hemo1d/core/types.hpp"
 
namespace hemo1d {

// Below this many loop iterations (elements, or degrees of freedom,
// depending on the loop), OpenMP's thread-pool wake/join overhead exceeds
// whatever time is saved by splitting the work, so the loop runs on a
// single thread instead. Every `#pragma omp parallel for` in this codebase
// is guarded with `if (count >= kOmpParallelThreshold)`, so the dispatch
// decision is made once per call at runtime rather than duplicating a
// sequential and a parallel code path.
inline constexpr Index kOmpParallelThreshold = 256;

// Separate, higher threshold for loops that do very little work per
// iteration (MinmodLimiter's two element passes and Solver's two RK
// stage-combine loops). These pay the same fixed per-region OpenMP
// wake/join cost as SpatialOperator's loop, but have much less per-element
// work to amortize it against, so the break-even element/DOF count is much
// higher
inline constexpr Index kOmpParallelThresholdLight = 8192;
 
} // namespace hemo1d