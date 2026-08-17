#pragma once
 
#include <cstddef>
#include <limits>
 
namespace hemo1d {
 
using Real = double;

// For containers positions
using Index = std::size_t;
 
// For (not necessarily contiguous) id coming from outside the code 
using Id = std::size_t;
 
inline constexpr Id kInvalidId = std::numeric_limits<Id>::max();
 
 
} // namespace hemo1d
 