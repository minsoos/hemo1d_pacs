#pragma once
 
#include <cstddef>
#include <limits>
 
namespace hemo1d {
 
using Real = double;
using Index = std::size_t;
 
using Id = std::size_t;
 
inline constexpr Id kInvalidId = std::numeric_limits<Id>::max();
 
 
} // namespace hemo1d
 