#pragma once

#include <utility>

#include "hemo1d/core/types.hpp"
#include "hemo1d/core/vessel.hpp"

namespace hemo1d::physics {

class BoundaryStateProvider {
public:
    virtual ~BoundaryStateProvider() = default;
    virtual std::pair<Real, Real> exteriorState(
        Id vesselId, VesselEnd end, Real time
    ) const = 0;
};

} // namespace hemo1d::physics