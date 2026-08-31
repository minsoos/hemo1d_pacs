#pragma once

#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/section_state.hpp"
namespace hemo1d::physics {

class BloodFlowModel;

class NumericalFlux {
public:
    virtual ~NumericalFlux() = default;

    virtual SectionState compute(
        SectionState left, SectionState right, const VesselParameters& p,
        const BloodFlowModel& model
    ) const = 0;
};


class LaxFriedrichsFlux : public NumericalFlux {
public:
    SectionState compute(
        SectionState left, SectionState right, const VesselParameters& p,
        const BloodFlowModel& model
    ) const override; 
};


class HllFlux : public NumericalFlux {
public:
    SectionState compute(
        SectionState left, SectionState right, const VesselParameters& p,
        const BloodFlowModel& model
    ) const override;
};

} // namespace hemo1d::physics