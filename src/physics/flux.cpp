#include "hemo1d/physics/flux.hpp"

#include <algorithm>
#include <cmath>

#include "hemo1d/physics/blood_flow_model.hpp"

namespace hemo1d::physics {

SectionState LaxFriedrichsFlux::compute(
    SectionState left, SectionState right, const VesselParameters& p,
    const BloodFlowModel& model
) const {

    const SectionState Fl = model.physicalFlux(left, p);
    const SectionState Fr = model.physicalFlux(right, p);

    const Characteristics cl = model.characteristics(left, p);
    const Characteristics cr = model.characteristics(right, p);
    const Real lambdaMax = std::max({
        std::abs(cl.lambdaMinus), std::abs(cl.lambdaPlus),
        std::abs(cr.lambdaMinus), std::abs(cr.lambdaPlus)
    });

    return {
        0.5 * (Fl.A + Fr.A) - 0.5 * lambdaMax * (right.A - left.A),
        0.5 * (Fl.Q + Fr.Q) - 0.5 * lambdaMax * (right.Q - left.Q)
    };
}


SectionState HllFlux::compute(
    SectionState left, SectionState right, const VesselParameters& p,
    const BloodFlowModel& model
) const {

    const Characteristics cl = model.characteristics(left, p);
    const Characteristics cr = model.characteristics(right, p);

    const Real sL = std::min(cl.lambdaMinus, cr.lambdaMinus);
    const Real sR = std::max(cl.lambdaPlus, cr.lambdaPlus);

    const SectionState Fl = model.physicalFlux(left, p);
    if (sL >= 0.0) return Fl;

    const SectionState Fr = model.physicalFlux(right, p);
    if (sR <= 0.0) return Fr;

    const Real denom = sR - sL;

    return {
        (sR * Fl.A - sL * Fr.A + sR * sL * (right.A - left.A)) / denom,
        (sR * Fl.Q - sL * Fr.Q + sR * sL * (right.Q - left.Q)) / denom
    };

}

} // namespace hemo1d::physics