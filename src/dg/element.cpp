#include "hemo1d/dg/element.hpp"

#include <utility>

namespace hemo1d::dg {

Element::Element(
    Index flatIndex, Id vesselId, Index localIndex, Real zLeft, Real zRight,
    const ReferenceElement& referenceElement, VesselParameters vesselParameters,
    std::optional<Index> leftNeighbor, std::optional<Index> rightNeighbor,
    Index dofOffset
)   : flatIndex_(flatIndex),
      vesselId_(vesselId),
      localIndex_(localIndex),
      zLeft_(zLeft),
      zRight_(zRight),
      referenceElement_(&referenceElement),
      vesselParameters_(std::move(vesselParameters)),
      leftNeighbor_(leftNeighbor),
      rightNeighbor_(rightNeighbor),
      dofOffset_(dofOffset)
{
    const std::vector<Real>& refNodes = referenceElement_->nodes();
    physicalNodes_.resize(refNodes.size());
    for (Index i = 0; i < refNodes.size(); ++i) {
        physicalNodes_[i] = zLeft_ + 0.5 * (1.0 + refNodes[i]) * length();
    }   
}

} // namespace hemo1d::dg