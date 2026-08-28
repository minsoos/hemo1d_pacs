#pragma once

#include <optional>
#include <vector>

#include "hemo1d/core/vessel.hpp"
#include "hemo1d/dg/reference_element.hpp"

namespace hemo1d::dg {

class Element {
public:
    Element(
        Index flatIndex, Id vesselId, Index localIndex, Real zLeft, Real zRight,
        const ReferenceElement& referenceElement, VesselParameters vesselParameters,
        std::optional<Index> leftNeighbor, std::optional<Index> rightNeighbor, 
        Index dofOffset
    );

    Index flatIndex() const noexcept { return flatIndex_; }
    Id vesselId() const noexcept { return vesselId_; }
    Index localIndex() const noexcept { return localIndex_; }

    Real zLeft() const noexcept { return zLeft_; }
    Real zRight() const noexcept { return zRight_; }
    Real length() const noexcept { return zRight_ - zLeft_; }

    Real jacobian() const noexcept { return 0.5 * length(); }

    const ReferenceElement& referenceElement() const noexcept { return *referenceElement_; }
    unsigned order() const noexcept { return referenceElement_->order(); }
    Index numNodes() const noexcept { return referenceElement_->numNodes(); }
    Index dofOffset() const noexcept { return dofOffset_; }

    const VesselParameters& vesselParameters() const noexcept { return vesselParameters_; }
    const std::vector<Real>& physicalNodes() const noexcept { return physicalNodes_; }

    std::optional<Index> leftNeighbor() const noexcept { return leftNeighbor_; }
    std::optional<Index> rightNeighbor() const noexcept { return rightNeighbor_; }

    bool isFirstInVessel( ) const noexcept { return !leftNeighbor_.has_value(); }
    bool isLastInVessel( ) const noexcept { return !rightNeighbor_.has_value(); }

private:
    Index flatIndex_;
    Id vesselId_;
    Index localIndex_;
    Real zLeft_;
    Real zRight_;
    const ReferenceElement* referenceElement_;
    VesselParameters vesselParameters_;
    std::vector<Real> physicalNodes_;
    std::optional<Index> leftNeighbor_;
    std::optional<Index> rightNeighbor_;
    Index dofOffset_;
};


}