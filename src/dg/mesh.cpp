#include "hemo1d/dg/mesh.hpp"

#include <stdexcept>
#include <string>

namespace hemo1d::dg {

Mesh::Mesh(
    const Network& network, const DgSettings& settings
) {
    for (const Vessel& vessel : network.vessels()) {
        const unsigned order = 
            vessel.parameters().polynomialOrder.value_or(settings.defaultPolynomialOrder);
        referenceElements_.try_emplace(order, order);
    }

    for (const Vessel& vessel : network.vessels()) {
        const VesselParameters& params = vessel.parameters();
        if (params.nElements == 0) {
            throw std::runtime_error("Mesh: vessel " + std::to_string(vessel.id()) + " has zero elements.");
        }

        const unsigned order = params.polynomialOrder.value_or(settings.defaultPolynomialOrder);
        const ReferenceElement& refElem = referenceElements_.at(order);

        const Index rangeBegin = elements_.size();
        const Real h = params.length / static_cast<Real>(params.nElements);

        for (Index local = 0; local < params.nElements; ++local) {
            const Real zLeft = static_cast<Real>(local) * h;
            const Real zRight = zLeft + h;
            
            const std::optional<Index> leftNeighbor = 
                local == 0 ? std::nullopt : std::optional<Index>(elements_.size() - 1);
            
            const std::optional<Index> rightNeighbor =
                local + 1 == params.nElements ? std::nullopt : std::optional<Index>(elements_.size() + 1);
            
            elements_.emplace_back(
                elements_.size(), vessel.id(), local, zLeft, zRight,
                refElem, params, leftNeighbor, rightNeighbor, totalDofs_
            );
            totalDofs_ += refElem.numNodes();
        }

        vesselElementRange_.emplace(vessel.id(), std::make_pair(rangeBegin, elements_.size()));    
    }
}


std::pair<Index, Index> Mesh::vesselElementRange(Id vesselId) const {
    const auto it = vesselElementRange_.find(vesselId);
    if (it == vesselElementRange_.end()) {
        throw std::runtime_error("Mesh::vesselElementRange: unknown vessel id " + std::to_string(vesselId));
    }
    return it->second;
}

const ReferenceElement& Mesh::referenceElement(unsigned order) const {
    const auto it = referenceElements_.find(order);
    if (it == referenceElements_.end()) {
        throw std::runtime_error("Mesh::referenceElement: no reference element cached for order " + std::to_string(order));
    }
    return it->second;
}

} // namespace hemo1d::dg