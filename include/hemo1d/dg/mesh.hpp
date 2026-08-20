#pragma once 

#include <unordered_map>
#include <utility>
#include <vector>

#include "hemo1d/core/network.hpp"
#include "hemo1d/dg/element.hpp"
#include "hemo1d/dg/reference_element.hpp"

namespace hemo1d::dg {

struct DgSettings {
    unsigned defaultPolynomialOrder = 1;
};


class Mesh {
public:
    explicit Mesh(const Network& network, const DgSettings& settings = {});

    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;
    Mesh(Mesh&&) = default;
    Mesh& operator=(Mesh&&) = default;

    const std::vector<Element>& elements() const noexcept { return elements_; }
    Index elementCount() const noexcept { return elements_.size(); }

    std::pair<Index, Index> vesselElementRange(Id vesselId) const;

    const ReferenceElement& referenceElement(unsigned order) const;

    Index totalDofs() const noexcept { return totalDofs_; }

private:
    std::unordered_map<unsigned, ReferenceElement> referenceElements_;
    std::vector<Element> elements_;
    std::unordered_map<Id, std::pair<Index, Index>> vesselElementRange_;
    Index totalDofs_ = 0;
};

} // namespace hemo1d::dg