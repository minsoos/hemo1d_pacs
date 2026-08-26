#pragma once

#include <optional>
#include <string>

#include "hemo1d/core/types.hpp"

namespace hemo1d {
    // Which extreme of the vessel a node is
    enum class VesselEnd {Proximal, Distal};

    // Discretization parameters:

    struct VesselParameters { 
        Real length = 0.0;
        Real A0 = 0.0;
        Real beta = 0.0;
        Real alpha = 4.0/3.0;
        Real frictionKr = 0.0;

        Index nElements = 1;
        std::optional<unsigned> polynomialOrder;
    };

    class Vessel {
        private:
            Id id_;
            std::string name_;
            VesselParameters params_;
        
        public:
            Vessel(Id id, std::string name, VesselParameters params);
 
            Id id() const noexcept { return id_;}
            const std::string& name() const noexcept { return name_;}
            const VesselParameters& parameters() const noexcept {return params_;}
    };

} //namespace hemo1d