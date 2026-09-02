// Small CLI utility: loads a network definition file and prints a summary.
// Usage: hemo1d_network_info <network.json>

#include <iomanip>
#include <iostream>

#include "hemo1d/io/network_parser.hpp"

namespace {

const char* endName(hemo1d::VesselEnd end) {
    return end == hemo1d::VesselEnd::Proximal ? "proximal" : "distal";
}

const char* nodeKindName(hemo1d::NodeKind kind) {
    return kind == hemo1d::NodeKind::Terminal ? "terminal" : "junction";
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " <network.json>\n";
        return 1;
    }

    try {
        const hemo1d::Network network = hemo1d::io::loadNetwork(argv[1]);

        std::cout << "Fluid: density=" << network.fluid().density
                  << " g/cm^3, viscosity=" << network.fluid().viscosity << " g/(cm s)\n\n";

        std::cout << "Vessels (" << network.vesselCount() << "):\n";
        for (const hemo1d::Vessel& v : network.vessels()) {
            const auto& p = v.parameters();
            std::cout << "  [" << v.id() << "] " << v.name() << ": L=" << p.length
                      << " cm, A0=" << p.A0 << " cm^2, beta=" << p.beta << " g/s^2, alpha=" << p.alpha
                      << ", n_elements=" << p.nElements << " -> proximal node "
                      << network.nodeIdAt(v.id(), hemo1d::VesselEnd::Proximal) << ", distal node "
                      << network.nodeIdAt(v.id(), hemo1d::VesselEnd::Distal) << '\n';
        }

        std::cout << "\nNodes (" << network.nodeCount() << "):\n";
        for (const hemo1d::Node& n : network.nodes()) {
            std::cout << "  [" << n.id() << "] " << n.name() << " (" << nodeKindName(n.kind())
                      << "), connections:";
            for (const auto& c : n.connections()) {
                std::cout << " vessel " << c.vesselId << " (" << endName(c.end) << ")";
            }
            if (n.boundaryCondition().has_value()) {
                std::cout << ", boundary_condition=";
                switch (n.boundaryCondition()->type) {
                    case hemo1d::BoundaryConditionType::Prescribed:
                        std::cout << "prescribed";
                        break;
                    case hemo1d::BoundaryConditionType::External:
                        std::cout << "external:" << n.boundaryCondition()->modelName;
                        break;
                    case hemo1d::BoundaryConditionType::NonReflecting:
                        std::cout << "non_reflecting";
                        break;
                }
            }
            std::cout << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
