
#include <iostream>
 
#include "hemo1d/core/types.hpp"
#include "hemo1d/io/network_parser.hpp"

namespace{

const char* endName(hemo1d::VesselEnd end) {
    return end == hemo1d::VesselEnd::Proximal ? "proximal" : "distal";
}

const char* nodeKindName(hemo1d::NodeKind kind) {
    return kind == hemo1d::NodeKind::Terminal ? "terminal" : "junction";
}

const char* boundaryConditionName(hemo1d::BoundaryConditionType typeBoundary) {
    return typeBoundary == hemo1d::BoundaryConditionType::Prescribed ? "prescribed" : "non_reflecting";
}

} // namespace

int main(int argc, char** argv){
    if (argc != 2){
        std::cerr << "usage: " << argv[0] << "<network.json>\n";
        return 1;
    }
    try{
        const hemo1d::Network network = hemo1d::io::loadNetwork(argv[1]);

        std::cout << "Fluid: Density=" << network.fluid().density << 
        " g/cm^3, viscosity" << network.fluid().density << " g/(cm s)\n\n";

        std::cout << "Vessels (" << network.vesselCount() << "):\n";
        for (const hemo1d::Vessel& v : network.vessels()){
            const auto& p = v.parameters();
            std::cout << " [" << v.id() << "] " << v.name() << ": L=" << p.length <<
                        " cm, A0=" <<p.A0 << " cm^2," << p.beta << " g/s^2, alpha="<<
                        p.alpha << ", n_elements=" << p.nElements << " -> proximal node " <<
                        network.nodeIdAt(v.id(), hemo1d::VesselEnd::Proximal) << ", distal node " <<
                        network.nodeIdAt(v.id(), hemo1d::VesselEnd::Distal) << "\n";
        }

        std::cout << "\n Nodes (" << network.nodeCount() << "):\n";
        for (const hemo1d::Node& n : network.nodes()){
            std::cout << " [" << n.id() << "] " << n.name() << " ()" << nodeKindName(n.kind()) <<
                        ") , connections:";
                        for (const auto& c : n.connections()){
                            std::cout << " vessel " << c.vesselId << " (" << endName(c.end) << ")";
                        }
                        if (n.boundaryCondition().has_value()){
                            std::cout << ", boundary_condition=" << boundaryConditionName(n.boundaryCondition()->type);
                        }
                        std::cout << "\n";
        }
        std::size_t totalDofs = 0, totalElements = 0;
        for (const auto& v : network.vessels()) {
            const auto& p = v.parameters();
            totalElements += p.nElements;
            totalDofs += p.nElements * (p.polynomialOrder.value_or(2) + 1);
        }
        std::cout << "\nTotals: " << totalElements << " elements, " << totalDofs << " DOFs\n";
        } catch (const std::exception& e){
            std::cerr << " error " << e.what() << "\n";
            return 1;
        }

        return 0;
    }