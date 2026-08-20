#include "hemo1d/io/network_parser.hpp"

#include <fstream>
#include <string>
#include <nlohmann/json.hpp>
#include <stdexcept>

namespace hemo1d{
namespace{

using json = nlohmann::json;

VesselEnd parseVesselEnd(const std::string& s){
    if (s=="proximal") return VesselEnd::Proximal;
    if (s=="distal") return VesselEnd::Distal;
    throw std::runtime_error("Unknown Vessel type: \"" + s + "\" (expected \"proximal\" or \"distal\")");

}

BoundaryConditionType parseBcType(const std::string& s){
    if (s=="non_reflecting") return BoundaryConditionType::NonReflecting;
    if (s=="prescribed") return BoundaryConditionType::Prescribed;
    throw std::runtime_error("Unknown boundary condition type: \"" + s + "\" ( \"non_reflecting\" or \"prescribed\")");
}

PrescribedQuantity parsePrescribedQuantity(const std::string& s){
    if (s=="flow_rate") return PrescribedQuantity::FlowRate;
    if (s=="velocity") return PrescribedQuantity::Velocity;
    if (s=="area") return PrescribedQuantity::Area;
    if (s=="pressure") return PrescribedQuantity::Pressure;
    throw std::runtime_error("Unknown prescribed condition type: \"" + s + 
        "\" ( \"flow_rate\", \"velocity\", \"area\", \"pressure\")");
}

Vessel parseVessel(const json& j){
    const Id id = j.at("id").get<Id>();
    const std::string name = j.value("name", std::string());

    VesselParameters params;
    params.length = j.at("length").get<Real>();
    params.A0 = j.at("A0").get<Real>();
    params.beta = j.at("beta").get<Real>();
    params.alpha = j.value("alpha", 4.0/3.0);
    params.frictionKr = j.value("friction_kr", 0.0);
    params.nElements = j.value("n_elements", Index{1});

    if (j.contains("polynomial_order")){
        params.polynomialOrder = j.at("polynomial_order").get<unsigned>();
    }
    
    return Vessel(id, name, params);
}

std::optional<BoundaryConditionSpec> parseBoundaryCondition(const json& j,
        const std::filesystem::path& baseDir) {
            if (!j.contains("boundary_condition")) return std::nullopt;
            const json& bc = j.at("boundary_condition");

            BoundaryConditionSpec spec;
            spec.type = parseBcType(bc.at("type").get<std::string>());
            if (spec.type == BoundaryConditionType::Prescribed){
                spec.quantity = parsePrescribedQuantity(bc.at("quantity").get<std::string>());
                const std::string csvFile = bc.at("csv_file").get<std::string>();
                spec.csvFile = baseDir / csvFile;
            }
            return spec;
}

Node parseNode(const json& j, const std::filesystem::path& baseDir){
    const Id id = j.at("id").get<Id>();
    const std::string name = j.value("name", std::string());

    std::vector<VesselConnection> connections;
    for (const json& c : j.at("connections")){
        connections.push_back(
            VesselConnection{c.at("vessel").get<Id>(), parseVesselEnd(c.at("end").get<std::string>())}
        );
    }

    std::optional<BoundaryConditionSpec> bc = parseBoundaryCondition(j, baseDir);

    std::vector<Real> angles;
    if (j.contains("bifurcation_angles_rad")){
        angles = j.at("bifurcation_angles_rad").get<std::vector<Real>>();
    }

    return Node(id, name, std::move(connections), std::move(bc), std::move(angles));
}
} // namespace
} // namespace hemo1d

namespace hemo1d::io{
Network loadNetwork(const std::filesystem::path& jsonFile){
    std::ifstream file(jsonFile);
    if (!file){
        throw std::runtime_error("loadNetwork: Unknown path: \"" + jsonFile.string() + "\"" );
    }

    json root;

    try{
        file >> root;
    } catch (const json::parse_error& e){
        throw std::runtime_error("loadNetwork: Path: \"" + jsonFile.string() + "\" has a parse error: " + e.what());
    }

    const std::filesystem::path baseDir = jsonFile.parent_path();

    // Get fluid properties
    FluidProperties fluid;
    if (root.contains("fluid")){
        const json& f = root.at("fluid");
        fluid.density = f.value("density", fluid.density);
        fluid.viscosity = f.value("viscosity", fluid.viscosity);
    }

    // Parse Vessels
    std::vector<Vessel> vessels;
    const json& vesselsJson = root.at("vessels");
    vessels.reserve(vesselsJson.size());
    for (Index i=0; i<vesselsJson.size(); ++i){
        try{
            vessels.push_back(parseVessel(vesselsJson[i]));
        } catch (const std::exception& e) {
            throw std::runtime_error("loadNetwork: Error at parsing vessel n : " + std::to_string(i) + ", due to:" + e.what());
        }
    }

    // Parse Nodes
    std::vector<Node> nodes;
    const json& nodesJson = root.at("nodes");
    vessels.reserve(nodesJson.size());
    for (Index i=0; i<nodesJson.size(); ++i){
        try{
            nodes.push_back(parseNode(nodesJson[i], baseDir));
        } catch (const std::exception& e) {
            throw std::runtime_error("loadNetwork: Error at parsing node n : " + std::to_string(i) + ", due to:" + e.what());
        }
    }

    return Network(fluid, std::move(vessels), std::move(nodes));
}
} // namespace hemo1d::io