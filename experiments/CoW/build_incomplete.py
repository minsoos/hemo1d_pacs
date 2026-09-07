import os
import json
import build_network

COMMUNICATING_ARTERIES = {4, 11, 16}
VESSEL_TO_ELIMINATE = 4

def remove_communicating_artery(network, vessel_id):
    if vessel_id not in COMMUNICATING_ARTERIES:
        raise ValueError(
            f"Expected a communicating artery ID in "
            f"{sorted(COMMUNICATING_ARTERIES)}, got {vessel_id}"
        )

    vessel_name = next(
        vessel["name"]
        for vessel in network["vessels"]
        if vessel["id"] == vessel_id
    )

    # Remove the vessel itself.
    network["vessels"] = [
        vessel
        for vessel in network["vessels"]
        if vessel["id"] != vessel_id
    ]

    affected_nodes = []

    # Remove both vessel-end connections from their junctions.
    for node in network["nodes"]:
        old_connections = node["connections"]

        new_connections = [
            connection
            for connection in old_connections
            if connection["vessel"] != vessel_id
        ]

        if len(new_connections) == len(old_connections):
            continue
        if len(new_connections) < 2:
            raise RuntimeError(
                f"Removing vessel {vessel_id} leaves node "
                f"{node['id']} with fewer than two connections"
            )

        node["connections"] = new_connections

        node.pop("bifurcation_angles_rad", None)

        affected_nodes.append(node["name"])

    if len(affected_nodes) != 2:
        raise RuntimeError(
            f"Expected vessel {vessel_id} to connect two junctions, "
            f"but modified {len(affected_nodes)} nodes"
        )

    network["_description"] = (
        "Incomplete Circle of Willis with Windkessel terminal models; "
        f"missing vessel {vessel_id}: {vessel_name}"
    )

    return network


def main():
    GENERATED_DIR = build_network.GENERATED_DIR
    os.makedirs(GENERATED_DIR, exist_ok=True)
    

    network, n_cycles = build_network.build_network_json()
    network = remove_communicating_artery(network, vessel_id=VESSEL_TO_ELIMINATE)

    build_network.build_input_csv(n_cycles)

    network_path = os.path.join(GENERATED_DIR, "cow_network_windkessel_incomplete.json")
    with open(network_path, "w") as f:
        json.dump(network, f, indent=2)
    print(f"Wrote {network_path}")

    total_elements = sum(v["n_elements"] for v in network["vessels"])
    print(f"\n{len(network['vessels'])} vessels, {len(network['nodes'])} nodes, "
          f"{total_elements} elements total (target h={build_network.H} cm, p={build_network.POLYNOMIAL_ORDER}).")

    
if __name__ == "__main__":
    main()