import os
import sys
import hemo1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import common as c


def main():
    os.makedirs(c.GENERATED_DIR, exist_ok=True)
    pulse_csv_path = os.path.join(c.GENERATED_DIR, c.PULSE_CSV_NAME)
    c.write_pulse_csv(pulse_csv_path)
    print(f"Wrote inlet pulse CSV: {pulse_csv_path}")

    network_path = os.path.join(c.GENERATED_DIR, "p1_n0008.json")
    c.write_network_json(network_path, 8, 1, c.PULSE_CSV_NAME)
    print(f"Wrote {network_path}")

    network = hemo1d.load_network(network_path)
    print(f"Loaded: {network.vessel_count} vessels, {network.node_count} nodes")


if __name__ == "__main__":
    main()