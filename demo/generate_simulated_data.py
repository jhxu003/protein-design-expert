"""Generate realistic simulated design results for EGFR demo."""
import random
import csv
import os

random.seed(42)

designs = []

def gen(design_id, profile):
    """Generate a design with noise around a profile."""
    return {
        "design_id": design_id,
        "ipsae_min": round(max(0, profile["ipsae"] + random.gauss(0, 0.05)), 3),
        "iptm": round(min(1, max(0, profile["iptm"] + random.gauss(0, 0.05))), 3),
        "plddt": round(min(1, max(0, profile["plddt"] + random.gauss(0, 0.04))), 3),
        "pae": round(max(0, profile["pae"] + random.gauss(0, 2.0)), 1),
        "interface_area": round(max(0, profile["area"] + random.gauss(0, 100)), 0),
        "ca_rmsd": round(max(0, profile["rmsd"] + random.gauss(0, 0.3)), 2),
        "hotspot_contact_rate": round(min(1, max(0, profile["hotspot"] + random.gauss(0, 0.05))), 3),
        "shape_complementarity": round(min(1, max(0, profile["sc"] + random.gauss(0, 0.05))), 3),
        "binder_length": profile["length"] + random.randint(-5, 5),
    }

# Simulate 500 designs with realistic distribution for EGFR binder campaign
# ~3% gold, ~8% silver, ~15% bronze, rest fail (typical for RFdiffusion)
profiles = {
    "gold": {"ipsae": 0.12, "iptm": 0.90, "plddt": 0.94, "pae": 6.5, "area": 1350, "rmsd": 0.8, "hotspot": 0.38, "sc": 0.77, "length": 75},
    "silver": {"ipsae": 0.22, "iptm": 0.82, "plddt": 0.88, "pae": 9.5, "area": 1150, "rmsd": 1.2, "hotspot": 0.30, "sc": 0.65, "length": 90},
    "bronze": {"ipsae": 0.35, "iptm": 0.72, "plddt": 0.83, "pae": 13.0, "area": 950, "rmsd": 1.7, "hotspot": 0.23, "sc": 0.55, "length": 105},
    "type_i_severe": {"ipsae": 0.58, "iptm": 0.42, "plddt": 0.52, "pae": 23.0, "area": 620, "rmsd": 3.8, "hotspot": 0.10, "sc": 0.28, "length": 80},
    "type_i_partial": {"ipsae": 0.45, "iptm": 0.55, "plddt": 0.75, "pae": 17.5, "area": 750, "rmsd": 2.8, "hotspot": 0.15, "sc": 0.40, "length": 88},
    "type_ii_no_bind": {"ipsae": 0.50, "iptm": 0.35, "plddt": 0.89, "pae": 20.0, "area": 1050, "rmsd": 1.0, "hotspot": 0.28, "sc": 0.50, "length": 78},
    "type_ii_small": {"ipsae": 0.38, "iptm": 0.68, "plddt": 0.86, "pae": 12.5, "area": 680, "rmsd": 1.3, "hotspot": 0.22, "sc": 0.52, "length": 58},
    "poor_hotspot": {"ipsae": 0.30, "iptm": 0.76, "plddt": 0.87, "pae": 10.0, "area": 1100, "rmsd": 1.0, "hotspot": 0.10, "sc": 0.62, "length": 85},
    "weak_binding": {"ipsae": 0.32, "iptm": 0.68, "plddt": 0.84, "pae": 13.0, "area": 900, "rmsd": 1.4, "hotspot": 0.25, "sc": 0.48, "length": 100},
}

# Distribution: realistic for EGFR with RFdiffusion
counts = {
    "gold": 15, "silver": 40, "bronze": 75,
    "type_i_severe": 60, "type_i_partial": 95,
    "type_ii_no_bind": 80, "type_ii_small": 45,
    "poor_hotspot": 40, "weak_binding": 50,
}

idx = 0
for profile_name, count in counts.items():
    for j in range(count):
        idx += 1
        d = gen(f"egfr_{idx:04d}", profiles[profile_name])
        designs.append(d)

# Shuffle
random.shuffle(designs)

# Write CSV
output = os.path.join(os.path.dirname(os.path.abspath(__file__)), "egfr_round1_results.csv")
with open(output, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(designs[0].keys()))
    writer.writeheader()
    writer.writerows(designs)

print(f"Generated {len(designs)} simulated EGFR binder designs → {output}")
