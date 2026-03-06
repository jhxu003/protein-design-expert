"""
Shared PDB parsing and structural analysis utilities for protein binder design.

Provides lightweight PDB parsing (no BioPython dependency), secondary structure
assignment, surface analysis, and epitope classification. All other scripts
in this package depend on the data structures defined here.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
import math
import json
import sys


# ── Amino acid constants ─────────────────────────────────────────────────────

AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}

# Polar/charged residues (DEKRHNQST) used for hydrophilicity assessment
POLAR_RESIDUES = set("DEKRHNQST")

# Hydrophobic residues
HYDROPHOBIC_RESIDUES = set("AVILMFYW")

# Van der Waals radii (Å) for SASA approximation
VDW_RADII = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "H": 1.2}


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class Atom:
    """Single atom from a PDB ATOM record."""
    serial: int
    name: str
    resname: str
    chain: str
    resnum: int
    x: float
    y: float
    z: float
    element: str
    bfactor: float = 0.0

    @property
    def coord(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass
class Residue:
    """Single amino acid residue."""
    chain: str
    resnum: int
    resname: str  # 3-letter code
    ca_coord: Tuple[float, float, float]
    atoms: List[Atom] = field(default_factory=list, repr=False)
    secondary_structure: Optional[str] = None  # "H" (helix), "E" (strand), "C" (coil)
    sasa: Optional[float] = None
    bfactor: Optional[float] = None

    @property
    def one_letter(self) -> str:
        return AA_3TO1.get(self.resname, "X")

    @property
    def cb_coord(self) -> Tuple[float, float, float]:
        """Return CB coordinate, or CA for glycine."""
        for atom in self.atoms:
            if atom.name.strip() == "CB":
                return atom.coord
        return self.ca_coord


@dataclass
class ContactPair:
    """A pair of contacting residues across an interface."""
    res1_chain: str
    res1_resnum: int
    res2_chain: str
    res2_resnum: int
    min_distance: float


@dataclass
class Interface:
    """Interface between two protein chains."""
    residues_chain1: List[Residue]
    residues_chain2: List[Residue]
    area: float  # Å²  (estimated)
    contact_pairs: List[ContactPair] = field(default_factory=list)

    @property
    def n_contacts(self) -> int:
        return len(self.contact_pairs)


@dataclass
class EpitopeRegion:
    """Classified epitope region on a target protein."""
    residues: List[Residue]
    classification: str  # "helix", "loop", "flat", "concave", "edge_strand", "mixed"
    hydrophilicity_score: float  # 0-1, fraction polar/charged
    curvature: str  # "convex", "flat", "concave"
    secondary_structure_composition: Dict[str, float] = field(default_factory=dict)
    # {"H": 0.4, "E": 0.3, "C": 0.3}

    def to_dict(self) -> dict:
        return {
            "residue_range": f"{self.residues[0].chain}{self.residues[0].resnum}-{self.residues[-1].resnum}",
            "classification": self.classification,
            "hydrophilicity_score": round(self.hydrophilicity_score, 3),
            "curvature": self.curvature,
            "secondary_structure_composition": {
                k: round(v, 3) for k, v in self.secondary_structure_composition.items()
            },
            "n_residues": len(self.residues),
        }


# ── Geometry helpers ─────────────────────────────────────────────────────────

def distance(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def centroid(coords: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    """Centroid of a list of 3D points."""
    n = len(coords)
    if n == 0:
        return (0.0, 0.0, 0.0)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    return (cx, cy, cz)


def angle_between_vectors(v1: Tuple[float, float, float],
                          v2: Tuple[float, float, float]) -> float:
    """Angle in degrees between two 3D vectors."""
    dot = sum(a * b for a, b in zip(v1, v2))
    mag1 = math.sqrt(sum(a ** 2 for a in v1))
    mag2 = math.sqrt(sum(a ** 2 for a in v2))
    if mag1 < 1e-9 or mag2 < 1e-9:
        return 0.0
    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))


# ── PDB Parser ───────────────────────────────────────────────────────────────

class PDBParser:
    """
    Lightweight PDB parser. Works with raw ATOM records — no BioPython needed.

    Usage:
        parser = PDBParser()
        chains = parser.parse("target.pdb")
        epitope = parser.classify_epitope(chains["A"][10:30])
    """

    def parse(self, pdb_path: str) -> Dict[str, List[Residue]]:
        """
        Parse a PDB file into a dict of chain_id -> List[Residue].
        Only ATOM records for standard amino acids are included.
        """
        atoms_by_chain_res: Dict[str, Dict[int, List[Atom]]] = {}
        resname_map: Dict[str, Dict[int, str]] = {}

        with open(pdb_path, "r") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue

                resname = line[17:20].strip()
                if resname not in AA_3TO1:
                    continue

                chain = line[21].strip() or "A"
                try:
                    resnum = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except (ValueError, IndexError):
                    continue

                atom_name = line[12:16].strip()
                serial = int(line[6:11].strip()) if line[6:11].strip().isdigit() else 0
                bfactor = float(line[60:66].strip()) if len(line) >= 66 and line[60:66].strip() else 0.0
                element = line[76:78].strip() if len(line) >= 78 else atom_name[0]

                atom = Atom(
                    serial=serial, name=atom_name, resname=resname,
                    chain=chain, resnum=resnum,
                    x=x, y=y, z=z, element=element, bfactor=bfactor,
                )

                atoms_by_chain_res.setdefault(chain, {}).setdefault(resnum, []).append(atom)
                resname_map.setdefault(chain, {})[resnum] = resname

        # Build Residue objects
        chains: Dict[str, List[Residue]] = {}
        for chain_id in sorted(atoms_by_chain_res):
            residues = []
            for resnum in sorted(atoms_by_chain_res[chain_id]):
                atoms = atoms_by_chain_res[chain_id][resnum]
                resname = resname_map[chain_id][resnum]

                ca_coord = None
                for a in atoms:
                    if a.name == "CA":
                        ca_coord = a.coord
                        break
                if ca_coord is None:
                    continue

                avg_bfactor = sum(a.bfactor for a in atoms) / len(atoms)

                residues.append(Residue(
                    chain=chain_id, resnum=resnum, resname=resname,
                    ca_coord=ca_coord, atoms=atoms, bfactor=avg_bfactor,
                ))
            chains[chain_id] = residues

        return chains

    def assign_secondary_structure(self, residues: List[Residue]) -> None:
        """
        Simple secondary structure assignment based on CA-CA distances and
        dihedral-like geometry. Not as accurate as DSSP but works without
        external dependencies.

        Heuristic:
        - Alpha helix: consecutive CA-CA ~ 3.8Å, i to i+3 ~ 5.1-5.5Å
        - Beta strand: consecutive CA-CA ~ 3.8Å, i to i+1 extended (~6.5-7.5Å pitch)
        - Coil: everything else
        """
        n = len(residues)
        for i in range(n):
            residues[i].secondary_structure = "C"  # default coil

        for i in range(n - 3):
            d_i_i3 = distance(residues[i].ca_coord, residues[i + 3].ca_coord)
            # Alpha helix signature: i to i+3 distance ~5.0-5.5 Å
            if 4.8 <= d_i_i3 <= 5.8:
                for j in range(i, min(i + 4, n)):
                    if residues[j].secondary_structure == "C":
                        residues[j].secondary_structure = "H"

        for i in range(n - 1):
            d_ca = distance(residues[i].ca_coord, residues[i + 1].ca_coord)
            if residues[i].secondary_structure == "C" and d_ca > 3.5:
                # Check for extended conformation
                if i + 2 < n:
                    d_i_i2 = distance(residues[i].ca_coord, residues[i + 2].ca_coord)
                    if d_i_i2 > 6.0:
                        residues[i].secondary_structure = "E"
                        residues[i + 1].secondary_structure = "E"

    def compute_interface(self, chain1_residues: List[Residue],
                          chain2_residues: List[Residue],
                          distance_cutoff: float = 8.0) -> Interface:
        """
        Compute the interface between two chains.

        A residue is at the interface if any of its heavy atoms are within
        `distance_cutoff` Å of any heavy atom on the other chain.
        Returns interface residues, contact pairs, and estimated area.
        """
        interface_res1 = set()
        interface_res2 = set()
        contact_pairs = []

        for r1 in chain1_residues:
            for r2 in chain2_residues:
                # Quick CB distance pre-filter
                cb_dist = distance(r1.cb_coord, r2.cb_coord)
                if cb_dist > distance_cutoff + 5.0:
                    continue

                # Check atom-level contacts
                min_dist = float("inf")
                for a1 in r1.atoms:
                    if a1.element == "H":
                        continue
                    for a2 in r2.atoms:
                        if a2.element == "H":
                            continue
                        d = distance(a1.coord, a2.coord)
                        min_dist = min(min_dist, d)
                        if d < distance_cutoff:
                            break
                    if min_dist < distance_cutoff:
                        break

                if min_dist < distance_cutoff:
                    interface_res1.add(r1.resnum)
                    interface_res2.add(r2.resnum)
                    contact_pairs.append(ContactPair(
                        res1_chain=r1.chain, res1_resnum=r1.resnum,
                        res2_chain=r2.chain, res2_resnum=r2.resnum,
                        min_distance=round(min_dist, 2),
                    ))

        iface_r1 = [r for r in chain1_residues if r.resnum in interface_res1]
        iface_r2 = [r for r in chain2_residues if r.resnum in interface_res2]

        # Rough area estimate: ~120 Å² per interface residue (empirical average)
        estimated_area = (len(iface_r1) + len(iface_r2)) * 60.0

        return Interface(
            residues_chain1=iface_r1,
            residues_chain2=iface_r2,
            area=estimated_area,
            contact_pairs=contact_pairs,
        )

    def classify_epitope(self, residues: List[Residue]) -> EpitopeRegion:
        """
        Classify a set of residues as an epitope region.

        Returns the epitope type (helix, loop, flat, concave, edge_strand, mixed),
        hydrophilicity score, curvature estimate, and secondary structure composition.
        """
        if not residues:
            return EpitopeRegion(
                residues=[], classification="unknown",
                hydrophilicity_score=0.0, curvature="flat",
                secondary_structure_composition={},
            )

        # Ensure secondary structure is assigned
        if residues[0].secondary_structure is None:
            self.assign_secondary_structure(residues)

        # Secondary structure composition
        ss_counts = {"H": 0, "E": 0, "C": 0}
        for r in residues:
            ss = r.secondary_structure or "C"
            ss_counts[ss] = ss_counts.get(ss, 0) + 1
        total = len(residues)
        ss_fractions = {k: v / total for k, v in ss_counts.items()}

        # Hydrophilicity
        hydro_score = hydrophilicity_fraction(residues)

        # Curvature estimate from CA coordinates
        curvature = self._estimate_curvature(residues)

        # Classification logic
        classification = self._classify_ss(ss_fractions, curvature, hydro_score)

        return EpitopeRegion(
            residues=residues,
            classification=classification,
            hydrophilicity_score=hydro_score,
            curvature=curvature,
            secondary_structure_composition=ss_fractions,
        )

    def _classify_ss(self, ss_frac: Dict[str, float], curvature: str,
                     hydro: float) -> str:
        """Determine epitope classification from secondary structure and geometry."""
        h_frac = ss_frac.get("H", 0)
        e_frac = ss_frac.get("E", 0)
        c_frac = ss_frac.get("C", 0)

        # Dominant helix (>50%)
        if h_frac >= 0.50:
            return "helix"

        # Edge strand: beta-dominant + hydrophilic
        if e_frac >= 0.40 and hydro >= 0.50:
            return "edge_strand"

        # Beta sheet (not edge strand = less hydrophilic)
        if e_frac >= 0.40:
            if curvature == "concave":
                return "concave"
            return "flat"

        # Loop-dominant
        if c_frac >= 0.60:
            return "loop"

        # Concave pocket
        if curvature == "concave":
            return "concave"

        # Mixed
        return "mixed"

    def _estimate_curvature(self, residues: List[Residue]) -> str:
        """
        Estimate local surface curvature from CA coordinates.
        Uses deviation from the plane fitted through the CA atoms.
        """
        if len(residues) < 5:
            return "flat"

        coords = [r.ca_coord for r in residues]
        center = centroid(coords)

        # Compute deviations from centroid
        deviations = [distance(c, center) for c in coords]
        mean_dev = sum(deviations) / len(deviations)
        max_dev = max(deviations)

        # Simple heuristic: if max deviation is much larger than mean,
        # the surface curves away (convex). If all are similar, it's flat.
        # For concave, check if points curve inward (edges farther than center).
        edge_devs = [deviations[0], deviations[-1]]
        mid_devs = deviations[len(deviations) // 3 : 2 * len(deviations) // 3]

        if mid_devs:
            avg_mid = sum(mid_devs) / len(mid_devs)
            avg_edge = sum(edge_devs) / len(edge_devs) if edge_devs else avg_mid

            if avg_edge > avg_mid * 1.3:
                return "concave"
            elif avg_mid > avg_edge * 1.3:
                return "convex"

        return "flat"


# ── Utility functions ─────────────────────────────────────────────────────────

def hydrophilicity_fraction(residues: List[Residue]) -> float:
    """Fraction of polar/charged residues (DEKRHNQST) in the given residues."""
    if not residues:
        return 0.0
    polar_count = sum(1 for r in residues if r.one_letter in POLAR_RESIDUES)
    return polar_count / len(residues)


def get_residues_by_numbers(chain_residues: List[Residue],
                            resnums: List[int]) -> List[Residue]:
    """Select residues by their residue numbers."""
    resnum_set = set(resnums)
    return [r for r in chain_residues if r.resnum in resnum_set]


def get_residue_range(chain_residues: List[Residue],
                      start: int, end: int) -> List[Residue]:
    """Select residues within a range [start, end] inclusive."""
    return [r for r in chain_residues if start <= r.resnum <= end]


def sequence_from_residues(residues: List[Residue]) -> str:
    """Get one-letter amino acid sequence from a list of residues."""
    return "".join(r.one_letter for r in residues)


def summarize_structure(chains: Dict[str, List[Residue]]) -> dict:
    """Generate a quick summary of a parsed PDB structure."""
    summary = {}
    for chain_id, residues in chains.items():
        summary[chain_id] = {
            "n_residues": len(residues),
            "sequence_length": len(residues),
            "first_resnum": residues[0].resnum if residues else None,
            "last_resnum": residues[-1].resnum if residues else None,
        }
    return summary


# ── CLI entry point ──────────────────────────────────────────────────────────

def main():
    """Quick PDB analysis from command line."""
    import argparse
    parser = argparse.ArgumentParser(description="Parse and analyze PDB structure")
    parser.add_argument("pdb", help="PDB file path")
    parser.add_argument("--chain", help="Analyze specific chain")
    parser.add_argument("--epitope", nargs="+", type=int,
                        help="Epitope residue numbers to classify")
    parser.add_argument("--interface", nargs=2, metavar=("CHAIN1", "CHAIN2"),
                        help="Compute interface between two chains")
    args = parser.parse_args()

    pdb_parser = PDBParser()
    chains = pdb_parser.parse(args.pdb)

    print(f"Parsed {args.pdb}:")
    for cid, residues in chains.items():
        pdb_parser.assign_secondary_structure(residues)
        ss_counts = {"H": 0, "E": 0, "C": 0}
        for r in residues:
            ss_counts[r.secondary_structure or "C"] += 1
        print(f"  Chain {cid}: {len(residues)} residues "
              f"(H:{ss_counts['H']} E:{ss_counts['E']} C:{ss_counts['C']})")

    if args.epitope and args.chain:
        if args.chain in chains:
            epitope_res = get_residues_by_numbers(chains[args.chain], args.epitope)
            if epitope_res:
                epitope = pdb_parser.classify_epitope(epitope_res)
                print(f"\nEpitope classification:")
                print(json.dumps(epitope.to_dict(), indent=2))

    if args.interface:
        c1, c2 = args.interface
        if c1 in chains and c2 in chains:
            iface = pdb_parser.compute_interface(chains[c1], chains[c2])
            print(f"\nInterface {c1}-{c2}:")
            print(f"  Residues: {len(iface.residues_chain1)} ({c1}) + "
                  f"{len(iface.residues_chain2)} ({c2})")
            print(f"  Contact pairs: {iface.n_contacts}")
            print(f"  Estimated area: {iface.area:.0f} Å²")


if __name__ == "__main__":
    main()
