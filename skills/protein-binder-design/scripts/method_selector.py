"""
Protein binder design method selection engine.

Encodes expert knowledge about five computational binder design methods,
their strengths/weaknesses, and when to use each. Provides automated
method ranking based on target characteristics.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict
import argparse
import json


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class MethodProfile:
    """Profile of a computational binder design method."""
    name: str
    success_rate: str
    speed: str
    needs_downstream: List[str]  # e.g. ["ProteinMPNN", "AF2"]
    best_for: List[str]
    strengths: List[str]
    weaknesses: List[str]
    reference: str = ""

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "success_rate": self.success_rate,
            "speed": self.speed,
            "needs_downstream": self.needs_downstream,
            "best_for": self.best_for,
            "strengths": self.strengths,
            "weaknesses": self.weaknesses,
        }


@dataclass
class MethodRecommendation:
    """A ranked method recommendation with rationale."""
    method: MethodProfile
    rank: int
    score: float  # 0-100 composite suitability score
    rationale: str
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "rank": self.rank,
            "method": self.method.name,
            "score": round(self.score, 1),
            "rationale": self.rationale,
            "warnings": self.warnings,
        }


# ── Method database ──────────────────────────────────────────────────────────

METHODS: Dict[str, MethodProfile] = {}

METHODS["RFdiffusion"] = MethodProfile(
    name="RFdiffusion",
    success_rate="0.4-7.5% baseline",
    speed="slow",
    needs_downstream=["ProteinMPNN", "AF2"],
    best_for=[
        "general targets",
        "helical epitopes",
        "motif scaffolding",
    ],
    strengths=[
        "modular architecture",
        "strong community support",
        "flexible hotspot specification",
    ],
    weaknesses=[
        "lower success on hydrophilic targets",
        "requires multi-step pipeline",
    ],
    reference="Watson et al., Nature 2023",
)

METHODS["BindCraft"] = MethodProfile(
    name="BindCraft",
    success_rate="10-100%",
    speed="fast (one-shot)",
    needs_downstream=[],
    best_for=[
        "rapid prototyping",
        "loop epitopes",
        "beginners",
    ],
    strengths=[
        "one-shot design (no separate sequence design step)",
        "integrated AF2 validation",
        "repredicts complex at each step",
    ],
    weaknesses=[
        "newer method with less community data",
    ],
    reference="Pacesa et al., 2024",
)

METHODS["Latent-X"] = MethodProfile(
    name="Latent-X",
    success_rate="90%+ macrocycles, 10-64% mini-binders",
    speed="10X faster than pipelines",
    needs_downstream=[],
    best_for=[
        "macrocyclic binders",
        "high-affinity requirements",
        "speed-critical projects",
    ],
    strengths=[
        "picomolar affinities achievable",
        "all-atom generation",
        "very fast inference",
    ],
    weaknesses=[
        "best for macrocycles specifically",
    ],
    reference="Chu et al., 2024",
)

METHODS["PXDesign"] = MethodProfile(
    name="PXDesign",
    success_rate="17-82%",
    speed="medium",
    needs_downstream=["AF2"],
    best_for=[
        "diverse targets",
        "mixed epitope types",
    ],
    strengths=[
        "dual modality (diffusion + hallucination)",
        "broad fold coverage",
    ],
    weaknesses=[
        "variable success rate",
        "needs parameter tuning",
    ],
    reference="Wang et al., 2024",
)

METHODS["beta-pairing-RFdiffusion"] = MethodProfile(
    name="beta-pairing-RFdiffusion",
    success_rate="9.2% (vs 0.98% standard RFdiffusion)",
    speed="medium",
    needs_downstream=["ProteinMPNN", "AF2"],
    best_for=[
        "edge-strand epitopes",
        "hydrophilic targets",
    ],
    strengths=[
        "~10X improvement for edge-strand targets",
        "specifically trained for beta-pairing interfaces",
    ],
    weaknesses=[
        "narrow applicability (edge-strand only)",
    ],
    reference="Goudy et al., 2024",
)


# ── Scoring and selection ────────────────────────────────────────────────────

def _score_method(
    method: MethodProfile,
    epitope_type: str,
    hydrophilicity: str,
    time_budget: str,
    affinity_requirement: str,
) -> float:
    """
    Score a method 0-100 for a given set of target characteristics.

    Parameters
    ----------
    epitope_type : str
        One of: helix, loop, flat, concave, edge_strand, mixed, unknown
    hydrophilicity : str
        One of: high, medium, low
    time_budget : str
        One of: fast, standard, unlimited
    affinity_requirement : str
        One of: standard, high, picomolar
    """
    score = 50.0  # baseline

    name = method.name

    # ── Epitope type matching ────────────────────────────────────────────
    if epitope_type == "edge_strand":
        if name == "beta-pairing-RFdiffusion":
            score += 35
        elif name == "PXDesign":
            score += 10
        elif name == "RFdiffusion":
            score -= 10  # standard RFdiffusion struggles here
        elif name == "BindCraft":
            score += 5
    elif epitope_type == "loop":
        if name == "BindCraft":
            score += 30
        elif name == "PXDesign":
            score += 10
        elif name == "RFdiffusion":
            score += 5
    elif epitope_type == "helix":
        if name == "RFdiffusion":
            score += 25
        elif name == "BindCraft":
            score += 20
        elif name == "PXDesign":
            score += 10
    elif epitope_type == "flat":
        if name == "RFdiffusion":
            score += 20
        elif name == "PXDesign":
            score += 20
        elif name == "BindCraft":
            score += 10
    elif epitope_type == "concave":
        if name == "RFdiffusion":
            score += 15
        elif name == "PXDesign":
            score += 20
        elif name == "BindCraft":
            score += 10
    elif epitope_type == "mixed":
        if name == "PXDesign":
            score += 25
        elif name == "BindCraft":
            score += 15
        elif name == "RFdiffusion":
            score += 10

    # ── Hydrophilicity ───────────────────────────────────────────────────
    if hydrophilicity == "high":
        if name == "beta-pairing-RFdiffusion":
            score += 15
        elif name == "RFdiffusion":
            score -= 15  # known weakness
        elif name == "BindCraft":
            score += 5
    elif hydrophilicity == "low":
        if name == "RFdiffusion":
            score += 10  # works well on hydrophobic
        elif name == "beta-pairing-RFdiffusion":
            score -= 10  # not its strength

    # ── Time budget ──────────────────────────────────────────────────────
    if time_budget == "fast":
        if name == "BindCraft":
            score += 20
        elif name == "Latent-X":
            score += 25
        elif name == "RFdiffusion":
            score -= 15  # multi-step pipeline is slow
        elif name == "beta-pairing-RFdiffusion":
            score -= 10
    elif time_budget == "unlimited":
        # All methods viable; slight bonus for thoroughness
        if name == "RFdiffusion":
            score += 5
        elif name == "PXDesign":
            score += 5

    # ── Affinity requirement ─────────────────────────────────────────────
    if affinity_requirement == "picomolar":
        if name == "Latent-X":
            score += 30
        elif name == "BindCraft":
            score += 5
        elif name == "RFdiffusion":
            score -= 5
    elif affinity_requirement == "high":
        if name == "Latent-X":
            score += 15
        elif name == "BindCraft":
            score += 5

    return max(0.0, min(100.0, score))


def select_method(
    epitope_type: str = "mixed",
    hydrophilicity: str = "medium",
    time_budget: str = "standard",
    affinity_requirement: str = "standard",
) -> List[MethodRecommendation]:
    """
    Rank binder design methods for a given target scenario.

    Parameters
    ----------
    epitope_type : str
        One of: helix, loop, flat, concave, edge_strand, mixed, unknown
    hydrophilicity : str
        One of: high, medium, low
    time_budget : str
        One of: fast, standard, unlimited
    affinity_requirement : str
        One of: standard, high, picomolar

    Returns
    -------
    List[MethodRecommendation]
        Methods ranked from best to worst with scores and rationale.
    """
    scored: List[tuple] = []
    for method in METHODS.values():
        score = _score_method(method, epitope_type, hydrophilicity,
                              time_budget, affinity_requirement)
        scored.append((method, score))

    # Sort by score descending
    scored.sort(key=lambda x: x[1], reverse=True)

    recommendations = []
    for rank, (method, score) in enumerate(scored, 1):
        rationale = _build_rationale(method, epitope_type, hydrophilicity,
                                     time_budget, affinity_requirement)
        warnings = _build_warnings(method, epitope_type, hydrophilicity)

        recommendations.append(MethodRecommendation(
            method=method,
            rank=rank,
            score=score,
            rationale=rationale,
            warnings=warnings,
        ))

    return recommendations


def _build_rationale(
    method: MethodProfile,
    epitope_type: str,
    hydrophilicity: str,
    time_budget: str,
    affinity_requirement: str,
) -> str:
    """Generate a human-readable rationale for why a method was ranked."""
    reasons = []
    name = method.name

    # Epitope fit
    epitope_fit = {
        ("RFdiffusion", "helix"): "RFdiffusion excels at helical epitopes",
        ("RFdiffusion", "flat"): "RFdiffusion handles flat surfaces well",
        ("RFdiffusion", "concave"): "RFdiffusion can scaffold into concave pockets",
        ("BindCraft", "loop"): "BindCraft is strongest on loop epitopes",
        ("BindCraft", "helix"): "BindCraft handles helical epitopes well",
        ("PXDesign", "mixed"): "PXDesign dual modality suits mixed epitopes",
        ("PXDesign", "flat"): "PXDesign covers broad surface types",
        ("PXDesign", "concave"): "PXDesign hallucination mode works for concave targets",
        ("beta-pairing-RFdiffusion", "edge_strand"): "specifically designed for edge-strand beta-pairing",
        ("Latent-X", "helix"): "Latent-X generates compact binders for helical targets",
    }
    key = (name, epitope_type)
    if key in epitope_fit:
        reasons.append(epitope_fit[key])

    # Hydrophilicity
    if hydrophilicity == "high" and name == "beta-pairing-RFdiffusion":
        reasons.append("strong on hydrophilic targets")
    elif hydrophilicity == "high" and name == "RFdiffusion":
        reasons.append("may struggle with hydrophilic surfaces")

    # Speed
    if time_budget == "fast":
        if name in ("BindCraft", "Latent-X"):
            reasons.append("fast turnaround matches time constraint")
        elif name == "RFdiffusion":
            reasons.append("multi-step pipeline may be too slow")

    # Affinity
    if affinity_requirement in ("high", "picomolar") and name == "Latent-X":
        reasons.append("can achieve picomolar affinities")

    if not reasons:
        reasons.append("general-purpose method with moderate suitability")

    return "; ".join(reasons) + "."


def _build_warnings(
    method: MethodProfile,
    epitope_type: str,
    hydrophilicity: str,
) -> List[str]:
    """Generate warnings about potential issues with a method choice."""
    warnings = []

    if method.name == "RFdiffusion" and hydrophilicity == "high":
        warnings.append(
            "RFdiffusion has lower success rates on hydrophilic targets; "
            "consider beta-pairing-RFdiffusion or BindCraft instead"
        )

    if method.name == "beta-pairing-RFdiffusion" and epitope_type != "edge_strand":
        warnings.append(
            "beta-pairing-RFdiffusion is specialized for edge-strand epitopes; "
            "may not provide benefits for other epitope types"
        )

    if method.name == "Latent-X" and epitope_type not in ("helix", "mixed", "unknown"):
        warnings.append(
            "Latent-X is best validated for macrocyclic binders; "
            "performance on this epitope type is less characterized"
        )

    if method.needs_downstream:
        warnings.append(
            f"requires downstream tools: {', '.join(method.needs_downstream)}"
        )

    return warnings


# ── Comparison table ─────────────────────────────────────────────────────────

def print_comparison_table() -> str:
    """
    Generate a markdown comparison table of all methods.

    Returns the table as a string and also prints it.
    """
    header = (
        "| Method | Success Rate | Speed | Downstream | Best For |"
    )
    separator = (
        "|--------|-------------|-------|------------|----------|"
    )

    rows = []
    for method in METHODS.values():
        downstream = ", ".join(method.needs_downstream) if method.needs_downstream else "None (self-contained)"
        best_for = ", ".join(method.best_for)
        row = f"| {method.name} | {method.success_rate} | {method.speed} | {downstream} | {best_for} |"
        rows.append(row)

    table = "\n".join([header, separator] + rows)
    print(table)
    return table


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Select optimal protein binder design method based on target characteristics"
    )
    parser.add_argument(
        "--epitope-type",
        choices=["helix", "loop", "flat", "concave", "edge_strand", "mixed", "unknown"],
        default="mixed",
        help="Epitope structural type (default: mixed)",
    )
    parser.add_argument(
        "--hydrophilicity",
        choices=["high", "medium", "low"],
        default="medium",
        help="Epitope hydrophilicity level (default: medium)",
    )
    parser.add_argument(
        "--time-budget",
        choices=["fast", "standard", "unlimited"],
        default="standard",
        help="Available compute time (default: standard)",
    )
    parser.add_argument(
        "--affinity",
        choices=["standard", "high", "picomolar"],
        default="standard",
        help="Required binding affinity (default: standard)",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Print method comparison table and exit",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output as JSON",
    )

    args = parser.parse_args()

    if args.compare:
        print_comparison_table()
        return

    recommendations = select_method(
        epitope_type=args.epitope_type,
        hydrophilicity=args.hydrophilicity,
        time_budget=args.time_budget,
        affinity_requirement=args.affinity,
    )

    if args.json:
        output = {
            "query": {
                "epitope_type": args.epitope_type,
                "hydrophilicity": args.hydrophilicity,
                "time_budget": args.time_budget,
                "affinity_requirement": args.affinity,
            },
            "rankings": [r.to_dict() for r in recommendations],
        }
        print(json.dumps(output, indent=2))
    else:
        print(f"Method ranking for: epitope={args.epitope_type}, "
              f"hydrophilicity={args.hydrophilicity}, "
              f"time={args.time_budget}, affinity={args.affinity}")
        print("=" * 72)
        for rec in recommendations:
            marker = ">>>" if rec.rank == 1 else "   "
            print(f"{marker} #{rec.rank}  {rec.method.name:<30s}  score={rec.score:.1f}")
            print(f"       Success rate: {rec.method.success_rate}")
            print(f"       Rationale: {rec.rationale}")
            if rec.warnings:
                for w in rec.warnings:
                    print(f"       WARNING: {w}")
            print()


if __name__ == "__main__":
    main()
