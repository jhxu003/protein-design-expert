"""
Phase 1 target analysis engine for protein binder design.

Analyzes a target PDB structure, classifies the epitope, recommends hotspot
residues, ranks design methods, and generates warnings about potential issues
(flexibility, glycosylation, disorder). Produces both structured (JSON) and
human-readable reports.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
import argparse
import json
import re
import sys

from pdb_utils import (
    PDBParser,
    EpitopeRegion,
    Residue,
    get_residues_by_numbers,
    hydrophilicity_fraction,
    distance,
    AA_3TO1,
)
from method_selector import select_method, MethodRecommendation


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class HotspotRecommendation:
    """Recommended hotspot residues for binder design."""
    residues: List[int]
    rationale: str
    count: int
    distribution_quality: str  # "good", "clustered", "sparse"

    def to_dict(self) -> dict:
        return {
            "residues": self.residues,
            "rationale": self.rationale,
            "count": self.count,
            "distribution_quality": self.distribution_quality,
        }


@dataclass
class TargetAnalysis:
    """Complete Phase 1 analysis of a binder design target."""
    pdb_id: str
    target_chain: str
    epitope: EpitopeRegion
    hotspot_recommendation: HotspotRecommendation
    method_rankings: List[MethodRecommendation]
    warnings: List[str]
    summary: str

    def to_dict(self) -> dict:
        return {
            "pdb_id": self.pdb_id,
            "target_chain": self.target_chain,
            "epitope": self.epitope.to_dict(),
            "hotspot_recommendation": self.hotspot_recommendation.to_dict(),
            "method_rankings": [r.to_dict() for r in self.method_rankings],
            "warnings": self.warnings,
            "summary": self.summary,
        }

    def to_report(self) -> str:
        """Generate a human-readable analysis report."""
        lines = []
        lines.append("=" * 72)
        lines.append("PROTEIN BINDER DESIGN — TARGET ANALYSIS REPORT")
        lines.append("=" * 72)
        lines.append("")

        # Target info
        lines.append(f"Target:  {self.pdb_id}, chain {self.target_chain}")
        epi = self.epitope
        if epi.residues:
            resnums = [r.resnum for r in epi.residues]
            lines.append(f"Epitope: residues {min(resnums)}-{max(resnums)} "
                         f"({len(resnums)} residues)")
        lines.append(f"Type:    {epi.classification}")
        lines.append(f"Hydrophilicity: {epi.hydrophilicity_score:.2f} "
                      f"({'hydrophilic' if epi.hydrophilicity_score > 0.5 else 'hydrophobic' if epi.hydrophilicity_score < 0.3 else 'moderate'})")
        lines.append(f"Curvature: {epi.curvature}")

        ss = epi.secondary_structure_composition
        if ss:
            ss_str = ", ".join(f"{k}={v:.0%}" for k, v in sorted(ss.items()))
            lines.append(f"Secondary structure: {ss_str}")

        # Hotspots
        lines.append("")
        lines.append("-" * 40)
        lines.append("HOTSPOT RECOMMENDATION")
        lines.append("-" * 40)
        hs = self.hotspot_recommendation
        lines.append(f"Residues: {', '.join(str(r) for r in hs.residues)}")
        lines.append(f"Count: {hs.count}")
        lines.append(f"Distribution: {hs.distribution_quality}")
        lines.append(f"Rationale: {hs.rationale}")

        # Method ranking
        lines.append("")
        lines.append("-" * 40)
        lines.append("METHOD RANKING")
        lines.append("-" * 40)
        for rec in self.method_rankings:
            marker = ">>>" if rec.rank == 1 else "   "
            lines.append(f"{marker} #{rec.rank}  {rec.method.name:<30s}  "
                         f"score={rec.score:.1f}")
            lines.append(f"       {rec.rationale}")
            for w in rec.warnings:
                lines.append(f"       WARNING: {w}")

        # Warnings
        if self.warnings:
            lines.append("")
            lines.append("-" * 40)
            lines.append("TARGET WARNINGS")
            lines.append("-" * 40)
            for w in self.warnings:
                lines.append(f"  * {w}")

        # Summary
        lines.append("")
        lines.append("-" * 40)
        lines.append("SUMMARY")
        lines.append("-" * 40)
        lines.append(self.summary)
        lines.append("")

        return "\n".join(lines)


# ── Hotspot selection ────────────────────────────────────────────────────────

def _select_hotspots_helix(residues: List[Residue]) -> HotspotRecommendation:
    """
    Helix epitope: 3-5 hotspots, every 3-4 residues along the helix face.
    Selects residues that face outward (approximated by picking every
    3rd-4th residue, matching the helix turn period of 3.6 residues).
    """
    if len(residues) < 3:
        resnums = [r.resnum for r in residues]
        return HotspotRecommendation(
            residues=resnums, rationale="helix too short; using all residues",
            count=len(resnums), distribution_quality="clustered",
        )

    # Pick every 3-4 residues to stay on the same helix face
    selected = []
    i = 0
    while i < len(residues) and len(selected) < 5:
        selected.append(residues[i].resnum)
        i += 3  # ~one helix turn (3.6 residues)

    # Ensure at least 3
    if len(selected) < 3 and len(residues) >= 3:
        all_nums = [r.resnum for r in residues]
        for rn in all_nums:
            if rn not in selected:
                selected.append(rn)
            if len(selected) >= 3:
                break

    selected.sort()
    quality = _assess_distribution(selected, residues)

    return HotspotRecommendation(
        residues=selected,
        rationale="every 3-4 residues along helix face for consistent binding surface",
        count=len(selected),
        distribution_quality=quality,
    )


def _select_hotspots_loop(residues: List[Residue]) -> HotspotRecommendation:
    """
    Loop epitope: 2-4 hotspots, flanking residues, avoid flexible tip.
    The tip of a loop (middle residues) is typically most flexible;
    we prefer residues near the anchoring ends.
    """
    if len(residues) < 2:
        resnums = [r.resnum for r in residues]
        return HotspotRecommendation(
            residues=resnums, rationale="loop too short; using all residues",
            count=len(resnums), distribution_quality="clustered",
        )

    n = len(residues)
    selected = []

    # Flanking residues (near the anchored ends)
    selected.append(residues[0].resnum)
    if n > 1:
        selected.append(residues[-1].resnum)

    # Add one or two residues from the inner region, but not the very tip
    if n >= 5:
        # Pick residues at ~25% and ~75% along the loop
        q1 = max(1, n // 4)
        q3 = min(n - 2, 3 * n // 4)
        selected.append(residues[q1].resnum)
        if q3 != q1 and len(selected) < 4:
            selected.append(residues[q3].resnum)
    elif n >= 3:
        # Avoid the exact middle
        idx = 1 if n == 3 else 1
        selected.append(residues[idx].resnum)

    selected = sorted(set(selected))
    quality = _assess_distribution(selected, residues)

    return HotspotRecommendation(
        residues=selected,
        rationale="flanking residues preferred; avoiding flexible loop tip",
        count=len(selected),
        distribution_quality=quality,
    )


def _select_hotspots_flat(residues: List[Residue]) -> HotspotRecommendation:
    """
    Flat epitope: 4-6 hotspots, distributed across the surface.
    Spread them evenly by spacing along the residue list.
    """
    target_count = min(6, max(4, len(residues)))
    if len(residues) <= target_count:
        resnums = [r.resnum for r in residues]
        return HotspotRecommendation(
            residues=resnums,
            rationale="flat surface; using all available residues for coverage",
            count=len(resnums),
            distribution_quality="good" if len(resnums) >= 4 else "sparse",
        )

    step = max(1, len(residues) // target_count)
    selected = []
    for i in range(0, len(residues), step):
        selected.append(residues[i].resnum)
        if len(selected) >= target_count:
            break

    # Ensure we include the last residue for coverage
    if residues[-1].resnum not in selected and len(selected) < target_count:
        selected.append(residues[-1].resnum)

    selected.sort()
    quality = _assess_distribution(selected, residues)

    return HotspotRecommendation(
        residues=selected,
        rationale="evenly distributed across flat surface for broad contact",
        count=len(selected),
        distribution_quality=quality,
    )


def _select_hotspots_edge_strand(residues: List[Residue]) -> HotspotRecommendation:
    """
    Edge strand epitope: 3-5 hotspots, alternating i, i+2 pattern.
    Beta strands have side chains alternating above/below the sheet;
    picking i, i+2 ensures hotspots are on the same face.
    """
    if len(residues) < 3:
        resnums = [r.resnum for r in residues]
        return HotspotRecommendation(
            residues=resnums,
            rationale="edge strand too short; using all residues",
            count=len(resnums),
            distribution_quality="clustered",
        )

    # Alternating i, i+2 pattern
    selected = []
    i = 0
    while i < len(residues) and len(selected) < 5:
        selected.append(residues[i].resnum)
        i += 2

    selected.sort()
    quality = _assess_distribution(selected, residues)

    return HotspotRecommendation(
        residues=selected,
        rationale="alternating i, i+2 pattern to match beta-strand side-chain orientation",
        count=len(selected),
        distribution_quality=quality,
    )


def _select_hotspots_concave(residues: List[Residue]) -> HotspotRecommendation:
    """
    Concave epitope: 3-4 hotspots, rim + pocket residues.
    Pick residues at the edges (rim) and one or two in the center (pocket).
    """
    if len(residues) < 3:
        resnums = [r.resnum for r in residues]
        return HotspotRecommendation(
            residues=resnums,
            rationale="concave region too small; using all residues",
            count=len(resnums),
            distribution_quality="clustered",
        )

    n = len(residues)
    selected = []

    # Rim residues (edges)
    selected.append(residues[0].resnum)
    selected.append(residues[-1].resnum)

    # Pocket residues (middle)
    mid = n // 2
    selected.append(residues[mid].resnum)
    if n > 6:
        selected.append(residues[mid + 1].resnum if mid + 1 < n else residues[mid - 1].resnum)

    selected = sorted(set(selected))
    # Trim to 4 max
    selected = selected[:4]

    quality = _assess_distribution(selected, residues)

    return HotspotRecommendation(
        residues=selected,
        rationale="rim + pocket residues for concave surface engagement",
        count=len(selected),
        distribution_quality=quality,
    )


def _select_hotspots_mixed(residues: List[Residue]) -> HotspotRecommendation:
    """Fallback for mixed/unknown epitope types: use flat surface strategy."""
    return _select_hotspots_flat(residues)


_HOTSPOT_SELECTORS = {
    "helix": _select_hotspots_helix,
    "loop": _select_hotspots_loop,
    "flat": _select_hotspots_flat,
    "edge_strand": _select_hotspots_edge_strand,
    "concave": _select_hotspots_concave,
    "mixed": _select_hotspots_mixed,
    "unknown": _select_hotspots_mixed,
}


def _assess_distribution(selected: List[int], residues: List[Residue]) -> str:
    """Assess whether hotspot distribution is good, clustered, or sparse."""
    if len(selected) < 2:
        return "clustered"

    # Get the span of the epitope
    all_resnums = [r.resnum for r in residues]
    epitope_span = max(all_resnums) - min(all_resnums) + 1

    # Compute gaps between consecutive hotspots
    sorted_sel = sorted(selected)
    gaps = [sorted_sel[i + 1] - sorted_sel[i] for i in range(len(sorted_sel) - 1)]

    if not gaps:
        return "clustered"

    avg_gap = sum(gaps) / len(gaps)
    max_gap = max(gaps)
    min_gap = min(gaps)

    # Ideal gap: epitope_span / (n_hotspots - 1)
    ideal_gap = epitope_span / max(len(selected) - 1, 1)

    # Check clustering: if all hotspots are bunched together
    hotspot_span = sorted_sel[-1] - sorted_sel[0]
    if epitope_span > 0 and hotspot_span < epitope_span * 0.3:
        return "clustered"

    # Check sparsity: if there are very large gaps
    if max_gap > ideal_gap * 2.5:
        return "sparse"

    return "good"


# ── Warning detection ────────────────────────────────────────────────────────

def _detect_warnings(residues: List[Residue], chain_residues: List[Residue]) -> List[str]:
    """
    Detect potential issues with the target epitope.

    Checks for:
    - Flexible loops (high B-factor)
    - Glycosylation sites (N-X-S/T sequon)
    - Disordered regions (very high B-factor)
    """
    warnings = []

    if not residues:
        return warnings

    # ── Flexible loops (B-factor > 60) ───────────────────────────────────
    high_bfactor_residues = []
    for r in residues:
        if r.bfactor is not None and r.bfactor > 60.0:
            high_bfactor_residues.append(r.resnum)
    if high_bfactor_residues:
        warnings.append(
            f"High B-factor (>60) at residues {', '.join(str(x) for x in high_bfactor_residues)}: "
            f"indicates flexibility. Binder contacts here may be unreliable."
        )

    # ── Very disordered regions (B-factor > 80) ─────────────────────────
    disordered = [r.resnum for r in residues if r.bfactor is not None and r.bfactor > 80.0]
    if disordered:
        warnings.append(
            f"Possible disorder at residues {', '.join(str(x) for x in disordered)} "
            f"(B-factor >80). Consider excluding from epitope."
        )

    # ── Glycosylation sites (N-X-S/T sequon) ────────────────────────────
    # Build a lookup for the full chain to check sequons
    chain_by_resnum = {r.resnum: r for r in chain_residues}
    epitope_resnums = set(r.resnum for r in residues)

    for r in residues:
        if r.one_letter != "N":
            continue
        # Check for N-X-S/T sequon (X != P)
        next1 = chain_by_resnum.get(r.resnum + 1)
        next2 = chain_by_resnum.get(r.resnum + 2)
        if next1 and next2:
            if next1.one_letter != "P" and next2.one_letter in ("S", "T"):
                warnings.append(
                    f"Potential N-glycosylation site at N{r.resnum} "
                    f"(sequon N-{next1.one_letter}-{next2.one_letter}). "
                    f"Glycan may occlude binder binding."
                )

    # ── Average B-factor check ───────────────────────────────────────────
    bfactors = [r.bfactor for r in residues if r.bfactor is not None]
    if bfactors:
        avg_bf = sum(bfactors) / len(bfactors)
        if avg_bf > 50.0:
            warnings.append(
                f"Overall epitope has elevated B-factors (avg={avg_bf:.1f}). "
                f"Consider a more rigid region for higher binder success."
            )

    return warnings


# ── Main analyzer ────────────────────────────────────────────────────────────

class TargetAnalyzer:
    """
    Phase 1 target analysis for protein binder design.

    Analyzes a target PDB, classifies the epitope, selects hotspot residues,
    ranks design methods, and generates a comprehensive report.

    Usage:
        analyzer = TargetAnalyzer()
        result = analyzer.analyze("target.pdb", "A", epitope_residues=[100, 105, 110])
        print(result.to_report())
    """

    def __init__(self):
        self.parser = PDBParser()

    def analyze(
        self,
        pdb_path: str,
        target_chain: str,
        epitope_residues: Optional[List[int]] = None,
        time_constraint: str = "standard",
    ) -> TargetAnalysis:
        """
        Run full Phase 1 target analysis.

        Parameters
        ----------
        pdb_path : str
            Path to the target PDB file.
        target_chain : str
            Chain ID of the target protein.
        epitope_residues : list of int, optional
            Specific residue numbers defining the epitope. If None, the full
            chain surface is analyzed.
        time_constraint : str
            One of: "fast", "standard", "unlimited". Passed to method selector.

        Returns
        -------
        TargetAnalysis
            Complete analysis including epitope classification, hotspot
            recommendation, method rankings, and warnings.
        """
        # Parse structure
        chains = self.parser.parse(pdb_path)

        if target_chain not in chains:
            available = ", ".join(sorted(chains.keys()))
            raise ValueError(
                f"Chain '{target_chain}' not found in {pdb_path}. "
                f"Available chains: {available}"
            )

        chain_residues = chains[target_chain]

        # Assign secondary structure for the full chain
        self.parser.assign_secondary_structure(chain_residues)

        # Get epitope residues
        if epitope_residues:
            epi_res = get_residues_by_numbers(chain_residues, epitope_residues)
            if not epi_res:
                raise ValueError(
                    f"None of the specified epitope residues "
                    f"{epitope_residues} found in chain {target_chain}"
                )
        else:
            # Use all surface residues (simplified: use all residues)
            epi_res = chain_residues

        # Classify epitope
        epitope = self.parser.classify_epitope(epi_res)

        # Select hotspots
        selector_fn = _HOTSPOT_SELECTORS.get(
            epitope.classification, _select_hotspots_mixed
        )
        hotspot_rec = selector_fn(epi_res)

        # Map hydrophilicity to method selector categories
        hydro = epitope.hydrophilicity_score
        if hydro > 0.5:
            hydro_cat = "high"
        elif hydro < 0.3:
            hydro_cat = "low"
        else:
            hydro_cat = "medium"

        # Rank methods
        method_rankings = select_method(
            epitope_type=epitope.classification,
            hydrophilicity=hydro_cat,
            time_budget=time_constraint,
            affinity_requirement="standard",
        )

        # Detect warnings
        warnings = _detect_warnings(epi_res, chain_residues)

        # Derive PDB ID from filename
        pdb_id = pdb_path.rsplit("/", 1)[-1].rsplit("\\", 1)[-1]
        pdb_id = pdb_id.replace(".pdb", "").replace(".PDB", "")

        # Build summary
        summary = self._build_summary(
            pdb_id, target_chain, epitope, hotspot_rec,
            method_rankings, warnings,
        )

        return TargetAnalysis(
            pdb_id=pdb_id,
            target_chain=target_chain,
            epitope=epitope,
            hotspot_recommendation=hotspot_rec,
            method_rankings=method_rankings,
            warnings=warnings,
            summary=summary,
        )

    def _build_summary(
        self,
        pdb_id: str,
        chain: str,
        epitope: EpitopeRegion,
        hotspots: HotspotRecommendation,
        rankings: List[MethodRecommendation],
        warnings: List[str],
    ) -> str:
        """Build a concise human-readable summary paragraph."""
        parts = []

        # Epitope description
        n_res = len(epitope.residues)
        parts.append(
            f"Target {pdb_id} chain {chain} has a {epitope.classification} epitope "
            f"({n_res} residues, hydrophilicity={epitope.hydrophilicity_score:.2f}, "
            f"curvature={epitope.curvature})."
        )

        # Top method
        if rankings:
            top = rankings[0]
            parts.append(
                f"Recommended method: {top.method.name} (score={top.score:.1f}). "
                f"{top.rationale}"
            )

        # Hotspots
        parts.append(
            f"Suggested hotspots: {', '.join(str(r) for r in hotspots.residues)} "
            f"({hotspots.distribution_quality} distribution)."
        )

        # Warnings count
        if warnings:
            parts.append(
                f"Note: {len(warnings)} warning(s) detected — review before proceeding."
            )

        return " ".join(parts)


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Analyze a protein target for binder design (Phase 1)"
    )
    parser.add_argument(
        "pdb",
        help="Path to target PDB file",
    )
    parser.add_argument(
        "--chain",
        required=True,
        help="Target chain ID",
    )
    parser.add_argument(
        "--epitope",
        nargs="+",
        type=int,
        help="Epitope residue numbers (space-separated)",
    )
    parser.add_argument(
        "--time",
        choices=["fast", "standard", "unlimited"],
        default="standard",
        help="Time constraint for method selection (default: standard)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output as JSON instead of text report",
    )

    args = parser.parse_args()

    analyzer = TargetAnalyzer()

    try:
        result = analyzer.analyze(
            pdb_path=args.pdb,
            target_chain=args.chain,
            epitope_residues=args.epitope,
            time_constraint=args.time,
        )
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
    else:
        print(result.to_report())


if __name__ == "__main__":
    main()
