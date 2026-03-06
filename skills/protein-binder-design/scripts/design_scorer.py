"""
Composite scoring engine for protein binder designs.

Combines multiple AF2/AF3 prediction metrics into a single composite score
using weights informed by large-scale experimental validation data. Assigns
designs to quality tiers (gold/silver/bronze/reject) and produces ranked
outputs.

Data source: Overath et al. 2025, "Scoring metrics for de novo protein
binder design across 15 targets" -- meta-analysis of 3,760 experimentally
tested designs (436 confirmed binders, 11.6% overall success rate) across
15 diverse targets. Per-target success rates range from 2.1% to 57.3%.
Repository: DigBioLab/de_novo_binder_scoring
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import csv
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from filter_engine import DesignMetrics, load_designs_csv


# ── Data-driven feature ranking ──────────────────────────────────────────────
#
# Median Average Precision (AP) across 15 targets from Overath et al. 2025.
# These values quantify each feature's ability to discriminate true binders
# from non-binders in a target-independent manner.
#
# Best single predictor: AF3 ipSAE_min (median AP = 0.5399)
# Best 3-feature logistic regression: ipSAE_min + RMSD_binder + interface_sc
#   -> median AP = 0.573, precision@F1 = 0.538
#   (vs. ipSAE_min alone: median AP = 0.540, precision@F1 = 0.465)
# Best interaction feature: AF3 LIS * input_interface_sc (median AP = 0.5752)
# Optimal filter thresholds: RMSD_binder < 3.73 A, input_interface_sc > 0.62

DATA_DRIVEN_FEATURE_RANKING = {
    # Feature name (AF3)                        Median AP
    "af3_ipSAE_min":                             0.5399,
    "af3_iptm_avg":                              0.5212,
    "af3_ipSAE_d0chn":                           0.5170,
    "af3_iptm_model_0":                          0.5098,
    "af3_min_pae_contact":                        0.4898,
    "af3_pDockQ2_max":                           0.4524,
    "af3_ipae":                                  0.4358,
    "RMSD_chA_aft_chB_align_input_af3":          0.3038,
    "af3_rosetta_interface_dG":                   0.3092,
    "af3_rosetta_interface_dG_SASA_ratio":        0.2849,
}


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class ScoringWeights:
    """
    Weights for composite score computation.

    Defaults are informed by the Overath et al. 2025 meta-analysis of 3,760
    de novo binder designs across 15 targets (436 confirmed binders):

    - ipSAE (0.30): Best single predictor, median AP = 0.5399
    - ipTM (0.25): Second best predictor, median AP = 0.5212
    - pLDDT (0.12): Catches misfolding but was NOT among top individual
      predictors in the meta-analysis; weight reduced accordingly
    - pAE (0.08): Related to af3_ipae (median AP = 0.4358) and
      af3_min_pae_contact (median AP = 0.4898)
    - interface_area (0.08): Interface area/energy metrics ranked lower
      (dG median AP = 0.3092); weight reduced
    - shape_complementarity (0.10): Part of the best 3-feature logistic
      regression model (ipSAE_min + RMSD_binder + interface_sc, median
      AP = 0.573); weight increased from 0.05
    - hotspot_contact (0.07): Design-specific contact quality; not directly
      benchmarked in the meta-analysis but captures target-specific intent
    """
    ipsae: float = 0.30
    iptm: float = 0.25
    plddt: float = 0.12
    pae: float = 0.08
    interface_area: float = 0.08
    shape_complementarity: float = 0.10
    hotspot_contact: float = 0.07

    def to_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items()}

    @classmethod
    def from_dict(cls, d: dict) -> "ScoringWeights":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})

    def validate(self) -> bool:
        """Check that weights sum to approximately 1.0."""
        total = sum(self.__dict__.values())
        return abs(total - 1.0) < 0.01

    @classmethod
    def data_driven(cls) -> "ScoringWeights":
        """
        Alternative weights proportional to median AP values from the
        Overath et al. 2025 meta-analysis (3,760 designs, 15 targets).

        Mapped features and their median AP:
          ipSAE  -> af3_ipSAE_min:       0.5399
          ipTM   -> af3_iptm_avg:        0.5212
          pLDDT  -> (not ranked; use af3_pDockQ2_max as proxy): 0.4524
          pAE    -> af3_min_pae_contact:  0.4898
          interface_area -> af3_rosetta_interface_dG: 0.3092
          shape_complementarity -> input_interface_sc (in best 3-feature
                                   model; use af3_ipSAE_d0chn proxy): 0.5170
          hotspot_contact -> (design-specific; use af3_ipae): 0.4358

        Raw APs are normalized to sum to 1.0.
        """
        raw = {
            "ipsae": 0.5399,
            "iptm": 0.5212,
            "plddt": 0.4524,
            "pae": 0.4898,
            "interface_area": 0.3092,
            "shape_complementarity": 0.5170,
            "hotspot_contact": 0.4358,
        }
        total = sum(raw.values())
        normalized = {k: round(v / total, 4) for k, v in raw.items()}
        # Adjust last weight so they sum exactly to 1.0
        diff = 1.0 - sum(normalized.values())
        normalized["hotspot_contact"] = round(normalized["hotspot_contact"] + diff, 4)
        return cls(**normalized)


@dataclass
class ScoredDesign:
    """A design with its composite score, tier, and rank."""
    design: DesignMetrics
    composite_score: float   # 0-1
    component_scores: Dict[str, float]
    tier: str                # "gold", "silver", "bronze", "reject"
    rank: int
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        d = self.design.to_dict()
        d["composite_score"] = round(self.composite_score, 4)
        d["tier"] = self.tier
        d["rank"] = self.rank
        d["component_scores"] = {
            k: round(v, 4) for k, v in self.component_scores.items()
        }
        d["notes"] = self.notes
        return d


# ── Normalization configuration ──────────────────────────────────────────────

# Each metric: (min_val, max_val, direction)
# direction: "higher_better" or "lower_better"
NORMALIZATION_RANGES = {
    "ipsae":                (0.0,   1.0,  "lower_better"),
    "iptm":                 (0.0,   1.0,  "higher_better"),
    "plddt":                (0.5,   1.0,  "higher_better"),
    "pae":                  (0.0,  30.0,  "lower_better"),
    "interface_area":       (500.0, 2500.0, "higher_better"),
    "shape_complementarity": (0.0,  1.0,  "higher_better"),
    "hotspot_contact":      (0.0,   0.6,  "higher_better"),
}

# Tier thresholds
TIER_THRESHOLDS = {
    "gold":   0.80,
    "silver": 0.65,
    "bronze": 0.50,
}


# ── Scoring Engine ───────────────────────────────────────────────────────────

class DesignScorer:
    """
    Composite scoring engine for protein binder designs.

    Normalizes each metric to 0-1, applies weights informed by the
    Overath et al. 2025 meta-analysis, and assigns quality tiers.
    Use ScoringWeights() for default weights or
    ScoringWeights.data_driven() for weights proportional to median AP.
    """

    def __init__(self, weights: Optional[ScoringWeights] = None):
        self.weights = weights or ScoringWeights()

    def score_designs(self, designs: List[DesignMetrics]) -> List[ScoredDesign]:
        """
        Score, rank, and tier-assign a list of designs.

        Returns a list of ScoredDesign sorted by composite_score descending.
        """
        scored = []
        for design in designs:
            composite, components, notes = self._compute_composite(design)
            tier = self._assign_tier(composite)
            scored.append(ScoredDesign(
                design=design,
                composite_score=composite,
                component_scores=components,
                tier=tier,
                rank=0,  # assigned after sorting
                notes=notes,
            ))

        # Sort by composite score descending
        scored.sort(key=lambda s: s.composite_score, reverse=True)

        # Assign ranks
        for i, s in enumerate(scored, 1):
            s.rank = i

        return scored

    def _normalize(self, value: float, metric: str) -> float:
        """
        Clip and normalize a metric value to 0-1.

        For higher_better metrics: (value - min) / (max - min)
        For lower_better metrics: (max - value) / (max - min)
        """
        if metric not in NORMALIZATION_RANGES:
            return 0.0

        min_val, max_val, direction = NORMALIZATION_RANGES[metric]

        # Clip to range
        clipped = max(min_val, min(max_val, value))

        span = max_val - min_val
        if span == 0:
            return 0.0

        if direction == "higher_better":
            return (clipped - min_val) / span
        else:
            return (max_val - clipped) / span

    def _compute_composite(
        self, metrics: DesignMetrics
    ) -> Tuple[float, Dict[str, float], List[str]]:
        """
        Compute composite score from individual metric scores.

        Returns (composite_score, component_scores_dict, notes).
        Handles missing optional metrics by redistributing their weight.
        """
        w = self.weights
        notes = []

        # Map metric names to (value, weight)
        metric_map = {
            "ipsae":                (metrics.ipsae_min, w.ipsae),
            "iptm":                 (metrics.iptm, w.iptm),
            "plddt":                (metrics.plddt, w.plddt),
            "pae":                  (metrics.pae, w.pae),
            "interface_area":       (metrics.interface_area, w.interface_area),
            "shape_complementarity": (metrics.shape_complementarity, w.shape_complementarity),
            "hotspot_contact":      (metrics.hotspot_contact_rate, w.hotspot_contact),
        }

        # Compute normalized scores and handle missing values
        component_scores = {}
        active_weight_sum = 0.0
        weighted_sum = 0.0

        for metric_name, (value, weight) in metric_map.items():
            if value is None:
                notes.append(f"{metric_name}: missing, weight redistributed")
                continue

            norm = self._normalize(value, metric_name)
            component_scores[metric_name] = norm
            weighted_sum += norm * weight
            active_weight_sum += weight

        # Normalize by active weight sum to handle missing metrics
        if active_weight_sum > 0:
            composite = weighted_sum / active_weight_sum
        else:
            composite = 0.0
            notes.append("No metrics available for scoring")

        return composite, component_scores, notes

    def _assign_tier(self, composite_score: float) -> str:
        """Assign a quality tier based on composite score."""
        if composite_score >= TIER_THRESHOLDS["gold"]:
            return "gold"
        elif composite_score >= TIER_THRESHOLDS["silver"]:
            return "silver"
        elif composite_score >= TIER_THRESHOLDS["bronze"]:
            return "bronze"
        else:
            return "reject"

    def export_rankings(
        self,
        scored: List[ScoredDesign],
        path: str,
        format: str = "csv",
    ):
        """
        Export ranked designs to CSV or JSON.

        Args:
            scored: List of ScoredDesign (already ranked).
            path: Output file path.
            format: "csv" or "json".
        """
        if format == "json":
            data = [s.to_dict() for s in scored]
            with open(path, "w") as f:
                json.dump(data, f, indent=2)
        else:
            if not scored:
                return
            rows = []
            for s in scored:
                row = {
                    "rank": s.rank,
                    "design_id": s.design.design_id,
                    "composite_score": round(s.composite_score, 4),
                    "tier": s.tier,
                    "ipsae_min": s.design.ipsae_min,
                    "iptm": s.design.iptm,
                    "plddt": s.design.plddt,
                    "pae": s.design.pae,
                    "interface_area": s.design.interface_area,
                    "ca_rmsd": s.design.ca_rmsd,
                    "hotspot_contact_rate": s.design.hotspot_contact_rate,
                }
                if s.design.shape_complementarity is not None:
                    row["shape_complementarity"] = s.design.shape_complementarity
                if s.design.binder_length is not None:
                    row["binder_length"] = s.design.binder_length
                # Add component scores
                for comp_name, comp_val in s.component_scores.items():
                    row[f"score_{comp_name}"] = round(comp_val, 4)
                if s.notes:
                    row["notes"] = "; ".join(s.notes)
                rows.append(row)

            fieldnames = list(rows[0].keys())
            # Ensure all rows have all keys
            all_keys = set()
            for r in rows:
                all_keys.update(r.keys())
            fieldnames = list(rows[0].keys())
            for k in sorted(all_keys - set(fieldnames)):
                fieldnames.append(k)

            with open(path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)

    def print_summary(self, scored: List[ScoredDesign]):
        """Print tier distribution and top designs to stdout."""
        tier_counts = {"gold": 0, "silver": 0, "bronze": 0, "reject": 0}
        for s in scored:
            tier_counts[s.tier] = tier_counts.get(s.tier, 0) + 1

        total = len(scored)
        print("=" * 60)
        print("DESIGN SCORING SUMMARY")
        print("=" * 60)
        print(f"Total designs scored: {total}")
        print()
        print("Tier distribution:")
        for tier in ["gold", "silver", "bronze", "reject"]:
            count = tier_counts[tier]
            pct = count / total * 100 if total > 0 else 0
            bar = "#" * int(pct / 2)
            print(f"  {tier:>7s}: {count:>4d} ({pct:5.1f}%) {bar}")

        print()
        top_n = min(10, total)
        print(f"Top {top_n} designs:")
        print(f"  {'Rank':>4s}  {'Design ID':<25s}  {'Score':>6s}  {'Tier':<7s}  "
              f"{'ipSAE':>6s}  {'ipTM':>5s}  {'pLDDT':>6s}")
        print("  " + "-" * 85)
        for s in scored[:top_n]:
            print(
                f"  {s.rank:>4d}  {s.design.design_id:<25s}  "
                f"{s.composite_score:>6.3f}  {s.tier:<7s}  "
                f"{s.design.ipsae_min:>6.3f}  {s.design.iptm:>5.3f}  "
                f"{s.design.plddt:>6.3f}"
            )
        print("=" * 60)


# ── CLI entry point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Score and rank protein binder designs using composite metrics"
    )
    parser.add_argument("input", help="CSV file with design metrics (e.g., filtered_designs.csv)")
    parser.add_argument(
        "--output", default="ranked.csv",
        help="Output file path (default: ranked.csv)"
    )
    parser.add_argument(
        "--format", choices=["csv", "json"], default=None,
        help="Output format (inferred from extension if not specified)"
    )
    parser.add_argument(
        "--weights", default=None,
        help="JSON file with custom scoring weights"
    )
    parser.add_argument(
        "--summary", action="store_true", default=True,
        help="Print summary to stdout (default: True)"
    )
    parser.add_argument(
        "--no-summary", action="store_false", dest="summary",
        help="Suppress summary output"
    )
    args = parser.parse_args()

    # Load designs
    designs = load_designs_csv(args.input)
    print(f"Loaded {len(designs)} designs from {args.input}")

    # Load custom weights if provided
    weights = None
    if args.weights:
        with open(args.weights, "r") as f:
            weights_dict = json.load(f)
        weights = ScoringWeights.from_dict(weights_dict)
        if not weights.validate():
            total = sum(weights.__dict__.values())
            print(f"WARNING: Custom weights sum to {total:.3f}, not 1.0")
        print(f"Using custom weights from {args.weights}")

    # Score designs
    scorer = DesignScorer(weights=weights)
    scored = scorer.score_designs(designs)

    # Print summary
    if args.summary:
        scorer.print_summary(scored)

    # Determine output format
    out_format = args.format
    if out_format is None:
        if args.output.endswith(".json"):
            out_format = "json"
        else:
            out_format = "csv"

    # Export
    scorer.export_rankings(scored, args.output, format=out_format)
    print(f"\nRankings saved to {args.output} ({out_format})")


if __name__ == "__main__":
    main()
