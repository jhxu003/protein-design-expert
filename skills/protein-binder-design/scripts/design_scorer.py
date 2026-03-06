"""
Composite scoring engine for protein binder designs.

Combines multiple AF2/AF3 prediction metrics into a single composite score
using expert-calibrated weights and normalization ranges. Assigns designs
to quality tiers (gold/silver/bronze/reject) and produces ranked outputs.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import csv
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from filter_engine import DesignMetrics, load_designs_csv


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class ScoringWeights:
    """
    Weights for composite score computation.

    Defaults reflect expert consensus:
    - ipSAE is the best single predictor of binding success
    - ipTM is the second best predictor
    - pLDDT catches misfolding (Type I failures)
    - Remaining metrics provide complementary signal
    """
    ipsae: float = 0.30
    iptm: float = 0.25
    plddt: float = 0.15
    pae: float = 0.10
    interface_area: float = 0.10
    shape_complementarity: float = 0.05
    hotspot_contact: float = 0.05

    def to_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items()}

    @classmethod
    def from_dict(cls, d: dict) -> "ScoringWeights":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})

    def validate(self) -> bool:
        """Check that weights sum to approximately 1.0."""
        total = sum(self.__dict__.values())
        return abs(total - 1.0) < 0.01


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

    Normalizes each metric to 0-1, applies expert-calibrated weights,
    and assigns quality tiers.
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
