"""
Failure mode diagnosis engine for protein binder designs.

Applies rule-based expert knowledge to classify design failures into
Type I (misfolding) and Type II (binding) categories, identify root causes,
and recommend actionable fixes.

Operates at both individual design and pool level to detect systematic
issues in design campaigns.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional
import csv
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from filter_engine import DesignMetrics, FilterResult, load_designs_csv


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class FailureDiagnosis:
    """Diagnosis for a single design that failed or underperformed."""
    failure_type: str       # "type_i", "type_ii", "mixed"
    subtype: str            # specific failure pattern name
    evidence: List[str]     # which metrics indicated this
    root_cause: str         # expert explanation
    recommended_fixes: List[str]  # actionable fixes
    priority: str           # "critical", "high", "medium", "low"

    def to_dict(self) -> dict:
        return {
            "failure_type": self.failure_type,
            "subtype": self.subtype,
            "evidence": self.evidence,
            "root_cause": self.root_cause,
            "recommended_fixes": self.recommended_fixes,
            "priority": self.priority,
        }

    def to_report(self) -> str:
        lines = [
            f"  Failure Type: {self.failure_type} ({self.subtype})",
            f"  Priority:     {self.priority}",
            f"  Root Cause:   {self.root_cause}",
            f"  Evidence:",
        ]
        for e in self.evidence:
            lines.append(f"    - {e}")
        lines.append("  Recommended Fixes:")
        for fix in self.recommended_fixes:
            lines.append(f"    - {fix}")
        return "\n".join(lines)


@dataclass
class PoolDiagnosis:
    """Diagnosis for an entire pool of designs."""
    total_designs: int
    pass_rate: float
    failure_distribution: Dict[str, int]   # failure_type -> count
    dominant_failure: str
    individual_diagnoses: List[FailureDiagnosis]
    pool_level_recommendations: List[str]
    systematic_issues: List[str]

    def to_dict(self) -> dict:
        return {
            "total_designs": self.total_designs,
            "pass_rate": round(self.pass_rate, 4),
            "failure_distribution": self.failure_distribution,
            "dominant_failure": self.dominant_failure,
            "individual_diagnoses": [d.to_dict() for d in self.individual_diagnoses],
            "pool_level_recommendations": self.pool_level_recommendations,
            "systematic_issues": self.systematic_issues,
        }

    def to_report(self) -> str:
        lines = [
            "=" * 70,
            "POOL DIAGNOSIS REPORT",
            "=" * 70,
            f"Total designs analyzed: {self.total_designs}",
            f"Pass rate:              {self.pass_rate:.1%}",
            f"Dominant failure mode:  {self.dominant_failure}",
            "",
            "Failure distribution:",
        ]
        for ftype, count in sorted(self.failure_distribution.items(),
                                    key=lambda x: -x[1]):
            pct = count / self.total_designs * 100 if self.total_designs > 0 else 0
            lines.append(f"  {ftype}: {count} ({pct:.1f}%)")

        if self.systematic_issues:
            lines.append("")
            lines.append("SYSTEMATIC ISSUES DETECTED:")
            for issue in self.systematic_issues:
                lines.append(f"  *** {issue}")

        if self.pool_level_recommendations:
            lines.append("")
            lines.append("Pool-level recommendations:")
            for rec in self.pool_level_recommendations:
                lines.append(f"  -> {rec}")

        lines.append("")
        lines.append("-" * 70)
        lines.append(f"Individual diagnoses ({len(self.individual_diagnoses)} failures):")
        lines.append("-" * 70)
        for i, diag in enumerate(self.individual_diagnoses, 1):
            lines.append(f"\n[{i}] {diag.subtype}")
            lines.append(diag.to_report())

        lines.append("")
        lines.append("=" * 70)
        return "\n".join(lines)


# ── Diagnosis rules ──────────────────────────────────────────────────────────

DIAGNOSIS_RULES = [
    {
        "name": "severe_misfolding",
        "failure_type": "type_i",
        "condition": lambda m: m.ca_rmsd > 2.0 and m.plddt < 0.70,
        "evidence_fn": lambda m: [
            f"ca_rmsd={m.ca_rmsd:.2f} (>2.0)",
            f"plddt={m.plddt:.3f} (<0.70)",
        ],
        "root_cause": "Designed backbone not achievable by any sequence",
        "recommended_fixes": [
            "Increase ProteinMPNN sampling (num_seqs 100->500)",
            "Lower design temperature (0.1->0.05)",
            "Simplify backbone topology",
            "Check for strained geometry",
        ],
        "priority": "critical",
    },
    {
        "name": "partial_misfolding",
        "failure_type": "type_i",
        "condition": lambda m: m.ca_rmsd > 2.0 and m.plddt >= 0.70,
        "evidence_fn": lambda m: [
            f"ca_rmsd={m.ca_rmsd:.2f} (>2.0)",
            f"plddt={m.plddt:.3f} (>=0.70, confident but wrong fold)",
        ],
        "root_cause": "Sequence folds to alternative stable state",
        "recommended_fixes": [
            "Add sequence constraints at key positions",
            "Use AF2 with more recycling iterations",
            "Consider disulfide bonds to stabilize fold",
        ],
        "priority": "high",
    },
    {
        "name": "no_binding_predicted",
        "failure_type": "type_ii",
        "condition": lambda m: m.ca_rmsd <= 2.0 and m.iptm < 0.6,
        "evidence_fn": lambda m: [
            f"ca_rmsd={m.ca_rmsd:.2f} (<=2.0, fold OK)",
            f"iptm={m.iptm:.3f} (<0.6, no binding signal)",
        ],
        "root_cause": "Interface not recognized as stable complex by AF2",
        "recommended_fixes": [
            "Increase interface area",
            "Improve shape complementarity",
            "Try different hotspot residues",
            "Consider switching to BindCraft",
        ],
        "priority": "high",
    },
    {
        "name": "interface_too_small",
        "failure_type": "type_ii",
        "condition": lambda m: (
            m.ca_rmsd <= 2.0
            and 0.6 <= m.iptm < 0.8
            and m.interface_area < 850
        ),
        "evidence_fn": lambda m: [
            f"iptm={m.iptm:.3f} (gray zone 0.6-0.8)",
            f"interface_area={m.interface_area:.0f} (<850 A^2)",
        ],
        "root_cause": "Insufficient buried surface area",
        "recommended_fixes": [
            "Extend binder length",
            "Use longer contigs",
            "Add flanking helices",
        ],
        "priority": "medium",
    },
    {
        "name": "poor_hotspot_engagement",
        "failure_type": "type_ii",
        "condition": lambda m: m.hotspot_contact_rate < 0.20,
        "evidence_fn": lambda m: [
            f"hotspot_contact_rate={m.hotspot_contact_rate:.2f} (<0.20)",
        ],
        "root_cause": "Designed binder not making intended contacts",
        "recommended_fixes": [
            "Increase number of hotspots (+2-3)",
            "Redistribute across epitope",
            "Check hotspot accessibility",
        ],
        "priority": "high",
    },
    {
        "name": "hotspot_over_constraint",
        "failure_type": "type_ii",
        "condition": lambda m: m.hotspot_contact_rate > 0.50,
        "evidence_fn": lambda m: [
            f"hotspot_contact_rate={m.hotspot_contact_rate:.2f} (>0.50)",
        ],
        "root_cause": "Too many constraints limiting design freedom",
        "recommended_fixes": [
            "Reduce hotspot count by 2-3",
            "Remove least important hotspots",
            "Use softer distance-based constraints",
        ],
        "priority": "medium",
    },
    {
        "name": "weak_binding",
        "failure_type": "type_ii",
        "condition": lambda m: (
            m.ca_rmsd <= 2.0
            and 0.6 <= m.iptm < 0.8
            and m.interface_area >= 850
        ),
        "evidence_fn": lambda m: [
            f"iptm={m.iptm:.3f} (gray zone 0.6-0.8)",
            f"interface_area={m.interface_area:.0f} (>=850, adequate)",
            f"ca_rmsd={m.ca_rmsd:.2f} (<=2.0, fold OK)",
        ],
        "root_cause": "Interface geometry suboptimal -- may bind weakly",
        "recommended_fixes": [
            "Optimize interface packing",
            "Add hydrophobic contacts",
            "Try affinity maturation with ProteinMPNN",
        ],
        "priority": "medium",
    },
]


# ── Diagnosis Engine ─────────────────────────────────────────────────────────

class DiagnosisEngine:
    """
    Rule-based diagnosis engine for protein binder design failures.

    Classifies individual design failures into Type I (misfolding) and
    Type II (binding) categories, and detects systematic patterns across
    design pools.
    """

    def __init__(self, rules: Optional[List[dict]] = None):
        self.rules = rules or DIAGNOSIS_RULES

    def diagnose_design(self, metrics: DesignMetrics) -> Optional[FailureDiagnosis]:
        """
        Diagnose a single design.

        Returns the first matching failure diagnosis, or None if the design
        does not match any failure pattern (i.e., it looks healthy).
        """
        for rule in self.rules:
            if rule["condition"](metrics):
                return FailureDiagnosis(
                    failure_type=rule["failure_type"],
                    subtype=rule["name"],
                    evidence=rule["evidence_fn"](metrics),
                    root_cause=rule["root_cause"],
                    recommended_fixes=list(rule["recommended_fixes"]),
                    priority=rule["priority"],
                )
        return None

    def diagnose_pool(
        self,
        designs: List[DesignMetrics],
        filter_results: Optional[Dict[str, List[FilterResult]]] = None,
    ) -> PoolDiagnosis:
        """
        Diagnose an entire pool of designs.

        Runs individual diagnosis on each design, then detects systematic
        patterns and generates pool-level recommendations.

        Args:
            designs: List of design metrics to diagnose.
            filter_results: Optional pre-computed filter results (used to
                determine pass rate if provided).
        """
        individual_diagnoses = []
        failure_distribution: Dict[str, int] = {}

        for metrics in designs:
            diag = self.diagnose_design(metrics)
            if diag is not None:
                individual_diagnoses.append(diag)
                failure_distribution[diag.subtype] = (
                    failure_distribution.get(diag.subtype, 0) + 1
                )

        # Compute pass rate
        if filter_results is not None:
            total = sum(len(v) for v in filter_results.values())
            n_pass = len(filter_results.get("pass", []))
            pass_rate = n_pass / total if total > 0 else 0.0
        else:
            n_failed = len(individual_diagnoses)
            total = len(designs)
            pass_rate = (total - n_failed) / total if total > 0 else 0.0

        # Determine dominant failure
        if failure_distribution:
            dominant_failure = max(failure_distribution, key=failure_distribution.get)
        else:
            dominant_failure = "none"

        # Detect systematic patterns
        systematic_issues = self._detect_systematic_patterns(
            individual_diagnoses, total
        )

        # Pool-level recommendations
        pool_recs = self._pool_recommendations(
            individual_diagnoses, total, pass_rate
        )

        return PoolDiagnosis(
            total_designs=total,
            pass_rate=pass_rate,
            failure_distribution=failure_distribution,
            dominant_failure=dominant_failure,
            individual_diagnoses=individual_diagnoses,
            pool_level_recommendations=pool_recs,
            systematic_issues=systematic_issues,
        )

    def _detect_systematic_patterns(
        self,
        diagnoses: List[FailureDiagnosis],
        total_designs: int,
    ) -> List[str]:
        """
        Detect systematic issues when >60% of designs share the same
        failure mode.
        """
        if total_designs == 0 or not diagnoses:
            return []

        issues = []

        # Count by failure type (type_i vs type_ii)
        type_counts: Dict[str, int] = {}
        subtype_counts: Dict[str, int] = {}
        for d in diagnoses:
            type_counts[d.failure_type] = type_counts.get(d.failure_type, 0) + 1
            subtype_counts[d.subtype] = subtype_counts.get(d.subtype, 0) + 1

        for ftype, count in type_counts.items():
            if count / total_designs > 0.60:
                pct = count / total_designs * 100
                issues.append(
                    f"{ftype} failures dominate pool ({pct:.0f}% of designs) "
                    f"-- indicates a systematic {ftype} problem"
                )

        for subtype, count in subtype_counts.items():
            if count / total_designs > 0.60:
                pct = count / total_designs * 100
                issues.append(
                    f"'{subtype}' affects {pct:.0f}% of designs -- "
                    f"this is a systematic pattern requiring campaign-level changes"
                )

        return issues

    def _pool_recommendations(
        self,
        diagnoses: List[FailureDiagnosis],
        total_designs: int,
        pass_rate: float,
    ) -> List[str]:
        """Generate pool-level recommendations based on failure patterns."""
        recs = []

        if total_designs == 0:
            return recs

        type_counts: Dict[str, int] = {}
        subtype_counts: Dict[str, int] = {}
        for d in diagnoses:
            type_counts[d.failure_type] = type_counts.get(d.failure_type, 0) + 1
            subtype_counts[d.subtype] = subtype_counts.get(d.subtype, 0) + 1

        type_i_frac = type_counts.get("type_i", 0) / total_designs
        type_ii_frac = type_counts.get("type_ii", 0) / total_designs
        hotspot_frac = subtype_counts.get("poor_hotspot_engagement", 0) / total_designs

        if type_i_frac > 0.60:
            recs.append(
                "Backbone generation parameters need adjustment -- "
                ">60% of designs show Type I (misfolding) failures"
            )

        if type_ii_frac > 0.60:
            recs.append(
                "Interface design strategy needs revision -- "
                ">60% of designs show Type II (binding) failures"
            )

        if hotspot_frac > 0.40:
            recs.append(
                "Hotspot selection needs reconsideration -- "
                ">40% of designs show poor hotspot engagement"
            )

        if pass_rate < 0.05:
            recs.append(
                "Consider switching design method entirely -- "
                "overall pass rate is below 5%"
            )

        return recs


# ── CLI entry point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Diagnose failure modes in protein binder design pools"
    )
    parser.add_argument("input", help="CSV file with design metrics")
    parser.add_argument(
        "--output", default=None,
        help="Output file (JSON). If not provided, prints report to stdout."
    )
    parser.add_argument(
        "--format", choices=["report", "json"], default="report",
        help="Output format (default: report)"
    )
    args = parser.parse_args()

    # Load designs
    designs = load_designs_csv(args.input)
    print(f"Loaded {len(designs)} designs from {args.input}")

    # Run diagnosis
    engine = DiagnosisEngine()
    pool_diag = engine.diagnose_pool(designs)

    # Output
    if args.format == "json" or (args.output and args.output.endswith(".json")):
        output_data = pool_diag.to_dict()
        if args.output:
            with open(args.output, "w") as f:
                json.dump(output_data, f, indent=2)
            print(f"Diagnosis saved to {args.output}")
        else:
            print(json.dumps(output_data, indent=2))
    else:
        report = pool_diag.to_report()
        if args.output:
            with open(args.output, "w") as f:
                f.write(report)
            print(f"Report saved to {args.output}")
        else:
            print(report)


if __name__ == "__main__":
    main()
