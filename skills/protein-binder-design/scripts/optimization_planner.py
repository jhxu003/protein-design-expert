"""
Phase 4 iterative optimization planning engine for protein binder design.

Analyzes pool diagnosis results and scored designs to produce an actionable
optimization plan: parameter adjustments, hotspot surgery, method switching,
and convergence assessment. Encodes expert heuristics from competitive binder
design campaigns (Adaptyv Bio, BenchBB).

Imports:
    - PoolDiagnosis from diagnosis_engine.py
    - ScoredDesign from design_scorer.py
"""

from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple
import json
import math
import sys
import os
import argparse

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from diagnosis_engine import PoolDiagnosis
from design_scorer import ScoredDesign


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class ParameterAdjustment:
    """A single recommended parameter change for the next optimization round."""
    parameter: str          # e.g., "num_seqs", "mpnn_temperature"
    current_value: str      # e.g., "100"
    recommended_value: str  # e.g., "500"
    rationale: str
    expected_impact: str    # e.g., "Reduce Type I failures by ~30%"

    def to_dict(self) -> Dict:
        return asdict(self)


@dataclass
class OptimizationPlan:
    """Complete optimization plan for the next round of binder design."""
    round_number: int
    convergence_status: str    # "improving", "plateaued", "diverging"
    should_continue: bool
    should_switch_method: bool
    switch_to: Optional[str]
    switch_rationale: Optional[str]
    parameter_adjustments: List[ParameterAdjustment]
    hotspot_changes: List[str]  # e.g., ["add res 45", "remove res 72"]
    designs_to_generate: int
    estimated_improvement: str
    convergence_criteria_met: List[str]
    convergence_criteria_remaining: List[str]

    def to_dict(self) -> Dict:
        d = asdict(self)
        return d

    def to_report(self) -> str:
        """Human-readable optimization plan report."""
        lines = []
        lines.append("=" * 72)
        lines.append(f"  OPTIMIZATION PLAN — Round {self.round_number}")
        lines.append("=" * 72)
        lines.append("")

        # Convergence status
        status_icon = {
            "improving": "[+]",
            "plateaued": "[=]",
            "diverging": "[-]",
        }.get(self.convergence_status, "[?]")
        lines.append(f"Convergence: {status_icon} {self.convergence_status}")
        lines.append(f"Continue: {'Yes' if self.should_continue else 'No — stopping'}")
        lines.append("")

        # Method switch
        if self.should_switch_method:
            lines.append(f"METHOD SWITCH RECOMMENDED: → {self.switch_to}")
            lines.append(f"  Rationale: {self.switch_rationale}")
            lines.append("")

        # Parameter adjustments
        if self.parameter_adjustments:
            lines.append("Parameter Adjustments:")
            lines.append("-" * 50)
            for adj in self.parameter_adjustments:
                lines.append(f"  {adj.parameter}:")
                lines.append(f"    Current:     {adj.current_value}")
                lines.append(f"    Recommended: {adj.recommended_value}")
                lines.append(f"    Rationale:   {adj.rationale}")
                lines.append(f"    Impact:      {adj.expected_impact}")
                lines.append("")

        # Hotspot changes
        if self.hotspot_changes:
            lines.append("Hotspot Surgery:")
            lines.append("-" * 50)
            for change in self.hotspot_changes:
                lines.append(f"  • {change}")
            lines.append("")

        # Design count
        lines.append(f"Designs to generate: {self.designs_to_generate}")
        lines.append(f"Estimated improvement: {self.estimated_improvement}")
        lines.append("")

        # Convergence criteria
        if self.convergence_criteria_met:
            lines.append("Convergence criteria MET:")
            for c in self.convergence_criteria_met:
                lines.append(f"  [x] {c}")
        if self.convergence_criteria_remaining:
            lines.append("Convergence criteria REMAINING:")
            for c in self.convergence_criteria_remaining:
                lines.append(f"  [ ] {c}")

        lines.append("")
        lines.append("=" * 72)
        return "\n".join(lines)


@dataclass
class RoundHistory:
    """Tracks results from a single optimization round for trend analysis."""
    round: int
    method: str
    params: Dict
    results: Dict  # {total, passed, gold, silver, bronze, pass_rate, best_score}
    dominant_failure: str
    convergence_status: str

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "RoundHistory":
        return cls(
            round=d["round"],
            method=d["method"],
            params=d["params"],
            results=d["results"],
            dominant_failure=d["dominant_failure"],
            convergence_status=d["convergence_status"],
        )


# ── Constants ────────────────────────────────────────────────────────────────

CONVERGENCE_CRITERIA = {
    "min_gold_tier": 5,           # At least 5 gold-tier designs
    "min_pass_rate": 0.10,        # At least 10% pass rate
    "max_rounds": 5,              # Stop after 5 rounds
    "improvement_threshold": 0.05, # Must improve >5% per round to continue
}

ADJUSTMENT_RULES = {
    "type_i_dominant": {
        "severe_misfolding": [
            ("num_seqs", "increase 2X", "More sequence diversity to find foldable variants"),
            ("mpnn_temperature", "decrease to 0.05", "More conservative sequence sampling"),
            ("backbone_noise", "decrease by 0.5", "Reduce backbone perturbation"),
            ("contig_topology", "simplify", "Remove complex topological features"),
        ],
        "partial_misfolding": [
            ("num_seqs", "increase 2X", "Sample more sequences"),
            ("sequence_constraints", "add at key positions", "Pin critical folding residues"),
            ("disulfide_bonds", "consider adding", "Stabilize fold with covalent constraints"),
            ("af2_recycling", "increase iterations", "Better structure prediction for validation"),
        ],
    },
    "type_ii_dominant": {
        "no_binding_predicted": [
            ("hotspot_residues", "change selection", "Current hotspots may not be optimal"),
            ("interface_design", "increase contact area", "Need more interfacial contacts"),
            ("method", "consider BindCraft", "BindCraft has higher binding success rates"),
        ],
        "interface_too_small": [
            ("binder_length", "increase by 30-50 residues", "Longer binders can make more contacts"),
            ("contig_range", "widen range", "Allow larger binder scaffolds"),
            ("hotspot_count", "increase by 2-3", "More hotspots drive larger interfaces"),
        ],
        "poor_hotspot_engagement": [
            ("hotspot_count", "increase by 2-3", "More hotspots increase chance of engagement"),
            ("hotspot_distribution", "redistribute across epitope", "Spread hotspots more evenly"),
            ("hotspot_accessibility", "verify surface exposure", "Ensure hotspots are solvent-accessible"),
        ],
        "hotspot_over_constraint": [
            ("hotspot_count", "reduce by 2-3", "Too many constraints limit design freedom"),
            ("hotspot_weight", "reduce constraint strength", "Softer constraints allow more exploration"),
        ],
        "weak_binding": [
            ("interface_optimization", "add hydrophobic contacts", "Hydrophobic packing improves affinity"),
            ("affinity_maturation", "run ProteinMPNN refinement", "Sequence optimization at interface"),
            ("scoring_function", "weight shape complementarity higher", "Better geometric fit needed"),
        ],
    },
}

METHOD_SWITCH_RULES = {
    # (current_method, condition) -> (recommended_method, rationale)
    ("RFdiffusion", "plateaued"): (
        "BindCraft",
        "BindCraft has 10X higher success rates and may overcome current plateau",
    ),
    ("RFdiffusion", "type_ii_dominant"): (
        "PXDesign",
        "PXDesign's dual modality may find better interfaces",
    ),
    ("RFdiffusion", "edge_strand_target"): (
        "beta-pairing-RFdiffusion",
        "Beta-pairing gives 9.2X improvement for edge-strand targets",
    ),
    ("BindCraft", "plateaued"): (
        "PXDesign",
        "PXDesign offers different fold exploration",
    ),
    ("BindCraft", "type_i_dominant"): (
        "RFdiffusion",
        "RFdiffusion backbone generation may produce more stable folds",
    ),
    ("PXDesign", "plateaued"): (
        "BindCraft",
        "Try BindCraft's hallucination approach",
    ),
}


# ── OptimizationPlanner ─────────────────────────────────────────────────────

class OptimizationPlanner:
    """
    Plans the next round of iterative binder optimization.

    Given a pool diagnosis, scored designs, and current parameters, produces
    an OptimizationPlan with concrete parameter adjustments, hotspot changes,
    method switching recommendations, and convergence assessment.
    """

    def __init__(
        self,
        convergence_criteria: Optional[Dict] = None,
        adjustment_rules: Optional[Dict] = None,
        method_switch_rules: Optional[Dict] = None,
    ):
        self.convergence_criteria = convergence_criteria or CONVERGENCE_CRITERIA
        self.adjustment_rules = adjustment_rules or ADJUSTMENT_RULES
        self.method_switch_rules = method_switch_rules or METHOD_SWITCH_RULES

    def plan_next_round(
        self,
        round_number: int,
        pool_diagnosis: PoolDiagnosis,
        scored_designs: List[ScoredDesign],
        current_params: Dict,
        history: Optional[List[RoundHistory]] = None,
    ) -> OptimizationPlan:
        """
        Produce a full optimization plan for the next design round.

        Args:
            round_number: The upcoming round number (1-indexed).
            pool_diagnosis: Diagnosis of the current design pool.
            scored_designs: List of scored designs from the current round.
            current_params: Current design parameters (method, hotspots, etc.).
            history: Optional list of previous round histories for trend analysis.

        Returns:
            OptimizationPlan with all recommendations.
        """
        history = history or []

        # Extract current round results
        current_results = self._extract_results(scored_designs)

        # Get previous round results for trend comparison
        previous_results = history[-1].results if history else None

        # Assess convergence trend
        convergence_status = self._assess_convergence(current_results, previous_results)

        # Check convergence criteria
        criteria_met, criteria_remaining = self._check_convergence_criteria(
            current_results, round_number, history
        )

        # Decide whether to continue
        should_continue = self._should_continue(
            round_number, convergence_status, criteria_met, criteria_remaining
        )

        # Check for method switch recommendation
        current_method = current_params.get("method", "RFdiffusion")
        should_switch, switch_to, switch_rationale = self._recommend_method_switch(
            current_method, pool_diagnosis, convergence_status
        )

        # Plan parameter adjustments
        parameter_adjustments = self._plan_parameter_adjustments(
            pool_diagnosis, current_params
        )

        # Plan hotspot surgery
        current_hotspots = current_params.get("hotspot_residues", [])
        hotspot_changes = self._plan_hotspot_surgery(pool_diagnosis, current_hotspots)

        # Estimate designs needed
        current_pass_rate = current_results.get("pass_rate", 0.0)
        target_gold = self.convergence_criteria["min_gold_tier"]
        designs_to_generate = self._estimate_designs_needed(
            current_pass_rate, target_gold
        )

        # Estimate improvement
        estimated_improvement = self._estimate_improvement(
            convergence_status, parameter_adjustments, should_switch
        )

        return OptimizationPlan(
            round_number=round_number,
            convergence_status=convergence_status,
            should_continue=should_continue,
            should_switch_method=should_switch,
            switch_to=switch_to,
            switch_rationale=switch_rationale,
            parameter_adjustments=parameter_adjustments,
            hotspot_changes=hotspot_changes,
            designs_to_generate=designs_to_generate,
            estimated_improvement=estimated_improvement,
            convergence_criteria_met=criteria_met,
            convergence_criteria_remaining=criteria_remaining,
        )

    def _assess_convergence(
        self,
        current_results: Dict,
        previous_results: Optional[Dict],
    ) -> str:
        """
        Assess the convergence trend between rounds.

        Returns:
            "improving" if pass rate improved by more than the threshold,
            "plateaued" if change is within threshold,
            "diverging" if pass rate decreased.
        """
        if previous_results is None:
            # First round — no trend available, assume improving
            return "improving"

        current_rate = current_results.get("pass_rate", 0.0)
        previous_rate = previous_results.get("pass_rate", 0.0)
        delta = current_rate - previous_rate
        threshold = self.convergence_criteria["improvement_threshold"]

        if delta > threshold:
            return "improving"
        elif delta < -threshold:
            return "diverging"
        else:
            return "plateaued"

    def _recommend_method_switch(
        self,
        current_method: str,
        diagnosis: PoolDiagnosis,
        convergence_status: str,
    ) -> Tuple[bool, Optional[str], Optional[str]]:
        """
        Decide whether to recommend switching design methods.

        Checks method switch rules against current method, convergence status,
        and dominant failure patterns.

        Returns:
            (should_switch, recommended_method, rationale)
        """
        # Check convergence-based switches (plateaued or diverging)
        if convergence_status in ("plateaued", "diverging"):
            key = (current_method, "plateaued")
            if key in self.method_switch_rules:
                method, rationale = self.method_switch_rules[key]
                return True, method, rationale

        # Check failure-pattern-based switches
        dominant_failure = getattr(diagnosis, "dominant_failure_type", None)
        if dominant_failure:
            key = (current_method, dominant_failure)
            if key in self.method_switch_rules:
                method, rationale = self.method_switch_rules[key]
                return True, method, rationale

        # Check target-specific switches (e.g., edge strand)
        target_features = getattr(diagnosis, "target_features", [])
        if isinstance(target_features, list):
            for feature in target_features:
                key = (current_method, feature)
                if key in self.method_switch_rules:
                    method, rationale = self.method_switch_rules[key]
                    return True, method, rationale

        return False, None, None

    def _plan_parameter_adjustments(
        self,
        diagnosis: PoolDiagnosis,
        current_params: Dict,
    ) -> List[ParameterAdjustment]:
        """
        Generate parameter adjustment recommendations based on diagnosis.

        Maps the dominant failure type and sub-pattern from the diagnosis
        to concrete parameter changes using ADJUSTMENT_RULES.
        """
        adjustments = []

        dominant_failure = getattr(diagnosis, "dominant_failure_type", None)
        sub_pattern = getattr(diagnosis, "failure_sub_pattern", None)

        if not dominant_failure:
            return adjustments

        # Look up rules for this failure pattern
        failure_rules = self.adjustment_rules.get(dominant_failure, {})
        if sub_pattern and sub_pattern in failure_rules:
            applicable_rules = failure_rules[sub_pattern]
        else:
            # If no specific sub-pattern match, apply rules from the first
            # matching sub-pattern as a fallback
            applicable_rules = []
            for _sub, rules in failure_rules.items():
                applicable_rules = rules
                break

        for param, recommendation, rationale in applicable_rules:
            current_value = str(current_params.get(param, "not set"))
            recommended_value = self._compute_recommended_value(
                param, current_value, recommendation
            )
            expected_impact = self._estimate_parameter_impact(
                param, dominant_failure, sub_pattern
            )
            adjustments.append(ParameterAdjustment(
                parameter=param,
                current_value=current_value,
                recommended_value=recommended_value,
                rationale=rationale,
                expected_impact=expected_impact,
            ))

        return adjustments

    def _plan_hotspot_surgery(
        self,
        diagnosis: PoolDiagnosis,
        current_hotspots: List,
    ) -> List[str]:
        """
        Plan hotspot residue additions/removals based on diagnosis.

        Examines interface engagement patterns from the diagnosis to
        recommend adding under-engaged residues or removing over-constrained ones.
        """
        changes = []

        # Check for hotspot-related issues in the diagnosis
        sub_pattern = getattr(diagnosis, "failure_sub_pattern", None)
        hotspot_analysis = getattr(diagnosis, "hotspot_analysis", {})

        if sub_pattern == "poor_hotspot_engagement":
            # Recommend adding residues that are surface-exposed near the
            # binding epitope but not currently in the hotspot list
            suggested_additions = hotspot_analysis.get("suggested_additions", [])
            for res in suggested_additions:
                changes.append(f"add res {res}")

        elif sub_pattern == "hotspot_over_constraint":
            # Recommend removing the least-engaged hotspot residues
            least_engaged = hotspot_analysis.get("least_engaged", [])
            for res in least_engaged:
                changes.append(f"remove res {res}")

        elif sub_pattern == "no_binding_predicted":
            # Recommend a complete hotspot redesign
            if current_hotspots:
                changes.append(
                    "redesign hotspot selection — current set ineffective"
                )
                suggested = hotspot_analysis.get("alternative_hotspots", [])
                for res in suggested:
                    changes.append(f"add res {res}")

        # If diagnosis has explicit hotspot recommendations, include them
        explicit_changes = getattr(diagnosis, "hotspot_recommendations", [])
        if isinstance(explicit_changes, list):
            for change in explicit_changes:
                if change not in changes:
                    changes.append(change)

        return changes

    def _estimate_designs_needed(
        self,
        current_pass_rate: float,
        target_gold: int,
    ) -> int:
        """
        Estimate the number of designs to generate in the next round.

        Uses the current pass rate to back-calculate how many total designs
        are needed to achieve the target number of gold-tier hits, with a
        safety margin. Assumes roughly 1 in 3 passing designs reaches gold tier.
        """
        if current_pass_rate <= 0:
            # No passing designs yet — start with a large batch
            return 1000

        # Approximate: gold ≈ pass_rate * gold_fraction * total
        # gold_fraction ~ 0.3 (roughly 30% of passing designs are gold-tier)
        gold_fraction = 0.30
        effective_rate = current_pass_rate * gold_fraction

        if effective_rate <= 0:
            return 1000

        # Designs needed = target_gold / effective_rate, with 2X safety margin
        raw_estimate = target_gold / effective_rate
        designs_needed = int(math.ceil(raw_estimate * 2.0))

        # Clamp to reasonable range
        designs_needed = max(100, min(designs_needed, 10000))

        return designs_needed

    # ── Helper methods ───────────────────────────────────────────────────────

    def _extract_results(self, scored_designs: List[ScoredDesign]) -> Dict:
        """Extract summary results from a list of scored designs."""
        total = len(scored_designs)
        if total == 0:
            return {
                "total": 0, "passed": 0, "gold": 0, "silver": 0,
                "bronze": 0, "pass_rate": 0.0, "best_score": 0.0,
            }

        tiers = {"gold": 0, "silver": 0, "bronze": 0}
        best_score = 0.0
        passed = 0

        for design in scored_designs:
            tier = getattr(design, "tier", None)
            score = getattr(design, "total_score", 0.0)
            is_pass = getattr(design, "passes", False)

            if tier and tier in tiers:
                tiers[tier] += 1
            if is_pass:
                passed += 1
            if score > best_score:
                best_score = score

        pass_rate = passed / total if total > 0 else 0.0

        return {
            "total": total,
            "passed": passed,
            "gold": tiers["gold"],
            "silver": tiers["silver"],
            "bronze": tiers["bronze"],
            "pass_rate": pass_rate,
            "best_score": best_score,
        }

    def _check_convergence_criteria(
        self,
        current_results: Dict,
        round_number: int,
        history: List[RoundHistory],
    ) -> Tuple[List[str], List[str]]:
        """Check which convergence criteria are met and which remain."""
        met = []
        remaining = []

        # Gold tier count
        gold = current_results.get("gold", 0)
        min_gold = self.convergence_criteria["min_gold_tier"]
        if gold >= min_gold:
            met.append(f"Gold-tier designs: {gold} >= {min_gold}")
        else:
            remaining.append(f"Gold-tier designs: {gold} < {min_gold}")

        # Pass rate
        pass_rate = current_results.get("pass_rate", 0.0)
        min_rate = self.convergence_criteria["min_pass_rate"]
        if pass_rate >= min_rate:
            met.append(f"Pass rate: {pass_rate:.1%} >= {min_rate:.0%}")
        else:
            remaining.append(f"Pass rate: {pass_rate:.1%} < {min_rate:.0%}")

        # Max rounds
        max_rounds = self.convergence_criteria["max_rounds"]
        if round_number >= max_rounds:
            met.append(f"Max rounds reached: {round_number} >= {max_rounds}")
        else:
            remaining.append(
                f"Rounds remaining: {max_rounds - round_number}"
            )

        return met, remaining

    def _should_continue(
        self,
        round_number: int,
        convergence_status: str,
        criteria_met: List[str],
        criteria_remaining: List[str],
    ) -> bool:
        """Decide whether to continue optimization."""
        max_rounds = self.convergence_criteria["max_rounds"]

        # Stop if max rounds reached
        if round_number >= max_rounds:
            return False

        # Stop if all quality criteria are met (gold tier + pass rate)
        quality_met = sum(
            1 for c in criteria_met
            if "Gold-tier" in c or "Pass rate" in c
        )
        if quality_met >= 2:
            return False

        # Stop if diverging for two consecutive rounds
        if convergence_status == "diverging":
            # Could be a temporary dip — continue unless consistently diverging
            # This is a simple heuristic; more sophisticated logic could check
            # multiple rounds of history
            return True

        return True

    def _compute_recommended_value(
        self,
        param: str,
        current_value: str,
        recommendation: str,
    ) -> str:
        """
        Compute a concrete recommended value from the adjustment rule.

        Handles patterns like "increase 2X", "decrease to 0.05", etc.
        """
        if recommendation.startswith("increase 2X"):
            try:
                val = float(current_value)
                return str(int(val * 2) if val == int(val) else val * 2)
            except (ValueError, TypeError):
                return f"{current_value} * 2"

        if recommendation.startswith("decrease to "):
            target = recommendation.replace("decrease to ", "")
            return target

        if recommendation.startswith("decrease by "):
            delta = recommendation.replace("decrease by ", "")
            try:
                val = float(current_value)
                d = float(delta)
                result = val - d
                return str(result)
            except (ValueError, TypeError):
                return f"{current_value} - {delta}"

        if recommendation.startswith("increase by "):
            delta = recommendation.replace("increase by ", "")
            try:
                parts = delta.split("-")
                if len(parts) == 2:
                    low, high = float(parts[0]), float(parts[1])
                    mid = (low + high) / 2
                    val = float(current_value)
                    return str(int(val + mid) if val == int(val) else val + mid)
                val = float(current_value)
                d = float(delta.split()[0])
                return str(val + d)
            except (ValueError, TypeError):
                return f"{current_value} + {delta}"

        if recommendation.startswith("reduce by "):
            delta = recommendation.replace("reduce by ", "")
            try:
                parts = delta.split("-")
                if len(parts) == 2:
                    low, high = float(parts[0]), float(parts[1])
                    mid = (low + high) / 2
                    val = float(current_value)
                    return str(int(val - mid) if val == int(val) else val - mid)
                val = float(current_value)
                d = float(delta.split()[0])
                return str(val - d)
            except (ValueError, TypeError):
                return f"{current_value} - {delta}"

        # For non-numeric recommendations, return the recommendation as-is
        return recommendation

    def _estimate_parameter_impact(
        self,
        param: str,
        dominant_failure: str,
        sub_pattern: Optional[str],
    ) -> str:
        """Estimate the expected impact of a parameter adjustment."""
        impact_estimates = {
            ("num_seqs", "type_i_dominant"): "Reduce Type I failures by ~30%",
            ("mpnn_temperature", "type_i_dominant"): "Improve fold success by ~20%",
            ("backbone_noise", "type_i_dominant"): "Reduce misfolding by ~15%",
            ("contig_topology", "type_i_dominant"): "Simplify topology, improve folding ~25%",
            ("sequence_constraints", "type_i_dominant"): "Stabilize key positions, ~20% fewer failures",
            ("disulfide_bonds", "type_i_dominant"): "Add stability, ~10-15% fewer misfolded",
            ("af2_recycling", "type_i_dominant"): "Better validation, filter earlier",
            ("hotspot_residues", "type_ii_dominant"): "Improve interface quality ~30%",
            ("interface_design", "type_ii_dominant"): "Increase binding contacts ~25%",
            ("method", "type_ii_dominant"): "Method switch may improve success rate 2-10X",
            ("binder_length", "type_ii_dominant"): "Larger interface area, ~20% more contacts",
            ("contig_range", "type_ii_dominant"): "More scaffold diversity, ~15% improvement",
            ("hotspot_count", "type_ii_dominant"): "Better epitope coverage, ~20% improvement",
            ("hotspot_distribution", "type_ii_dominant"): "Better engagement distribution, ~15%",
            ("hotspot_accessibility", "type_ii_dominant"): "Ensure productive contacts, ~10%",
            ("hotspot_weight", "type_ii_dominant"): "More design freedom, ~15% improvement",
            ("interface_optimization", "type_ii_dominant"): "Stronger hydrophobic core, ~20%",
            ("affinity_maturation", "type_ii_dominant"): "Sequence-level optimization, ~25%",
            ("scoring_function", "type_ii_dominant"): "Better geometric selection, ~10%",
        }

        key = (param, dominant_failure)
        return impact_estimates.get(key, "Expected incremental improvement")

    def _estimate_improvement(
        self,
        convergence_status: str,
        adjustments: List[ParameterAdjustment],
        method_switch: bool,
    ) -> str:
        """Estimate overall improvement from planned changes."""
        if method_switch:
            return "Potentially 2-10X improvement with method switch"
        if convergence_status == "diverging":
            return "Stabilize results; aim to recover previous pass rate"
        if not adjustments:
            return "Minimal changes; monitor for convergence"

        num_adj = len(adjustments)
        if num_adj >= 3:
            return f"~30-50% improvement from {num_adj} parameter adjustments"
        elif num_adj >= 1:
            return f"~15-30% improvement from {num_adj} parameter adjustment(s)"
        return "Incremental improvement expected"


# ── History I/O ──────────────────────────────────────────────────────────────

def save_history(history_list: List[RoundHistory], path: str) -> None:
    """Save optimization history to a JSON file."""
    data = [h.to_dict() for h in history_list]
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def load_history(path: str) -> List[RoundHistory]:
    """Load optimization history from a JSON file."""
    with open(path, "r") as f:
        data = json.load(f)
    return [RoundHistory.from_dict(d) for d in data]


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Phase 4: Iterative optimization planner for protein binder design"
    )
    parser.add_argument(
        "diagnosis",
        help="Path to pool diagnosis JSON file (from diagnosis_engine.py)",
    )
    parser.add_argument(
        "scores",
        help="Path to scored designs JSON file (from design_scorer.py)",
    )
    parser.add_argument(
        "--round", type=int, required=True,
        help="Current optimization round number (1-indexed)",
    )
    parser.add_argument(
        "--params",
        help="Path to current parameters JSON file",
    )
    parser.add_argument(
        "--history",
        help="Path to optimization history JSON file",
    )
    parser.add_argument(
        "--output", default="plan.json",
        help="Output path for the optimization plan (default: plan.json)",
    )
    args = parser.parse_args()

    # Load diagnosis
    with open(args.diagnosis, "r") as f:
        diag_data = json.load(f)
    # Reconstruct PoolDiagnosis from dict — use a simple namespace for
    # attribute access if PoolDiagnosis cannot be constructed from a dict
    try:
        pool_diagnosis = PoolDiagnosis(**diag_data)
    except TypeError:
        # Fallback: create a namespace object with dict attributes
        class _Namespace:
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)
        pool_diagnosis = _Namespace(**diag_data)

    # Load scored designs
    with open(args.scores, "r") as f:
        scores_data = json.load(f)
    scored_designs = []
    for sd in scores_data:
        try:
            scored_designs.append(ScoredDesign(**sd))
        except TypeError:
            class _Namespace:
                def __init__(self, **kwargs):
                    self.__dict__.update(kwargs)
            scored_designs.append(_Namespace(**sd))

    # Load current params
    current_params = {}
    if args.params:
        with open(args.params, "r") as f:
            current_params = json.load(f)

    # Load history
    history = []
    if args.history and os.path.exists(args.history):
        history = load_history(args.history)

    # Plan
    planner = OptimizationPlanner()
    plan = planner.plan_next_round(
        round_number=args.round,
        pool_diagnosis=pool_diagnosis,
        scored_designs=scored_designs,
        current_params=current_params,
        history=history,
    )

    # Output
    with open(args.output, "w") as f:
        json.dump(plan.to_dict(), f, indent=2)

    # Print human-readable report
    print(plan.to_report())
    print(f"\nPlan written to: {args.output}")


if __name__ == "__main__":
    main()
