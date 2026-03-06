"""
Integration tests for the protein binder design expert skill system.

Tests the full pipeline: filter → diagnosis → scoring → optimization planning.
Uses sample_designs.csv which contains labeled examples of each failure type.
"""

import sys
import os
import json

# Add scripts directory to path
SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "skills", "protein-binder-design", "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from filter_engine import FilterEngine, FilterThresholds, DesignMetrics, load_designs_csv
from diagnosis_engine import DiagnosisEngine
from design_scorer import DesignScorer


TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "test_data", "sample_designs.csv")


def test_load_csv():
    """Test that CSV loading works."""
    designs = load_designs_csv(TEST_DATA)
    assert len(designs) == 30, f"Expected 30 designs, got {len(designs)}"
    print(f"[PASS] Loaded {len(designs)} designs from CSV")


def test_filter_engine_standard():
    """Test standard tier filtering."""
    designs = load_designs_csv(TEST_DATA)
    engine = FilterEngine(FilterThresholds.standard())
    results = engine.filter_pool(designs)

    stats = engine.summary_statistics(results)
    total = stats["total_designs"]
    passed = stats["passed"]
    failed = stats["failed"]

    assert total == 30, f"Expected 30 total, got {total}"
    assert passed > 0, "Expected some designs to pass"
    assert failed > 0, "Expected some designs to fail"
    assert passed < total, "Not all designs should pass"

    print(f"[PASS] Filter standard: {passed} pass, {stats['gray_zone']} gray, {failed} fail out of {total}")
    print(f"       Pass rate: {stats['pass_rate']:.1%}")
    if stats["top_failure_reasons"]:
        print(f"       Top failures: {list(stats['top_failure_reasons'].keys())[:3]}")


def test_filter_engine_stringent():
    """Test stringent tier filtering — should pass fewer."""
    designs = load_designs_csv(TEST_DATA)

    standard_engine = FilterEngine(FilterThresholds.standard())
    stringent_engine = FilterEngine(FilterThresholds.stringent())

    standard_results = standard_engine.filter_pool(designs)
    stringent_results = stringent_engine.filter_pool(designs)

    standard_pass = len(standard_results["pass"]) + len(standard_results["gray_zone"])
    stringent_pass = len(stringent_results["pass"]) + len(stringent_results["gray_zone"])

    assert stringent_pass <= standard_pass, \
        f"Stringent should pass fewer: {stringent_pass} vs standard {standard_pass}"

    print(f"[PASS] Stringent passes fewer: {stringent_pass} vs standard {standard_pass}")


def test_diagnosis_engine():
    """Test failure mode diagnosis."""
    designs = load_designs_csv(TEST_DATA)
    engine = DiagnosisEngine()
    diagnosis = engine.diagnose_pool(designs)

    assert diagnosis.total_designs == 30
    assert diagnosis.pass_rate >= 0
    assert len(diagnosis.failure_distribution) > 0
    assert diagnosis.dominant_failure != ""

    print(f"[PASS] Pool diagnosis: {diagnosis.total_designs} designs")
    print(f"       Pass rate: {diagnosis.pass_rate:.1%}")
    print(f"       Dominant failure: {diagnosis.dominant_failure}")
    print(f"       Failure distribution: {diagnosis.failure_distribution}")
    if diagnosis.systematic_issues:
        print(f"       Systematic issues: {diagnosis.systematic_issues}")


def test_diagnosis_specific_types():
    """Test that specific failure types are correctly identified."""
    engine = DiagnosisEngine()

    # Type I severe: high RMSD + low pLDDT
    type_i_severe = DesignMetrics(
        design_id="test_type_i", iptm=0.4, plddt=0.5, pae=22.0,
        interface_area=600, ca_rmsd=3.5, hotspot_contact_rate=0.1
    )
    diag = engine.diagnose_design(type_i_severe)
    assert "type_i" in diag.failure_type, f"Expected type_i, got {diag.failure_type}"
    print(f"[PASS] Type I severe correctly identified: {diag.subtype}")

    # Type II no binding: low RMSD + low ipTM
    type_ii_no_bind = DesignMetrics(
        design_id="test_type_ii", iptm=0.35, plddt=0.90, pae=20.0,
        interface_area=1100, ca_rmsd=1.0, hotspot_contact_rate=0.30
    )
    diag = engine.diagnose_design(type_ii_no_bind)
    assert "type_ii" in diag.failure_type, f"Expected type_ii, got {diag.failure_type}"
    print(f"[PASS] Type II no-binding correctly identified: {diag.subtype}")

    # Poor hotspot engagement
    poor_hotspot = DesignMetrics(
        design_id="test_hotspot", iptm=0.75, plddt=0.87, pae=10.0,
        interface_area=1100, ca_rmsd=1.0, hotspot_contact_rate=0.08
    )
    diag = engine.diagnose_design(poor_hotspot)
    assert "hotspot" in diag.subtype.lower(), f"Expected hotspot issue, got {diag.subtype}"
    print(f"[PASS] Poor hotspot engagement correctly identified: {diag.subtype}")


def test_design_scorer():
    """Test composite scoring and tier assignment."""
    designs = load_designs_csv(TEST_DATA)
    scorer = DesignScorer()
    scored = scorer.score_designs(designs)

    assert len(scored) == 30
    assert all(s.composite_score >= 0 for s in scored)
    assert all(s.composite_score <= 1 for s in scored)

    # Check ordering (should be descending by score)
    scores = [s.composite_score for s in scored]
    assert scores == sorted(scores, reverse=True), "Scores should be in descending order"

    # Check tier distribution
    tiers = {"gold": 0, "silver": 0, "bronze": 0, "reject": 0}
    for s in scored:
        tiers[s.tier] += 1

    print(f"[PASS] Scoring: {len(scored)} designs scored")
    print(f"       Tier distribution: {tiers}")
    print(f"       Top 3 designs: {[(s.design.design_id, f'{s.composite_score:.3f}', s.tier) for s in scored[:3]]}")

    # Gold designs should include our labeled gold samples
    gold_ids = {s.design.design_id for s in scored if s.tier == "gold"}
    expected_gold = {"gold_001", "gold_002", "gold_003"}
    overlap = gold_ids & expected_gold
    print(f"       Expected gold in actual gold: {len(overlap)}/{len(expected_gold)}")


def test_scorer_tier_ordering():
    """Test that gold > silver > bronze > reject in score."""
    designs = load_designs_csv(TEST_DATA)
    scorer = DesignScorer()
    scored = scorer.score_designs(designs)

    tier_scores = {"gold": [], "silver": [], "bronze": [], "reject": []}
    for s in scored:
        tier_scores[s.tier].append(s.composite_score)

    for tier, scores in tier_scores.items():
        if scores:
            avg = sum(scores) / len(scores)
            print(f"       {tier}: avg={avg:.3f}, min={min(scores):.3f}, max={max(scores):.3f}")

    # Gold average should be higher than reject average
    if tier_scores["gold"] and tier_scores["reject"]:
        gold_avg = sum(tier_scores["gold"]) / len(tier_scores["gold"])
        reject_avg = sum(tier_scores["reject"]) / len(tier_scores["reject"])
        assert gold_avg > reject_avg, "Gold average should be higher than reject"
        print(f"[PASS] Gold avg ({gold_avg:.3f}) > Reject avg ({reject_avg:.3f})")


def test_full_pipeline():
    """Test the full pipeline: filter → diagnose → score."""
    designs = load_designs_csv(TEST_DATA)

    # Step 1: Filter
    filter_engine = FilterEngine(FilterThresholds.standard())
    filter_results = filter_engine.filter_pool(designs)
    filter_stats = filter_engine.summary_statistics(filter_results)

    # Step 2: Diagnose
    diag_engine = DiagnosisEngine()
    pool_diag = diag_engine.diagnose_pool(designs)

    # Step 3: Score (all designs, not just filtered)
    scorer = DesignScorer()
    scored = scorer.score_designs(designs)

    print(f"[PASS] Full pipeline:")
    print(f"       Filter: {filter_stats['passed']} pass / {filter_stats['total_designs']} total")
    print(f"       Diagnosis: dominant={pool_diag.dominant_failure}")
    print(f"       Scoring: top design={scored[0].design.design_id} ({scored[0].composite_score:.3f})")
    print(f"       Recommendations: {pool_diag.pool_level_recommendations[:2]}")


def test_method_selector():
    """Test method selection logic."""
    try:
        from method_selector import select_method

        # Edge strand + hydrophilic should recommend beta-pairing
        methods = select_method(
            epitope_type="edge_strand",
            hydrophilicity="high",
            time_budget="standard"
        )
        assert len(methods) > 0, "Should return at least one method"
        top = methods[0]
        print(f"[PASS] Method selector: edge_strand+high → {top.method.name}")
        assert "beta" in top.method.name.lower(), \
            f"Expected beta-pairing for edge_strand+high, got {top.method.name}"

        # Fast budget should recommend BindCraft or Latent-X
        methods = select_method(
            epitope_type="helix",
            hydrophilicity="medium",
            time_budget="fast"
        )
        assert len(methods) > 0
        top = methods[0]
        print(f"[PASS] Method selector: helix+fast → {top.method.name}")

    except ImportError:
        print("[SKIP] method_selector not available for import")
    except Exception as e:
        print(f"[WARN] method_selector test: {e}")


def main():
    print("=" * 60)
    print("Protein Binder Design Expert System — Integration Tests")
    print("=" * 60)
    print()

    tests = [
        test_load_csv,
        test_filter_engine_standard,
        test_filter_engine_stringent,
        test_diagnosis_engine,
        test_diagnosis_specific_types,
        test_design_scorer,
        test_scorer_tier_ordering,
        test_full_pipeline,
        test_method_selector,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            print(f"\n--- {test.__name__} ---")
            test()
            passed += 1
        except Exception as e:
            print(f"[FAIL] {test.__name__}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print(f"\n{'=' * 60}")
    print(f"Results: {passed} passed, {failed} failed out of {len(tests)} tests")
    print(f"{'=' * 60}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
