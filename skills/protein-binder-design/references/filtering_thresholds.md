# Filtering Thresholds — Expert Reference

## Overview

Filtering thresholds for protein binder design candidates, derived from analysis of 3,766 experimentally characterized binders. Two preset tiers are provided:
- **Standard**: Balanced sensitivity/specificity, retains more candidates
- **Stringent**: High-confidence candidates only

## Metric Hierarchy (by predictive power)

1. **ipSAE** (best single predictor) — 1.4X better average precision than ipAE
2. **ipTM** (interface confidence) — widely used but has target-dependent gray zone
3. **pLDDT** (structural confidence) — good for Type I failure detection
4. **PAE** (alignment error) — complementary to ipTM
5. **Interface area** — geometric constraint
6. **Hotspot contact rate** — design-specific quality metric

## Threshold Tables

### Standard Tier

| Metric | Pass | Warning Zone | Fail | Direction |
|--------|------|-------------|------|-----------|
| ipSAE rank | Top 600 | — | > 600 | Lower rank = better |
| ipTM | ≥ 0.80 | 0.60 – 0.80 | < 0.60 | Higher = better |
| pLDDT | ≥ 0.85 | 0.70 – 0.85 | < 0.70 | Higher = better |
| PAE | ≤ 12 Å | 12 – 18 Å | > 18 Å | Lower = better |
| Interface area | 1000–1600 Å² | 850–1000 Å² | < 850 Å² | In range = best |
| Cα RMSD | ≤ 2.0 Å | 2.0 – 3.0 Å | > 3.0 Å | Lower = better |
| Hotspot contact | 0.25 – 0.50 | 0.20 – 0.25 | < 0.20 | In range = best |
| Binder size | ≤ 150 aa | 150 – 200 aa | > 200 aa | Smaller = better |

### Stringent Tier

| Metric | Pass | Warning Zone | Fail |
|--------|------|-------------|------|
| ipSAE rank | Top 300 | — | > 300 |
| ipTM | ≥ 0.85 | 0.70 – 0.85 | < 0.70 |
| pLDDT | ≥ 0.90 | 0.80 – 0.90 | < 0.80 |
| PAE | ≤ 10 Å | 10 – 14 Å | > 14 Å |
| Interface area | 1100–1500 Å² | 950–1100 Å² | < 950 Å² |
| Cα RMSD | ≤ 1.5 Å | 1.5 – 2.5 Å | > 2.5 Å |
| Hotspot contact | 0.30 – 0.50 | 0.25 – 0.30 | < 0.25 |
| Binder size | ≤ 120 aa | 120 – 150 aa | > 150 aa |

## Filter Application Order

Apply filters in this order for maximum efficiency:

```
Step 1: ipSAE top-N selection
  └── Removes ~40-70% of pool (cheapest, most discriminative)

Step 2: Structural quality (pLDDT + Cα RMSD)
  └── Catches Type I failures (misfolding)

Step 3: Interface quality (ipTM + PAE + interface area)
  └── Catches Type II failures (non-binding)

Step 4: Hotspot engagement
  └── Design-specific quality check

Step 5: Size filter
  └── Shorter binders correlate with higher experimental success
```

## Target-Specific Calibration

**Critical**: ipTM thresholds should be calibrated per target:

| Target Type | ipTM Cutoff | Rationale |
|-------------|-------------|-----------|
| Easy (helical groove) | 0.80 standard | High-confidence predictions possible |
| Medium (mixed epitope) | 0.70 | Lower cutoff to retain diversity |
| Hard (flat surface) | 0.60 | Very few designs exceed 0.70 |
| Edge-strand | 0.65 | Beta-pairing designs have different score distributions |

**When to adjust thresholds**:
- If pass rate < 1%: Relax thresholds one step (stringent → standard)
- If pass rate > 30%: Tighten thresholds (standard → stringent)
- If >50% of passes are in gray zone: Target may need method change

## ipSAE vs ipTM: When to Use Which

- **ipSAE preferred**: When you have AlphaFold2/3 PAE matrices available
- **ipTM sufficient**: When only structure prediction confidence is available
- **Both together**: Best discrimination (ROC AUC ~0.75+)
- **ipTM unreliable**: When full-length sequences include unrelated domains (truncate first!)

## Compound Scoring

Individual metrics are weak predictors (ROC AUC 0.64-0.66). Combining improves discrimination:

```
Composite AUC ≈ 0.75+ (vs 0.64-0.66 individual)
```

Use the design_scorer.py for weighted composite scoring with expert-derived weights.
