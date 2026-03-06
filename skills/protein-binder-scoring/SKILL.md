---
name: protein-binder-scoring
description: Composite scoring engine for protein binder designs. Compute weighted multi-metric scores combining ipSAE, ipTM, pLDDT, PAE, interface area, shape complementarity, and hotspot contact rate. Assign quality tiers (gold/silver/bronze). Rank candidates for experimental testing with diversity selection. Use after filtering to prioritize designs.
---

# Protein Binder Composite Scoring

Score, rank, and tier-assign binder designs for experimental prioritization.

## When to Use

- After filtering designs (Phase 2), to rank survivors
- To decide which designs to send for experimental testing
- To compare designs across different methods or rounds
- To assess overall campaign quality

## Usage

```bash
python skills/protein-binder-design/scripts/design_scorer.py <filtered.csv> \
    --output ranked.csv \
    --weights custom_weights.json  # optional
```

## Scoring System

### Weights (Expert-Derived)

Based on meta-analysis of 3,766 experimentally characterized binders:

| Metric | Weight | Rationale |
|--------|--------|-----------|
| ipSAE | 0.30 | Best single predictor of binding (1.4X better than ipAE) |
| ipTM | 0.25 | Strong interface confidence signal |
| pLDDT | 0.15 | Structural confidence (catches misfolding) |
| PAE | 0.10 | Complementary to ipTM |
| Interface area | 0.10 | Geometric constraint (1000-1600 Å² ideal) |
| Shape complementarity | 0.05 | Packing quality |
| Hotspot contact rate | 0.05 | Design-specific engagement |

### Quality Tiers

| Tier | Score | Expected Outcome |
|------|-------|-----------------|
| Gold | ≥ 0.80 | High confidence — prioritize for testing |
| Silver | ≥ 0.65 | Good candidates — test if budget allows |
| Bronze | ≥ 0.50 | Marginal — include for diversity only |
| Reject | < 0.50 | Do not test |

### Experimental Selection Strategy

1. Rank by composite score (primary criterion)
2. Cluster by sequence similarity (70% identity cutoff)
3. Pick top design from each cluster (ensure diversity)
4. Include 2-3 negative controls
5. Exclude designs with sequence liabilities (NG deamidation, polybasic clusters)

## Custom Weights

Create a JSON file to override default weights:

```json
{
    "ipsae": 0.30,
    "iptm": 0.25,
    "plddt": 0.15,
    "pae": 0.10,
    "interface_area": 0.10,
    "shape_complementarity": 0.05,
    "hotspot_contact": 0.05
}
```
