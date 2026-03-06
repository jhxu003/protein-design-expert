---
name: protein-binder-diagnosis
description: Diagnose protein binder design failures and filter design pools. Classify Type I (misfolding) and Type II (non-binding) failures, apply tiered filtering thresholds (ipSAE, ipTM, pLDDT, PAE, interface area), analyze hotspot contact rates, and recommend specific fixes. Use after generating binder candidates with RFdiffusion, BindCraft, or similar tools.
---

# Protein Binder Design Diagnosis

Diagnose why binder designs are failing and filter design pools using expert thresholds.

## When to Use

- User has generated binder designs and wants to filter them
- Designs are failing and user wants to understand why
- User needs to know which designs to prioritize for experiments
- User has low pass rates and needs actionable fixes

## Usage

### Filter a Design Pool

```bash
python skills/protein-binder-design/scripts/filter_engine.py <results.csv> \
    --tier <standard|stringent> \
    --output filtered.json \
    --summary
```

### Diagnose Failures

```bash
python skills/protein-binder-design/scripts/diagnosis_engine.py <results.csv> \
    --output diagnosis.json
```

### Expected CSV Columns

The input CSV should contain design metrics (columns are flexible — uses what's available):

| Column | Description | Required |
|--------|-------------|----------|
| design_id | Unique identifier | Yes |
| iptm | Interface predicted TM-score (0-1) | Yes |
| plddt | Mean pLDDT of binder (0-1) | Yes |
| pae | Interface PAE (Å) | Recommended |
| ca_rmsd | Cα RMSD design vs prediction (Å) | Recommended |
| interface_area | Buried surface area (Å²) | Recommended |
| hotspot_contact_rate | Fraction of hotspots contacted | Recommended |
| ipsae_min | ipSAE score | Optional (best predictor if available) |
| shape_complementarity | Shape complementarity (0-1) | Optional |
| binder_length | Number of residues | Optional |

## Failure Mode Quick Reference

| Failure | Signature | Top Fix |
|---------|-----------|---------|
| Severe misfolding | RMSD > 2, pLDDT < 0.7 | Simplify backbone topology |
| Partial misfolding | RMSD > 2, pLDDT ≥ 0.7 | Add sequence constraints |
| No binding signal | RMSD ≤ 2, ipTM < 0.6 | Switch to BindCraft |
| Interface too small | ipTM 0.6-0.8, area < 850 | Extend binder length |
| Poor hotspot engagement | contact rate < 0.2 | Add 2-3 more hotspots |
| Over-constrained | contact rate > 0.5 | Reduce hotspot count |

See `references/failure_modes.md` for the complete failure taxonomy and diagnostic flowchart.
