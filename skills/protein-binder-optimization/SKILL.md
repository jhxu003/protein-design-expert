---
name: protein-binder-optimization
description: Plan iterative optimization rounds for protein binder design campaigns. Analyze current round results, recommend parameter adjustments (hotspot count/placement, backbone noise, sequence design temperature), decide whether to switch methods, and define convergence criteria. Use after scoring to plan the next design round.
---

# Protein Binder Iterative Optimization

Plan the next round of binder design based on current results.

## When to Use

- After scoring results from a design round
- When pass rate is below target and user wants to improve
- When deciding whether to continue current method or switch
- When adjusting hotspot residues or design parameters

## Usage

```bash
python skills/protein-binder-design/scripts/optimization_planner.py \
    <diagnosis.json> <scores.json> \
    --round <round_number> \
    --params current_params.json \
    --history history.json \
    --output plan.json
```

## Convergence Criteria

The optimization loop should STOP when ANY of these are met:

| Criterion | Threshold | Action |
|-----------|-----------|--------|
| Gold-tier designs | ≥ 5 | **Success** — proceed to experiments |
| Pass rate | ≥ 10% | One more round for diversity, then stop |
| Max rounds | 5 | Stop — switch to entirely different approach |
| No improvement | < 5% for 2 rounds | Switch method or epitope |

## Parameter Adjustment Rules

### When Type I Failures Dominate (>60%)
- Increase ProteinMPNN sampling: num_seqs × 2
- Lower ProteinMPNN temperature: → 0.05
- Reduce backbone noise: -0.3 to -0.5
- Simplify backbone topology

### When Type II Failures Dominate (>60%)
- Add 2-3 hotspot residues
- Extend binder length (+30-50 residues)
- Consider switching to BindCraft (10X higher success)

### When Plateaued
- Switch design method (see decision tree)
- Change epitope region
- Try multi-method parallel approach

## Method Switch Recommendations

| Current Method | Condition | Switch To |
|---------------|-----------|-----------|
| RFdiffusion | Plateaued | BindCraft |
| RFdiffusion | Edge-strand target | Beta-pairing RFdiffusion |
| RFdiffusion | Type II dominant | PXDesign |
| BindCraft | Plateaued | PXDesign |
| BindCraft | Type I dominant | RFdiffusion |
| PXDesign | Plateaued | BindCraft |

## History Tracking

The planner maintains a round-by-round history in JSON format, tracking:
- Method and parameters used per round
- Results (pass rate, tier distribution, best score)
- Dominant failure mode
- Convergence trend

This enables trend analysis across rounds to detect plateaus and recommend method switches.

See `references/decision_trees.md` for the complete convergence and method-switch decision trees.
