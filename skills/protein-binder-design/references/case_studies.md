# Case Studies — Protein Binder Design Expert System

## Case 1: Helical Epitope — EGFR Binder Design

### Target
- **Protein**: EGFR (Epidermal Growth Factor Receptor)
- **PDB**: 7KZB
- **Epitope**: Domain III, residues 340-370 (helical groove)
- **Difficulty**: Medium

### Round 1: RFdiffusion Pipeline

**Target Analysis**:
- Epitope classification: **helix** (65% helical, 20% coil, 15% strand)
- Hydrophilicity: 0.42 (moderate)
- Curvature: concave (groove between two helices)
- Recommended method: RFdiffusion (good for helical epitopes)

**Parameters**:
```
hotspot_res: [A345, A348, A352]  (3 hotspots, ~3-4 residue spacing)
contig: B1-80
num_designs: 1000
noise_scale: 1.0
ProteinMPNN: num_seqs=100, temp=0.1
```

**Results**:
- Total designs: 1000
- Pass filter: 87 (8.7%)
- Gold: 2, Silver: 12, Bronze: 73
- Dominant failure: Type I partial misfolding (48%)
- Second failure: Type II weak binding (28%)

**Diagnosis**:
- High Type I rate suggests backbone topology too complex for some designs
- Reasonable pass rate (>5%) — worth iterating

### Round 2: Optimized Parameters

**Adjustments** (based on diagnosis):
- Increased ProteinMPNN sampling: 100 → 500 (address Type I)
- Lowered temperature: 0.1 → 0.05
- Added 2 hotspots: [A345, A348, A352, A356, A360] (address weak binding)
- Increased designs: 1000 → 2000

**Results**:
- Total designs: 2000
- Pass filter: 234 (11.7%)
- Gold: 8, Silver: 31, Bronze: 195
- Improvement: +3.0% pass rate, Gold 2 → 8

**Outcome**: Converged. 8 Gold-tier designs selected for experimental testing.
Top design: composite score 0.87, ipTM 0.89, pLDDT 0.92, interface area 1340 Å².

---

## Case 2: Hydrophilic Edge-Strand — PD-L1 Binder

### Target
- **Protein**: PD-L1
- **Epitope**: Beta-sheet face, residues 54-72
- **Difficulty**: Hard (flat, hydrophilic)

### Round 1: Standard RFdiffusion (Failed)

**Target Analysis**:
- Epitope classification: **edge_strand** (55% strand, 30% coil)
- Hydrophilicity: 0.68 (high — many polar/charged residues)
- Warning: "Hydrophilic edge-strand epitope. Standard RFdiffusion expected to underperform. Recommend beta-pairing RFdiffusion."

**Result with standard RFdiffusion**:
- 1000 designs → 9 passed (0.9%)
- Gold: 0, Silver: 1, Bronze: 8
- Dominant failure: Type II no binding (72%)

### Round 2: Beta-pairing RFdiffusion

**Method switch**: Standard RFdiffusion → Beta-pairing RFdiffusion
(Rationale: edge_strand + hydrophilic → 9.2X expected improvement)

**Parameters**:
```
beta_pairing: enabled
hotspot_res: [A56, A58, A60, A64, A68]  (alternating i, i+2 pattern)
contig: B1-100
num_designs: 2000
```

**Results**:
- 2000 designs → 184 passed (9.2%)
- Gold: 5, Silver: 22, Bronze: 157
- Major improvement from 0.9% → 9.2%

**Outcome**: Converged in 2 rounds (with method switch). Beta-pairing was critical for this target.

---

## Case 3: Loop Epitope — BindCraft One-Shot

### Target
- **Protein**: IL-7Rα
- **Epitope**: Loop region, residues 85-98
- **Difficulty**: Medium-Easy

### Single Round: BindCraft

**Target Analysis**:
- Epitope classification: **loop** (75% coil, 15% helix, 10% strand)
- Hydrophilicity: 0.52 (moderate)
- Recommended method: BindCraft (best for loops, one-shot design)

**Parameters**:
```
protocol: default
target_chain: A
binder_length: 70
num_designs: 100  (BindCraft needs fewer due to higher success rate)
```

**Results**:
- 100 designs → 43 passed (43%)
- Gold: 12, Silver: 18, Bronze: 13
- No dominant failure mode (diverse failures in remaining 57)

**Outcome**: Converged in 1 round. BindCraft's one-shot approach was ideal.
Total time: ~8 hours (vs 2-3 days for RFdiffusion pipeline).

---

## Case 4: Difficult Flat Target — Multi-Round Optimization

### Target
- **Protein**: TNF-α (flat trimer interface)
- **Epitope**: Flat surface, residues 70-110
- **Difficulty**: Hard

### Round 1: RFdiffusion

**Target Analysis**:
- Epitope classification: **flat** (40% strand, 35% coil, 25% helix)
- Hydrophilicity: 0.45
- Curvature: flat
- Warning: "Flat epitope with no obvious pocket. Expected low success rate. Consider multi-method approach."

**Results**:
- 2000 designs → 32 passed (1.6%)
- Gold: 0, Silver: 3, Bronze: 29
- Dominant failure: Type II no binding (65%)

### Round 2: Parameter Optimization

**Adjustments**:
- More hotspots: 3 → 6 (distributed across flat surface)
- Larger binder: contig B1-80 → B1-120
- More designs: 2000 → 5000

**Results**:
- 5000 designs → 135 passed (2.7%)
- Gold: 1, Silver: 8, Bronze: 126
- Improvement: +1.1% (marginal)

### Round 3: Method Switch to PXDesign

**Diagnosis**: Plateaued (<5% improvement). Switch method.

**Results with PXDesign**:
- 1000 designs → 58 passed (5.8%)
- Gold: 3, Silver: 11, Bronze: 44

### Round 4: BindCraft Refinement

Used top PXDesign backbones as starting points for BindCraft refinement.

**Results**:
- 200 refined designs → 34 passed (17%)
- Gold: 7, Silver: 15, Bronze: 12

**Outcome**: Converged after 4 rounds with 2 method switches.
Key lesson: Difficult targets benefit from method diversity and multi-round iteration.

---

## Key Lessons Across Cases

1. **Match method to epitope type**: Edge-strand → beta-pairing; loop → BindCraft; helix → RFdiffusion
2. **Don't waste rounds on a bad method**: If pass rate <2% after Round 1, switch methods
3. **BindCraft needs fewer designs**: 100-200 vs 1000-5000 for RFdiffusion
4. **Hotspot count matters**: Too few → poor engagement; too many → over-constraint
5. **Method combinations work**: Use one method's outputs as starting points for another
6. **Convergence is usually 1-3 rounds** for easy/medium targets, 3-5 for hard targets
