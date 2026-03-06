# Failure Modes in Protein Binder Design — Expert Taxonomy

## Overview

Protein binder design failures fall into two main categories:
- **Type I**: Sequence doesn't fold to the intended structure
- **Type II**: Monomer folds correctly but doesn't bind the target

Understanding the failure mode is critical for choosing the right fix.

---

## Type I Failures — Folding Problems

### Type I-A: Severe Misfolding
**Signature**: Cα RMSD > 2.0 Å AND pLDDT < 0.70

**Root cause**: The designed backbone geometry is not achievable by any amino acid sequence. The backbone has internal strain, impossible torsion angles, or topology that no sequence can stabilize.

**Evidence**:
- AlphaFold predicts a completely different fold
- Large Cα RMSD between design model and AF2 prediction
- Low overall pLDDT (AF2 is not confident in any structure)

**Expert fixes** (priority order):
1. Simplify backbone topology (reduce loops, use more regular secondary structure)
2. Increase ProteinMPNN sampling: num_seqs 100 → 500
3. Lower ProteinMPNN temperature: 0.1 → 0.05
4. Check for strained geometry in the designed backbone
5. If using RFdiffusion, reduce noise_scale by 0.3

### Type I-B: Partial Misfolding
**Signature**: Cα RMSD > 2.0 Å AND pLDDT ≥ 0.70

**Root cause**: The sequence folds into a stable but unintended alternative conformation. AF2 finds a lower-energy state different from the design.

**Evidence**:
- AF2 prediction is confident (high pLDDT) but wrong (high RMSD)
- The alternative fold may be a known domain or common fold

**Expert fixes** (priority order):
1. Add sequence constraints at positions critical for the intended fold
2. Use AF2 with more recycling iterations (3 → 12) for better prediction
3. Add disulfide bonds to lock the intended conformation
4. Consider the AF2-predicted fold as a new starting point
5. Run ProteinMPNN with --fixed_residues at core positions

### Type I-C: Sequence Collapse
**Signature**: Very low sequence diversity across designs (>90% identity)

**Root cause**: ProteinMPNN sampling temperature too low, or backbone too constrained, yielding only one viable sequence family.

**Expert fixes**:
1. Increase ProteinMPNN temperature: 0.05 → 0.15
2. Relax backbone constraints (allow more flexibility)
3. Try multiple independent backbone starting points

---

## Type II Failures — Binding Problems

### Type II-A: No Binding Signal
**Signature**: Cα RMSD ≤ 2.0 Å AND ipTM < 0.60

**Root cause**: The binder folds correctly but AF2 sees no evidence of complex formation. The interface is not stable or specific.

**Evidence**:
- Good monomer structure (low RMSD, high pLDDT)
- But AF2 does not predict the binder-target complex
- ipTM < 0.6 means AF2 essentially says "no interaction"

**Expert fixes** (priority order):
1. Re-evaluate hotspot selection (current hotspots may not be at a druggable epitope)
2. Increase interface area (larger binder, more contact residues)
3. Improve shape complementarity at the interface
4. **Consider switching to BindCraft** (10X higher success rate)
5. Try a completely different epitope on the target

### Type II-B: Interface Too Small
**Signature**: Cα RMSD ≤ 2.0 Å AND 0.60 ≤ ipTM < 0.80 AND interface_area < 850 Å²

**Root cause**: The binder makes some contacts but the buried surface area is insufficient for stable binding. Typical natural protein-protein interfaces are 1000-2000 Å².

**Expert fixes**:
1. Extend binder length by 30-50 residues
2. Use wider contig range (e.g., B1-80 → B1-120)
3. Add flanking helices to increase buried surface
4. Increase hotspot count by 2-3 to drive more contacts

### Type II-C: Poor Hotspot Engagement
**Signature**: hotspot_contact_rate < 0.20

**Root cause**: The binder is not making contacts with the intended hotspot residues on the target. Only ~20-50% of specified hotspots are typically contacted.

**Evidence**:
- Design geometry positions binder away from hotspots
- Hotspots may be buried or inaccessible

**Expert fixes**:
1. Increase number of hotspots by 2-3 (more chances for engagement)
2. Redistribute hotspots across the epitope (not clustered)
3. Verify hotspot accessibility (check SASA — must be surface-exposed)
4. Use within 10 Å Cβ distance criterion for hotspot placement

### Type II-D: Hotspot Over-constraint
**Signature**: hotspot_contact_rate > 0.50

**Root cause**: Too many hotspot constraints prevent the design algorithm from exploring optimal interface geometries. The binder is forced into a suboptimal configuration.

**Expert fixes**:
1. Reduce hotspot count by 2-3 (keep the most important ones)
2. Remove hotspots at the epitope periphery (keep core contacts)
3. Use softer distance-based constraints instead of hard hotspot requirements

### Type II-E: Weak Binding
**Signature**: Cα RMSD ≤ 2.0 Å AND 0.60 ≤ ipTM < 0.80 AND interface_area ≥ 850 Å²

**Root cause**: Interface geometry is suboptimal. The binder makes contacts but the packing, electrostatics, or hydrophobic complementarity is poor.

**Expert fixes**:
1. Optimize interface packing (reduce cavities)
2. Add hydrophobic contacts at the core of the interface
3. Run affinity maturation with ProteinMPNN (redesign interface residues only)
4. Check for buried unsatisfied hydrogen bond donors/acceptors

---

## Systematic Failure Patterns (Pool-Level)

When analyzing an entire design pool, look for these systematic patterns:

### >60% Type I Failures
**Diagnosis**: Backbone generation parameters are fundamentally wrong.
**Action**: Don't generate more designs with same parameters. Fix backbone first.

### >60% Type II Failures
**Diagnosis**: Interface design strategy needs complete revision.
**Action**: Re-evaluate epitope selection and hotspot placement.

### >40% Poor Hotspot Engagement
**Diagnosis**: Hotspot selection is problematic.
**Action**: Comprehensive hotspot redesign before next round.

### Pass Rate < 1%
**Diagnosis**: Current method may not be suitable for this target.
**Action**: Consider switching design method entirely.

### Pass Rate < 5% for 2+ Consecutive Rounds
**Diagnosis**: Target is difficult. Current approach has plateaued.
**Action**: Switch method OR change epitope OR both.

---

## Scoring Model Complexity (Overath et al. 2025)

A key finding from the Overath et al. 2025 meta-analysis (3,760 designs, 15 targets) is that **complex scoring models are unnecessary**. Simple logistic regression with just 2-3 features (ipSAE_min ranking + RMSD binder < 3.73 A + shape_complementarity > 0.62) generalizes across targets as well as XGBoost trained on hundreds of features. This means:

- **Overfitting risk**: Complex multi-feature scoring models (random forests, gradient boosting with many features) tend to overfit to specific targets and do not generalize
- **Actionable implication**: When building custom scoring pipelines, prefer simple filters with a small number of well-validated metrics over elaborate ML models
- **Avoid stacking weak predictors**: Adding more low-quality metrics does not improve discrimination and can hurt generalization

---

## Diagnostic Flowchart

```
Start: Design fails filter
│
├── Cα RMSD > 2.0 Å?
│   ├── YES → Type I Failure
│   │   ├── pLDDT < 0.70? → Type I-A (Severe Misfolding)
│   │   └── pLDDT ≥ 0.70? → Type I-B (Partial Misfolding)
│   │
│   └── NO → Check for Type II
│       ├── ipTM < 0.60? → Type II-A (No Binding Signal)
│       ├── interface_area < 850? → Type II-B (Interface Too Small)
│       ├── hotspot_contact < 0.20? → Type II-C (Poor Hotspot)
│       ├── hotspot_contact > 0.50? → Type II-D (Over-constraint)
│       └── else → Type II-E (Weak Binding)
```
