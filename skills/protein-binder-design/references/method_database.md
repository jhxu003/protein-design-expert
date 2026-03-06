# Method Database — Protein Binder Design Expert Knowledge

## Overview

This document catalogs the 5 major computational methods for de novo protein binder design, with expert-curated selection criteria, parameters, and expected outcomes.

---

## 1. RFdiffusion

**Type:** Diffusion-based backbone generation
**Paper:** Watson et al., Nature 2023
**Pipeline:** RFdiffusion → ProteinMPNN → AlphaFold2

### Success Rates
- Baseline: 0.4-7.5% (target-dependent)
- With optimized hotspots: up to 15%
- EGFR competition: ~14% binding (53/378 expressed)

### Best For
- General targets with well-defined binding sites
- Helical epitopes
- Motif scaffolding (grafting known binding motifs)
- When modularity and control are needed

### Key Parameters
| Parameter | Default | Expert Guidance |
|-----------|---------|-----------------|
| `hotspot_res` | Required | 3-6 residues, format: `[A30,A33,A34]` |
| `contig` | Target-specific | `B1-80` for ~80-residue binder |
| `num_designs` | 1000 | Start with 1000, increase if pass rate <5% |
| `noise_scale` | 1.0 | Reduce to 0.5-0.7 for more conservative backbones |

### Downstream Steps
1. **ProteinMPNN** for inverse folding (num_seqs=100, temperature=0.1)
2. **AlphaFold2** for structure validation (5 models, 3 recycles)
3. **ipSAE/ipTM** scoring for binding prediction

### Strengths
- Most established method with large community
- Flexible contigmap for complex topologies
- Motif scaffolding capability
- Works well for helical interfaces

### Weaknesses
- Multi-step pipeline (more failure points)
- Lower success rate than newer methods
- Poor on hydrophilic/polar targets
- Requires separate sequence design step

---

## 2. BindCraft

**Type:** Gradient-based hallucination with integrated AF2
**Paper:** Pacesa et al., Nature 2025
**Pipeline:** Self-contained (no downstream tools needed)

### Success Rates
- 10-100% experimental binding (target-dependent)
- 10X higher than RFdiffusion pipeline
- One-shot design (design + validate in single run)

### Best For
- Rapid prototyping and iteration
- Loop epitopes
- Time-constrained projects
- Beginners (simpler workflow)

### Key Parameters
| Parameter | Default | Expert Guidance |
|-----------|---------|-----------------|
| `protocol` | "default" | "fast" for screening, "slow" for quality |
| `target_chain` | Required | Chain ID of target protein |
| `binder_length` | 60-100 | Shorter = higher success, longer = more contacts |
| `num_designs` | 100 | Smaller pool needed due to higher hit rate |

### Optimization Stages (Internal)
1. Continuous logit space (75 iterations) — explore diverse amino acids
2. Softmax annealing (45 iterations) — sharpen toward discrete
3. Straight-through estimator (5 iterations) — discrete with gradient flow
4. Greedy mutation (15 iterations) — accept improving point mutations

### Strengths
- Highest success rate among general-purpose methods
- Self-contained (no ProteinMPNN/AF2 pipeline needed)
- Repredicts complex at each optimization step
- Good for diverse epitope types

### Weaknesses
- Newer, less community data and tooling
- May plateau on very difficult targets
- Less control over backbone topology
- Resource-intensive per design

---

## 3. Latent-X

**Type:** Latent diffusion in learned protein space
**Paper:** Khalil et al., 2025

### Success Rates
- 90%+ for macrocyclic binders
- 10-64% for mini-binders
- Picomolar to nanomolar affinities achieved

### Best For
- Macrocyclic binder design
- Maximum affinity requirements
- Speed-critical projects (10X faster than pipelines)

### Key Advantages
- All-atom generation (side chains included)
- Structurally diverse outputs (not biased toward alpha-helices)
- Very high success rates on tested targets

### Limitations
- Best validated for macrocycles specifically
- Less general applicability data than RFdiffusion/BindCraft
- Newer method with smaller user base

---

## 4. PXDesign

**Type:** Dual-modality (diffusion + hallucination)
**Paper:** PXDesign, 2025

### Success Rates
- 17-82% across diverse targets
- Two modes: PXDesign-d (diffusion) and PXDesign-h (hallucination)

### Best For
- Diverse targets where single methods plateau
- Mixed epitope types
- When broader fold coverage is needed

### Key Advantages
- Dual approach covers more design space
- Broader fold coverage than RFdiffusion (less α-helix bias)
- Good balance of speed and quality

### Limitations
- Variable success rate across targets
- Requires careful parameter tuning
- Still needs AF2 validation downstream

---

## 5. Beta-pairing RFdiffusion

**Type:** Modified RFdiffusion with β-sheet conditioning
**Paper:** Sahtoe et al., Nature Communications 2026

### Success Rates
- 9.2% high-quality designs (vs 0.98% unconditioned)
- ~10X improvement for edge-strand targets specifically

### Best For
- Hydrophilic polar surfaces
- Edge-strand epitopes
- Targets that fail with standard RFdiffusion

### Mechanism
- Geometrically matched extended β-sheets
- Hydrogen bonding networks at interface
- Specifically designed for polar/charged interfaces

### When to Use
```
IF epitope_type == "edge_strand" AND hydrophilicity > 0.50:
    USE beta-pairing-RFdiffusion
    EXPECTED 9.2X improvement over standard RFdiffusion
```

---

## Method Selection Quick Reference

| Criterion | Recommended Method |
|-----------|-------------------|
| Fastest turnaround | BindCraft or Latent-X |
| Highest general success rate | BindCraft (10-100%) |
| Macrocyclic binders | Latent-X (90%+) |
| Maximum affinity | Latent-X (picomolar) |
| Hydrophilic/edge-strand target | Beta-pairing RFdiffusion |
| Modularity/control needed | RFdiffusion |
| Multiple methods plateau | PXDesign (dual modality) |
| Beginner-friendly | BindCraft |
| Motif scaffolding | RFdiffusion |
