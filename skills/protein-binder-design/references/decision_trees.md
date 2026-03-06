# Decision Trees — Protein Binder Design Expert System

## 1. Method Selection Decision Tree

```
START: New binder design campaign
│
├── Is this a macrocyclic binder?
│   └── YES → Latent-X (90%+ success, picomolar affinity)
│
├── Is the epitope an edge-strand with high hydrophilicity (>50% polar)?
│   └── YES → Beta-pairing RFdiffusion (9.2X improvement)
│
├── What is the time budget?
│   ├── FAST (<1 day) → BindCraft (one-shot, 10-100% success)
│   │
│   ├── STANDARD (1-3 days)
│   │   ├── Epitope type?
│   │   │   ├── Helix → RFdiffusion or BindCraft
│   │   │   ├── Loop → BindCraft (best for loops)
│   │   │   ├── Flat → RFdiffusion or PXDesign
│   │   │   ├── Concave → RFdiffusion (pocket-filling)
│   │   │   └── Mixed → PXDesign (dual modality)
│   │   │
│   │   └── Need maximum affinity?
│   │       └── YES → Latent-X
│   │
│   └── EXTENDED (days-weeks)
│       └── Multi-method parallel:
│           1. BindCraft (fast screening)
│           2. RFdiffusion (diverse backbones)
│           3. PXDesign (alternative folds)
│           → Pool and score together
```

## 2. Threshold Selection Decision Tree

```
START: Choose filtering tier
│
├── How many designs in the pool?
│   ├── < 500 → Standard tier (preserve candidates)
│   ├── 500-5000 → Standard tier
│   └── > 5000 → Stringent tier (reduce to manageable set)
│
├── Target difficulty?
│   ├── Easy (helical groove) → Stringent tier
│   ├── Medium → Standard tier
│   └── Hard (flat surface) → Relaxed thresholds
│       └── Relax ipTM to 0.50, pLDDT to 0.70
│
├── Previous round pass rate?
│   ├── > 30% → Tighten to stringent
│   ├── 5-30% → Keep current
│   ├── 1-5% → Relax one step
│   └── < 1% → Consider method switch before relaxing
│
└── How many candidates needed for experiment?
    ├── < 50 → Stringent (top quality only)
    ├── 50-200 → Standard
    └── > 200 → Standard + diversity selection
```

## 3. Method Switch Decision Tree

```
START: Current method not producing results
│
├── Convergence status?
│   ├── PLATEAUED (2+ rounds, <5% improvement)
│   │   ├── Current = RFdiffusion → Switch to BindCraft
│   │   ├── Current = BindCraft → Switch to PXDesign
│   │   ├── Current = PXDesign → Switch to BindCraft
│   │   └── Current = Beta-pairing → Switch to BindCraft
│   │
│   ├── DIVERGING (pass rate declining)
│   │   └── STOP. Re-evaluate target and epitope selection.
│   │       Consider: Is this target designable at all?
│   │
│   └── IMPROVING (keep current method)
│       └── Continue with parameter adjustments only
│
├── Dominant failure mode?
│   ├── Type I dominant (>60%)
│   │   ├── Current = BindCraft → Switch to RFdiffusion
│   │   │   (RFdiffusion backbones may be more stable)
│   │   └── Current = RFdiffusion → Fix backbone parameters first
│   │
│   └── Type II dominant (>60%)
│       ├── Edge-strand target? → Switch to Beta-pairing RFdiffusion
│       └── General target → Switch to BindCraft or PXDesign
│
└── Tried multiple methods already?
    └── YES → Consider:
        1. Change epitope entirely
        2. Try peptide binder instead of mini-protein
        3. Use experimental screening (display libraries)
```

## 4. Hotspot Count Decision Tree

```
START: How many hotspots to specify?
│
├── Epitope type?
│   ├── Helix → 3-5 hotspots
│   │   └── Pattern: Every 3-4 residues along one face
│   │       (same side of helix, i, i+3 or i+4)
│   │
│   ├── Loop → 2-4 hotspots
│   │   └── Pattern: Flanking residues (start+end of loop)
│   │       AVOID: Flexible tip residues (high B-factor)
│   │
│   ├── Flat surface → 4-6 hotspots
│   │   └── Pattern: Distributed across surface
│   │       Space at least 8 Å apart (Cβ distance)
│   │
│   ├── Edge-strand → 3-5 hotspots
│   │   └── Pattern: Alternating i, i+2 (beta-strand register)
│   │
│   └── Concave pocket → 3-4 hotspots
│       └── Pattern: 1-2 at rim + 1-2 in pocket interior
│
├── First generation or subsequent?
│   ├── First → Start with 3 hotspots (conservative)
│   │   └── Pilot run with 200 designs to test
│   └── Subsequent → Adjust based on contact rate:
│       ├── Contact rate < 20% → Add 2 more hotspots
│       ├── Contact rate 20-50% → Keep current
│       └── Contact rate > 50% → Remove 2 hotspots
│
└── Quality checks for hotspot selection:
    ├── All hotspots surface-exposed? (SASA > 20 Å²)
    ├── No hotspots in flexible loops? (B-factor < 60)
    ├── Not too clustered? (min pairwise Cβ distance > 6 Å)
    └── Not at glycosylation sites? (avoid NxS/T motifs)
```

## 5. Convergence / Termination Decision Tree

```
START: Should we continue iterating?
│
├── Check convergence criteria:
│   ├── ≥ 5 Gold-tier designs? → CONVERGED. Stop.
│   ├── Pass rate ≥ 10%? → LIKELY CONVERGED. One more round for diversity.
│   └── Neither met → Continue...
│
├── Check iteration limits:
│   ├── Round > 5? → STOP. Switch approach entirely.
│   └── Round ≤ 5 → Check improvement...
│
├── Check improvement trend:
│   ├── Pass rate improved > 5% from last round? → IMPROVING. Continue.
│   ├── Pass rate changed < 5%? → PLATEAUED.
│   │   ├── Round 2-3 → Try parameter adjustments
│   │   └── Round 4-5 → Switch method
│   └── Pass rate declined? → DIVERGING.
│       └── Revert to best parameters and reconsider strategy.
│
└── Check diminishing returns:
    ├── Gold-tier count increasing? → Continue.
    ├── Gold-tier count stable? → One more round with diversity focus.
    └── Gold-tier count decreasing? → STOP. Use best results so far.
```

## 6. Experimental Prioritization Decision Tree

```
START: Which designs to test experimentally?
│
├── Budget: How many designs can be tested?
│   ├── < 20 → Gold-tier only, diverse sequences
│   ├── 20-100 → Gold + top Silver, diversity-selected
│   └── > 100 → Gold + Silver + top Bronze
│
├── Selection criteria:
│   1. Rank by composite score (primary)
│   2. Cluster by sequence similarity (70% identity cutoff)
│   3. Pick top design from each cluster (diversity)
│   4. Ensure epitope coverage (designs targeting different subregions)
│
├── Red flags to exclude:
│   ├── Cysteines at interface (unless disulfide intended)
│   ├── Deamidation-prone NG motifs at interface
│   ├── Polybasic clusters (RKRK sequences → aggregation)
│   ├── Very hydrophobic surface (GRAVY > 0.5)
│   └── Extreme pI (< 4 or > 10)
│
└── Include controls:
    ├── 2-3 negative controls (shuffled sequences)
    ├── 1 positive control (if known binder exists)
    └── 1-2 "borderline" designs (Gray zone) for calibration
```
