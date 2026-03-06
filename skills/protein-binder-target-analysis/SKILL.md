---
name: protein-binder-target-analysis
description: Analyze protein targets for binder design. Classify epitope type (helix, loop, flat, concave, edge-strand), assess druggability, compute surface properties, identify hotspot residues, and recommend design methods with expected success rates. Use before starting any binder design campaign.
---

# Protein Binder Target Analysis

Analyze a protein target to determine the optimal binder design strategy.

## When to Use

- User wants to design a binder and has a target PDB structure
- User needs to choose hotspot residues
- User needs method recommendation for a specific target
- User wants to assess target difficulty before committing resources

## Usage

### With PDB File

```bash
python skills/protein-binder-design/scripts/target_analyzer.py <pdb_file> \
    --chain <target_chain_id> \
    --epitope <residue_numbers...> \
    --time <fast|standard> \
    --output analysis.json
```

### Without PDB (Verbal Description)

```bash
python skills/protein-binder-design/scripts/method_selector.py \
    --epitope-type <helix|loop|flat|concave|edge_strand|mixed> \
    --hydrophilicity <low|medium|high> \
    --time <fast|standard> \
    --affinity <standard|high>
```

## Output Interpretation

The analysis provides:

1. **Epitope Classification**: helix, loop, flat, concave, edge_strand, or mixed
   - This is the single most important factor for method selection

2. **Hydrophilicity Score** (0-1):
   - < 0.35: Hydrophobic (most design methods work well)
   - 0.35-0.55: Mixed (standard approaches)
   - > 0.55: Hydrophilic (consider beta-pairing RFdiffusion for edge-strands)

3. **Hotspot Recommendation**: 3-6 residues optimized for the epitope type
   - Helix: every 3-4 residues along one face
   - Loop: flanking residues, avoid flexible tip
   - Flat: distributed, minimum 8 Å apart
   - Edge-strand: alternating i, i+2 pattern
   - Concave: rim + pocket interior

4. **Method Rankings**: Top 3 methods with expected success rates

5. **Warnings**: Flexible regions, glycosylation sites, difficulty assessment

## Expert Knowledge

Consult `references/method_database.md` for detailed method profiles and `references/decision_trees.md` for the complete method selection decision tree.
