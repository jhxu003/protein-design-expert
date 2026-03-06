# Meta-Analysis Dataset Reference

## Source

- **Paper**: "Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders"
- **Authors**: Overath, Rygaard, Jacobsen, Brasas, Morell, Sormanni, Jenkins (DTU Denmark)
- **bioRxiv DOI**: 10.1101/2025.08.14.670059
- **GitHub**: https://github.com/DigBioLab/de_novo_binder_scoring
- **Dataset**: https://doi.org/10.5281/zenodo.15722219

## Dataset Summary

- 3,760 experimentally tested binder designs
- 436 confirmed binders (11.6% overall success rate)
- 15 structurally diverse targets
- 200+ structural, energetic, and confidence features per design
- Prediction tools evaluated: AF2, AF3, ColabFold, Boltz-1

## Per-Target Statistics

| Target | n_designs | n_binders | Success Rate | Source |
|--------|-----------|-----------|-------------|--------|
| EGFR | 434 | 28 | 6.5% | Cao et al, Adaptyv Bio competition |
| FGFR2 | 2123 | 193 | 9.1% | Cao et al |
| IL10Ra | 22 | 2 | 9.1% | Bennet et al prospective |
| IL2Ra | 66 | 6 | 9.1% | Bennet et al prospective |
| IL7Ra | 171 | 38 | 22.2% | Watson et al (RFdiffusion) |
| InsulinR | 117 | 20 | 17.1% | Watson et al |
| LTK | 33 | 3 | 9.1% | Bennet et al prospective |
| Mdm2 | 96 | 55 | 57.3% | Watson et al |
| Pdl1 | 95 | 12 | 12.6% | Watson et al |
| SARS_CoV2_RBD | 99 | 9 | 9.1% | Cao et al |
| TrkA | 128 | 9 | 7.0% | Cao et al, Watson et al |
| VirB8 | 99 | 9 | 9.1% | Cao et al |
| pMHC_NY1 | 132 | 43 | 32.6% | pMHC screening |
| pMHC_SILSY1 | 96 | 2 | 2.1% | pMHC screening |
| sntx | 49 | 7 | 14.3% | Vazquez-Torres et al |

## Feature Ranking (AF3, by Median Average Precision)

### Individual Features (Top 20)

| Rank | Feature | Median AP |
|------|---------|-----------|
| 1 | af3_ipSAE_min | 0.5399 |
| 2 | af3_iptm_avg | 0.5212 |
| 3 | af3_ipSAE_d0chn | 0.5170 |
| 4 | af3_ipSAE_d0dom | 0.5102 |
| 5 | af3_iptm_model_0 | 0.5098 |
| 6 | af3_min_pae_contact | 0.4898 |
| 7 | af3_ipSAE_max | 0.4803 |
| 8 | af3_ipSAE_avg | 0.4756 |
| 9 | af3_LIS | 0.4587 |
| 10 | af3_ptm_model_0 | 0.4557 |
| 11 | af3_pDockQ2_max | 0.4524 |
| 12 | af3_ipae | 0.4358 |
| 13 | af3_ptm_avg | 0.4285 |
| 14 | af3_pDockQ2_min | 0.4242 |
| 15 | af3_pDockQ_min | 0.3897 |
| 16 | af3_pDockQ_max | 0.3897 |
| 17 | af3_rosetta_interface_dG | 0.3092 |
| 18 | RMSD_chA_aft_chB_align_input_af3 | 0.3038 |
| 19 | af3_rosetta_interface_dG_SASA_ratio | 0.2849 |
| 20 | af3_pymol_percent_nonpolar | 0.2178 |

**af3_ipSAE_min** is the single best predictor across all tools and features.

### Interaction Features (Top 10)

| Rank | Feature Interaction | Median AP |
|------|-------------------|-----------|
| 1 | af3_LIS * input_rosetta_interface_sc | 0.5752 |
| 2 | af3_rosetta_interface_dG_SASA_ratio * af3_ipSAE_avg | 0.5674 |
| 3 | af3_ipSAE_d0chn * input_rosetta_interface_sc | 0.5638 |
| 4 | af3_rosetta_interface_dG_SASA_ratio * af3_ipSAE_min | 0.5574 |
| 5 | af3_rosetta_interface_sc * af3_ipSAE_min | 0.5554 |
| 6 | af3_ipSAE_min * input_rosetta_interface_sc | 0.5530 |
| 7 | af3_ipSAE_d0dom * input_rosetta_interface_sc | 0.5530 |
| 8 | af3_ipSAE_max * input_rosetta_interface_sc | 0.5302 |
| 9 | af3_ipSAE_min * input_rosetta_interface_dG_SASA_ratio | 0.5298 |
| 10 | af3_ipSAE_min * input_sap_delta | 0.5215 |

**af3_LIS * input_rosetta_interface_sc** is the single best interaction feature, exceeding all individual features.

## Best Logistic Regression Models (Nested LOCO-CV)

### Baseline: ipSAE_min only

- Median AP: 0.540 (IQR 0.347)
- Median Precision@F1: 0.465 (IQR 0.365)
- Median F1: 0.462

### + RMSD binder

- Median AP: 0.546 (IQR 0.325)
- Median Precision@F1: 0.500

### + input_interface_sc (shape complementarity)

- Median AP: 0.545 (IQR 0.342)
- Median Precision@F1: 0.500

### Both (ipSAE_min + RMSD + interface_sc)

- Median AP: 0.573 (IQR 0.303)
- Median Precision@F1: 0.538

### Optimal Filter Thresholds (from logistic regression)

- **RMSD binder** < 3.73 Angstrom
- **input_interface_sc** > 0.62

## Key Findings

1. **ipSAE_min is the best single predictor** -- 1.4x better than ipAE by average precision.
2. **Simple 2-3 feature logistic regression generalizes as well as XGBoost** with hundreds of features.
3. **AF3 features generally outperform AF2, ColabFold, and Boltz-1** across all evaluation metrics.
4. **Feature interactions** (e.g., LIS x interface_sc) can exceed single features in predictive power.
5. **Per-target success rates vary enormously** (2.1% to 57.3%), suggesting target-specific calibration is important.
6. **Rosetta dG and shape complementarity** provide complementary signal to AF3 confidence metrics.

## Comparison: AF3 vs AF2 vs ColabFold vs Boltz-1

The paper found AF3 features consistently outperform other predictors. The logistic regression greedy feature selection with AF3 features achieves the highest median AP. When combining features from all predictors ("All"), only marginal improvement is observed over AF3-only models.

## Practical Implications for Binder Design Scoring

1. **Rank designs by AF3 ipSAE_min** as the primary sort criterion.
2. **Apply RMSD binder < 3.73 Angstrom filter** to remove misfolded designs.
3. **Apply shape complementarity > 0.62 filter** to ensure good interface packing.
4. **For top-k selection**: filtering by RMSD + SC before ranking by ipSAE_min improves precision at all k values.
5. **Complex multi-feature scoring models do NOT significantly outperform 2-3 feature models** -- simplicity generalizes better across targets.
