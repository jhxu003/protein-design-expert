"""
Fit scoring models on the Overath et al. 2025 meta-analysis dataset.

Uses the 3,760 experimentally tested binder designs across 15 targets to:
1. Fit logistic regression with Leave-One-Group-Out CV (LOGO-CV by target)
2. Validate our filter_engine thresholds against real experimental outcomes
3. Export per-target calibration parameters
4. Compare our composite scorer against data-fitted models

Data source: Overath et al. 2025, "Predicting Experimental Success in
De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally
Characterised Binders", bioRxiv 10.1101/2025.08.14.670059
GitHub: https://github.com/DigBioLab/de_novo_binder_scoring

Usage:
    python scripts/fit_scoring_model.py [--output-dir outputs/model_fitting]
"""

import csv
import json
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Tuple, Optional


# ── Data loading ────────────────────────────────────────────────────────────

DATASET_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "binder_dataset_slim.csv"
)


@dataclass
class BinderRecord:
    """One row from the meta-analysis dataset."""
    binder_id: str
    target_id: str
    is_binder: int  # 1 = confirmed binder, 0 = non-binder
    source: str
    binder_length: int
    # AF3 confidence
    af3_ipSAE_min: float
    af3_iptm_avg: float
    af3_ipae: float
    af3_min_pae_contact: float
    af3_pDockQ2_max: float
    af3_LIS: float
    # RMSD
    rmsd_binder: float  # RMSD_chA_aft_chB_align_input_af3
    # Rosetta
    af3_interface_sc: float
    af3_interface_dG: float
    af3_interface_dG_SASA_ratio: float
    af3_interface_dSASA: float
    input_interface_sc: float
    input_interface_dG: float
    # PyMOL
    af3_delta_SASA: float
    af3_num_interface_residues: float


def load_dataset(path: str = None) -> List[BinderRecord]:
    """Load the slim binder dataset."""
    path = path or DATASET_PATH
    records = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            def gf(key, default=float('nan')):
                v = row.get(key, "")
                try:
                    return float(v) if v else default
                except ValueError:
                    return default

            records.append(BinderRecord(
                binder_id=row.get("binder_id", ""),
                target_id=row.get("target_id", ""),
                is_binder=int(float(row.get("binder", "0"))),
                source=row.get("source", ""),
                binder_length=int(gf("A_length", 0)),
                af3_ipSAE_min=gf("af3_ipSAE_min"),
                af3_iptm_avg=gf("af3_iptm_avg"),
                af3_ipae=gf("af3_ipae"),
                af3_min_pae_contact=gf("af3_min_pae_contact"),
                af3_pDockQ2_max=gf("af3_pDockQ2_max"),
                af3_LIS=gf("af3_LIS"),
                rmsd_binder=gf("RMSD_chA_aft_chB_align_input_af3"),
                af3_interface_sc=gf("af3_rosetta_interface_sc"),
                af3_interface_dG=gf("af3_rosetta_interface_dG"),
                af3_interface_dG_SASA_ratio=gf("af3_rosetta_interface_dG_SASA_ratio"),
                af3_interface_dSASA=gf("af3_rosetta_interface_dSASA"),
                input_interface_sc=gf("input_rosetta_interface_sc"),
                input_interface_dG=gf("input_rosetta_interface_dG"),
                af3_delta_SASA=gf("af3_pymol_delta_SASA"),
                af3_num_interface_residues=gf("af3_pymol_num_interface_residues"),
            ))
    return records


# ── Pure-Python logistic regression ─────────────────────────────────────────
# No sklearn dependency — implements from scratch for zero-dependency policy.

def sigmoid(x: float) -> float:
    """Numerically stable sigmoid."""
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    else:
        ez = math.exp(x)
        return ez / (1.0 + ez)


def standardize(values: List[float]) -> Tuple[List[float], float, float]:
    """Standardize a list of values to zero mean, unit variance."""
    clean = [v for v in values if not math.isnan(v)]
    if not clean:
        return values, 0.0, 1.0
    mean = sum(clean) / len(clean)
    var = sum((v - mean) ** 2 for v in clean) / max(len(clean) - 1, 1)
    std = math.sqrt(var) if var > 0 else 1.0
    result = [(v - mean) / std if not math.isnan(v) else 0.0 for v in values]
    return result, mean, std


def impute_median(values: List[float]) -> List[float]:
    """Replace NaN with median."""
    clean = sorted(v for v in values if not math.isnan(v))
    if not clean:
        return [0.0] * len(values)
    n = len(clean)
    median = clean[n // 2] if n % 2 == 1 else (clean[n // 2 - 1] + clean[n // 2]) / 2
    return [v if not math.isnan(v) else median for v in values]


def fit_logistic_regression(
    X: List[List[float]],  # n_samples x n_features
    y: List[int],          # n_samples (0/1)
    C: float = 1.0,       # regularization (higher = less regularization)
    lr: float = 0.1,
    max_iter: int = 200,
) -> Tuple[List[float], float]:
    """
    Fit L2-regularized logistic regression via batch gradient descent.
    Uses vectorized inner loops and adaptive learning rate for speed.
    Returns (weights, bias).
    """
    n = len(y)
    if n == 0:
        return [0.0] * len(X[0]), 0.0
    d = len(X[0])
    w = [0.0] * d
    b = 0.0
    lam = 1.0 / C

    # Pre-transpose X for faster column access
    X_T = [[X[i][j] for i in range(n)] for j in range(d)]

    for iteration in range(max_iter):
        # Compute all predictions at once
        preds = [0.0] * n
        for i in range(n):
            z = b
            for j in range(d):
                z += w[j] * X[i][j]
            preds[i] = sigmoid(z)

        # Compute errors
        errors = [preds[i] - y[i] for i in range(n)]

        # Gradient for each weight using pre-transposed X
        for j in range(d):
            gj = sum(errors[i] * X_T[j][i] for i in range(n)) / n + lam * w[j]
            w[j] -= lr * gj

        # Gradient for bias
        gb = sum(errors) / n
        b -= lr * gb

    return w, b


def predict_proba(X: List[List[float]], w: List[float], b: float) -> List[float]:
    """Predict probabilities."""
    return [sigmoid(sum(w[j] * x[j] for j in range(len(w))) + b) for x in X]


def average_precision(y_true: List[int], y_scores: List[float]) -> float:
    """Compute average precision (area under precision-recall curve)."""
    # Sort by score descending
    pairs = sorted(zip(y_scores, y_true), key=lambda x: -x[0])
    tp = 0
    fp = 0
    ap_sum = 0.0
    n_pos = sum(y_true)
    if n_pos == 0:
        return 0.0

    for score, label in pairs:
        if label == 1:
            tp += 1
            precision = tp / (tp + fp)
            ap_sum += precision
        else:
            fp += 1

    return ap_sum / n_pos


def precision_at_k(y_true: List[int], y_scores: List[float], k: int) -> float:
    """Precision in top-k ranked by score."""
    pairs = sorted(zip(y_scores, y_true), key=lambda x: -x[0])
    top_k = pairs[:k]
    if not top_k:
        return 0.0
    return sum(label for _, label in top_k) / len(top_k)


def f1_optimal_threshold(y_true: List[int], y_scores: List[float]) -> Tuple[float, float]:
    """Find threshold that maximizes F1 score. Returns (threshold, f1)."""
    pairs = sorted(zip(y_scores, y_true), key=lambda x: -x[0])
    n_pos = sum(y_true)
    if n_pos == 0:
        return 0.5, 0.0

    best_f1 = 0.0
    best_thr = 0.5
    tp = 0
    fp = 0
    for i, (score, label) in enumerate(pairs):
        if label == 1:
            tp += 1
        else:
            fp += 1
        fn = n_pos - tp
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
        if f1 > best_f1:
            best_f1 = f1
            best_thr = score

    return best_thr, best_f1


# ── Model fitting with LOGO-CV ─────────────────────────────────────────────

@dataclass
class FeatureSet:
    """A named set of features for model comparison."""
    name: str
    feature_names: List[str]
    extractor: object  # callable: BinderRecord -> List[float]


def extract_ipsae_only(r: BinderRecord) -> List[float]:
    return [r.af3_ipSAE_min]

def extract_ipsae_rmsd(r: BinderRecord) -> List[float]:
    return [r.af3_ipSAE_min, r.rmsd_binder]

def extract_ipsae_sc(r: BinderRecord) -> List[float]:
    return [r.af3_ipSAE_min, r.input_interface_sc]

def extract_3feat(r: BinderRecord) -> List[float]:
    return [r.af3_ipSAE_min, r.rmsd_binder, r.input_interface_sc]

def extract_full(r: BinderRecord) -> List[float]:
    return [
        r.af3_ipSAE_min, r.af3_iptm_avg, r.rmsd_binder,
        r.input_interface_sc, r.af3_interface_dG_SASA_ratio,
        r.af3_pDockQ2_max, r.af3_min_pae_contact,
    ]


FEATURE_SETS = [
    FeatureSet("ipSAE_min_only", ["af3_ipSAE_min"], extract_ipsae_only),
    FeatureSet("ipSAE+RMSD", ["af3_ipSAE_min", "RMSD_binder"], extract_ipsae_rmsd),
    FeatureSet("ipSAE+SC", ["af3_ipSAE_min", "input_interface_sc"], extract_ipsae_sc),
    FeatureSet("ipSAE+RMSD+SC", ["af3_ipSAE_min", "RMSD_binder", "input_interface_sc"], extract_3feat),
    FeatureSet("7_feature_full", [
        "af3_ipSAE_min", "af3_iptm_avg", "RMSD_binder",
        "input_interface_sc", "af3_dG_SASA_ratio",
        "af3_pDockQ2_max", "af3_min_pae_contact"
    ], extract_full),
]


def logo_cv(
    records: List[BinderRecord],
    feature_set: FeatureSet,
    C: float = 1.0,
) -> Dict:
    """
    Leave-One-Group-Out cross-validation by target_id.
    Returns per-target and aggregate metrics.
    """
    # Group by target
    targets = defaultdict(list)
    for r in records:
        targets[r.target_id].append(r)

    target_ids = sorted(targets.keys())
    per_target = {}
    all_y_true = []
    all_y_scores = []

    for held_out in target_ids:
        # Train on all other targets
        train_records = []
        for tid in target_ids:
            if tid != held_out:
                train_records.extend(targets[tid])
        test_records = targets[held_out]

        # Extract features
        X_train_raw = [feature_set.extractor(r) for r in train_records]
        X_test_raw = [feature_set.extractor(r) for r in test_records]
        y_train = [r.is_binder for r in train_records]
        y_test = [r.is_binder for r in test_records]

        n_features = len(X_train_raw[0])

        # Impute + standardize per feature
        stats = []  # (mean, std) per feature
        X_train = [[] for _ in range(len(train_records))]
        X_test = [[] for _ in range(len(test_records))]

        for j in range(n_features):
            col_train = [X_train_raw[i][j] for i in range(len(train_records))]
            col_test = [X_test_raw[i][j] for i in range(len(test_records))]

            # Impute
            col_train = impute_median(col_train)
            # Use train median for test imputation
            clean_train = sorted(v for v in col_train if not math.isnan(v))
            train_median = clean_train[len(clean_train) // 2] if clean_train else 0.0
            col_test = [v if not math.isnan(v) else train_median for v in col_test]

            # Standardize using train stats
            col_train_std, mean, std = standardize(col_train)
            col_test_std = [(v - mean) / std if std > 0 else 0.0 for v in col_test]

            stats.append((mean, std))
            for i, v in enumerate(col_train_std):
                X_train[i].append(v)
            for i, v in enumerate(col_test_std):
                X_test[i].append(v)

        # Fit model
        w, b = fit_logistic_regression(X_train, y_train, C=C, lr=0.1, max_iter=200)

        # Predict
        scores = predict_proba(X_test, w, b)

        # Evaluate
        n_pos = sum(y_test)
        ap = average_precision(y_test, scores) if n_pos > 0 else float('nan')
        p5 = precision_at_k(y_test, scores, 5) if len(scores) >= 5 else float('nan')
        p10 = precision_at_k(y_test, scores, 10) if len(scores) >= 10 else float('nan')
        p20 = precision_at_k(y_test, scores, 20) if len(scores) >= 20 else float('nan')
        thr, f1 = f1_optimal_threshold(y_train, predict_proba(X_train, w, b))

        per_target[held_out] = {
            "n_designs": len(test_records),
            "n_binders": n_pos,
            "success_rate": round(n_pos / len(test_records), 4) if test_records else 0,
            "average_precision": round(ap, 4) if not math.isnan(ap) else None,
            "precision_at_5": round(p5, 4) if not math.isnan(p5) else None,
            "precision_at_10": round(p10, 4) if not math.isnan(p10) else None,
            "precision_at_20": round(p20, 4) if not math.isnan(p20) else None,
            "f1_threshold": round(thr, 4),
            "weights": [round(wi, 4) for wi in w],
            "bias": round(b, 4),
            "standardization": [(round(m, 4), round(s, 4)) for m, s in stats],
        }

        all_y_true.extend(y_test)
        all_y_scores.extend(scores)

    # Aggregate
    aps = [v["average_precision"] for v in per_target.values()
           if v["average_precision"] is not None]
    aps.sort()
    median_ap = aps[len(aps) // 2] if aps else 0.0
    mean_ap = sum(aps) / len(aps) if aps else 0.0

    global_ap = average_precision(all_y_true, all_y_scores)
    global_p10 = precision_at_k(all_y_true, all_y_scores, 10)
    global_p20 = precision_at_k(all_y_true, all_y_scores, 20)

    return {
        "feature_set": feature_set.name,
        "feature_names": feature_set.feature_names,
        "n_features": len(feature_set.feature_names),
        "n_targets": len(target_ids),
        "n_total": len(records),
        "n_binders": sum(r.is_binder for r in records),
        "median_AP": round(median_ap, 4),
        "mean_AP": round(mean_ap, 4),
        "global_AP": round(global_ap, 4),
        "global_precision_at_10": round(global_p10, 4),
        "global_precision_at_20": round(global_p20, 4),
        "per_target": per_target,
    }


# ── Filter validation ──────────────────────────────────────────────────────

def validate_filters(records: List[BinderRecord]) -> Dict:
    """
    Validate our filter_engine thresholds against real experimental outcomes.
    Tests how well our filtering rules separate binders from non-binders.
    """
    # Define filter scenarios
    filter_scenarios = {
        "ipSAE_top100": lambda r: True,  # rank-based, handled separately
        "RMSD<2.0": lambda r: r.rmsd_binder < 2.0 if not math.isnan(r.rmsd_binder) else False,
        "RMSD<3.0": lambda r: r.rmsd_binder < 3.0 if not math.isnan(r.rmsd_binder) else False,
        "RMSD<3.73_data_validated": lambda r: r.rmsd_binder < 3.73 if not math.isnan(r.rmsd_binder) else False,
        "SC>0.62_data_validated": lambda r: not math.isnan(r.input_interface_sc) and r.input_interface_sc > 0.62,
        "ipTM>0.60": lambda r: r.af3_iptm_avg > 0.60 if not math.isnan(r.af3_iptm_avg) else False,
        "ipTM>0.80": lambda r: r.af3_iptm_avg > 0.80 if not math.isnan(r.af3_iptm_avg) else False,
        "combined_RMSD+SC": lambda r: (
            (not math.isnan(r.rmsd_binder) and r.rmsd_binder < 3.73) and
            (not math.isnan(r.input_interface_sc) and r.input_interface_sc > 0.62)
        ),
    }

    results = {}
    for name, filter_fn in filter_scenarios.items():
        if name == "ipSAE_top100":
            # Rank-based: per target, take top 100 by ipSAE_min
            targets = defaultdict(list)
            for r in records:
                targets[r.target_id].append(r)

            total_selected = 0
            total_binders_selected = 0
            total_binders = 0
            for tid, recs in targets.items():
                valid = [(r, r.af3_ipSAE_min) for r in recs if not math.isnan(r.af3_ipSAE_min)]
                valid.sort(key=lambda x: x[1])  # lower ipSAE = better
                top = valid[:min(100, len(valid))]
                total_selected += len(top)
                total_binders_selected += sum(1 for r, _ in top if r.is_binder)
                total_binders += sum(1 for r in recs if r.is_binder)

            results[name] = {
                "n_selected": total_selected,
                "n_binders_selected": total_binders_selected,
                "n_total_binders": total_binders,
                "precision": round(total_binders_selected / total_selected, 4) if total_selected > 0 else 0,
                "recall": round(total_binders_selected / total_binders, 4) if total_binders > 0 else 0,
            }
        else:
            passed = [r for r in records if filter_fn(r)]
            binders_passed = sum(1 for r in passed if r.is_binder)
            total_binders = sum(1 for r in records if r.is_binder)

            results[name] = {
                "n_selected": len(passed),
                "n_binders_selected": binders_passed,
                "n_total_binders": total_binders,
                "precision": round(binders_passed / len(passed), 4) if passed else 0,
                "recall": round(binders_passed / total_binders, 4) if total_binders > 0 else 0,
            }

    return results


# ── Per-target calibration ──────────────────────────────────────────────────

def compute_per_target_calibration(records: List[BinderRecord]) -> Dict:
    """
    Compute per-target statistics for threshold calibration.
    For each target: binder/non-binder distribution of key metrics.
    """
    targets = defaultdict(list)
    for r in records:
        targets[r.target_id].append(r)

    calibration = {}
    for tid in sorted(targets.keys()):
        recs = targets[tid]
        binders = [r for r in recs if r.is_binder]
        non_binders = [r for r in recs if not r.is_binder]

        def metric_stats(values):
            clean = sorted(v for v in values if not math.isnan(v))
            if not clean:
                return {"median": None, "q25": None, "q75": None, "n": 0}
            n = len(clean)
            return {
                "median": round(clean[n // 2], 4),
                "q25": round(clean[n // 4], 4),
                "q75": round(clean[3 * n // 4], 4),
                "n": n,
            }

        calibration[tid] = {
            "n_designs": len(recs),
            "n_binders": len(binders),
            "success_rate": round(len(binders) / len(recs), 4),
            "binder_metrics": {
                "ipSAE_min": metric_stats([r.af3_ipSAE_min for r in binders]),
                "iptm_avg": metric_stats([r.af3_iptm_avg for r in binders]),
                "rmsd_binder": metric_stats([r.rmsd_binder for r in binders]),
                "interface_sc": metric_stats([r.input_interface_sc for r in binders]),
                "interface_dG": metric_stats([r.af3_interface_dG for r in binders]),
                "delta_SASA": metric_stats([r.af3_delta_SASA for r in binders]),
            },
            "non_binder_metrics": {
                "ipSAE_min": metric_stats([r.af3_ipSAE_min for r in non_binders]),
                "iptm_avg": metric_stats([r.af3_iptm_avg for r in non_binders]),
                "rmsd_binder": metric_stats([r.rmsd_binder for r in non_binders]),
                "interface_sc": metric_stats([r.input_interface_sc for r in non_binders]),
                "interface_dG": metric_stats([r.af3_interface_dG for r in non_binders]),
                "delta_SASA": metric_stats([r.af3_delta_SASA for r in non_binders]),
            },
        }

    return calibration


# ── Export fitted weights for design_scorer ─────────────────────────────────

def export_fitted_weights(cv_result: Dict, output_path: str):
    """
    Export the globally-fitted logistic regression weights as JSON,
    ready to be loaded by design_scorer.py.
    """
    # Average weights across all per-target folds
    all_weights = []
    all_biases = []
    all_stats = []
    for target_data in cv_result["per_target"].values():
        all_weights.append(target_data["weights"])
        all_biases.append(target_data["bias"])
        all_stats.append(target_data["standardization"])

    n_folds = len(all_weights)
    n_feats = len(all_weights[0])

    avg_weights = [sum(all_weights[f][j] for f in range(n_folds)) / n_folds
                   for j in range(n_feats)]
    avg_bias = sum(all_biases) / n_folds
    avg_stats = [
        (sum(all_stats[f][j][0] for f in range(n_folds)) / n_folds,
         sum(all_stats[f][j][1] for f in range(n_folds)) / n_folds)
        for j in range(n_feats)
    ]

    output = {
        "model": "logistic_regression",
        "feature_set": cv_result["feature_set"],
        "feature_names": cv_result["feature_names"],
        "weights": [round(w, 6) for w in avg_weights],
        "bias": round(avg_bias, 6),
        "standardization_params": [
            {"mean": round(m, 6), "std": round(s, 6)}
            for m, s in avg_stats
        ],
        "performance": {
            "median_AP": cv_result["median_AP"],
            "mean_AP": cv_result["mean_AP"],
            "global_AP": cv_result["global_AP"],
        },
        "training_data": {
            "source": "Overath et al. 2025",
            "n_designs": cv_result["n_total"],
            "n_binders": cv_result["n_binders"],
            "n_targets": cv_result["n_targets"],
        },
    }

    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Fit and validate scoring models on the Overath et al. 2025 binder dataset"
    )
    parser.add_argument(
        "--dataset", default=DATASET_PATH,
        help="Path to binder_dataset_slim.csv"
    )
    parser.add_argument(
        "--output-dir", default="outputs/model_fitting",
        help="Directory for output files"
    )
    parser.add_argument(
        "--C", type=float, default=1.0,
        help="Regularization parameter for logistic regression"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    print("=" * 70)
    print("BINDER SCORING MODEL FITTING")
    print("Data: Overath et al. 2025 Meta-Analysis")
    print("=" * 70)

    records = load_dataset(args.dataset)
    print(f"\nLoaded {len(records)} designs")
    print(f"  Binders: {sum(r.is_binder for r in records)}")
    print(f"  Non-binders: {sum(not r.is_binder for r in records)}")
    print(f"  Targets: {len(set(r.target_id for r in records))}")

    # 1. Fit models with LOGO-CV
    print("\n" + "=" * 70)
    print("LOGISTIC REGRESSION — LOGO-CV RESULTS")
    print("=" * 70)

    all_results = []
    for fs in FEATURE_SETS:
        print(f"\n--- {fs.name} ({len(fs.feature_names)} features) ---")
        result = logo_cv(records, fs, C=args.C)
        all_results.append(result)

        print(f"  Median AP: {result['median_AP']:.4f}")
        print(f"  Mean AP:   {result['mean_AP']:.4f}")
        print(f"  Global AP: {result['global_AP']:.4f}")
        print(f"  P@10: {result['global_precision_at_10']:.4f}")
        print(f"  P@20: {result['global_precision_at_20']:.4f}")

        print(f"\n  Per-target AP:")
        for tid, td in sorted(result["per_target"].items()):
            ap_str = f"{td['average_precision']:.3f}" if td['average_precision'] else "N/A"
            print(f"    {tid:<20s}: AP={ap_str}  "
                  f"({td['n_binders']}/{td['n_designs']} binders)")

    # Save model comparison
    comparison_path = os.path.join(args.output_dir, "model_comparison.json")
    with open(comparison_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nModel comparison saved to {comparison_path}")

    # 2. Export best model weights
    best_3feat = [r for r in all_results if r["feature_set"] == "ipSAE+RMSD+SC"][0]
    weights_path = os.path.join(args.output_dir, "fitted_weights.json")
    export_fitted_weights(best_3feat, weights_path)
    print(f"Fitted weights saved to {weights_path}")

    best_full = [r for r in all_results if r["feature_set"] == "7_feature_full"][0]
    weights_full_path = os.path.join(args.output_dir, "fitted_weights_full.json")
    export_fitted_weights(best_full, weights_full_path)
    print(f"Full model weights saved to {weights_full_path}")

    # 3. Validate filter thresholds
    print("\n" + "=" * 70)
    print("FILTER THRESHOLD VALIDATION")
    print("=" * 70)

    filter_results = validate_filters(records)
    for name, stats in filter_results.items():
        print(f"\n  {name}:")
        print(f"    Selected: {stats['n_selected']}/{len(records)}")
        print(f"    Binders in selection: {stats['n_binders_selected']}/{stats['n_total_binders']}")
        print(f"    Precision: {stats['precision']:.4f}")
        print(f"    Recall: {stats['recall']:.4f}")

    filter_path = os.path.join(args.output_dir, "filter_validation.json")
    with open(filter_path, "w") as f:
        json.dump(filter_results, f, indent=2)
    print(f"\nFilter validation saved to {filter_path}")

    # 4. Per-target calibration
    print("\n" + "=" * 70)
    print("PER-TARGET CALIBRATION")
    print("=" * 70)

    calibration = compute_per_target_calibration(records)
    for tid, cal in sorted(calibration.items()):
        bm = cal["binder_metrics"]["ipSAE_min"]
        nm = cal["non_binder_metrics"]["ipSAE_min"]
        bm_iptm = cal["binder_metrics"]["iptm_avg"]
        print(f"\n  {tid} ({cal['n_binders']}/{cal['n_designs']} = {cal['success_rate']:.1%}):")
        if bm["median"] is not None and nm["median"] is not None:
            print(f"    ipSAE_min:  binder median={bm['median']:.3f}  non-binder median={nm['median']:.3f}")
        if bm_iptm["median"] is not None:
            print(f"    iptm_avg:   binder median={bm_iptm['median']:.3f}")

    cal_path = os.path.join(args.output_dir, "per_target_calibration.json")
    with open(cal_path, "w") as f:
        json.dump(calibration, f, indent=2)
    print(f"\nPer-target calibration saved to {cal_path}")

    # 5. Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nBest 3-feature model (ipSAE+RMSD+SC):")
    print(f"  Median AP = {best_3feat['median_AP']:.4f}")
    print(f"  Global AP = {best_3feat['global_AP']:.4f}")
    print(f"\nBest 7-feature model:")
    print(f"  Median AP = {best_full['median_AP']:.4f}")
    print(f"  Global AP = {best_full['global_AP']:.4f}")
    print(f"\nOutput files:")
    print(f"  {comparison_path}")
    print(f"  {weights_path}")
    print(f"  {weights_full_path}")
    print(f"  {filter_path}")
    print(f"  {cal_path}")
    print()


if __name__ == "__main__":
    main()
