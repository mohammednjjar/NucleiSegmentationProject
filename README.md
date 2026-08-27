# Supervised vs. Unsupervised Cell-Type Recovery on the Moffitt (2018) MERFISH Hypothalamic Preoptic Dataset

**A Feature-Contribution Benchmark**

M.Sc. Thesis — Computer Engineering, Bahçeşehir University
Author: Mohammed Abdulwahhab M Hamdi Al Najjar
Advisor: Assist. Prof. Dr. Lavdie Rada Ülgen

---

## Overview

Moffitt et al. (2018) published cell-type labels for ~875,000 MERFISH cells in the mouse hypothalamic preoptic region, derived by unsupervised Louvain clustering. They never trained a classifier on those labels, never reported an accuracy or F1 score for them, and never tested which input features carry the cell-type signal.

This work fills that gap. It places two supervised classifiers alongside three unsupervised clusterers on the same cells, the same labels and the same metric, while systematically varying the input feature block and the number of retained classes.

**Three questions:**

1. How well can the published labels be recovered by a supervised classifier, measured on an external holdout?
2. How does the best unsupervised method compare on an identical scale?
3. Which feature blocks — gene expression, spatial position, neighbourhood composition, boundary morphology — actually carry the signal?

---

## Key Findings

| Finding | Evidence |
|---|---|
| **Gene expression alone is sufficient** | No representation without genes exceeds 0.19 macro-F1 anywhere. Adding spatial position or neighbourhood to genes never improves recovery. |
| **Unsupervised nearly ties supervised at coarse granularity** | Gap of 0.0004 at K = 8 and 0.0028 at K = 9, widening to 0.1655 at K = 15 where rare classes enter. |
| **Over-segment-then-merge beats forcing exactly K** | Forcing the true class count costs 0.183–0.360 macro-F1, consistently across both graph algorithms and all granularities. |
| **Leiden ≥ Louvain ≫ GMM** | Leiden leads on macro-F1 throughout; both graph methods far ahead of the Gaussian Mixture Model. |
| **Failures are biological, not random** | The four weakest classes are oligodendrocyte maturation stages — a continuous process cut into discrete labels. |

---

## Dataset

| Constant | Value |
|---|---|
| Cells (after cleaning) | 874,768 |
| Genes | 155 |
| Published cell classes | 15 |
| Bregma sections | 12 (+0.26 to −0.29 mm) |
| Boundary-annotated subset | 35,522 |
| External holdout | up to 839,246 |

Source: Moffitt et al. (2018), released via Dryad (`doi:10.5061/dryad.8t8s248`) and GEO (`GSE113576`). The 2018 analysis code was never released, so the Louvain baseline is re-implemented rather than re-run.

---

## Experimental Design

**Splits.** The boundary-annotated subset is split 80/20, stratified by class, seed 42. The remaining cells form a fully external holdout, never touched during training or tuning.

**Granularities.** K = 8, 9, and 15. K = 15 is the full published set; K = 8 matches the eight clusters Dries et al. (2021) recovered from these data; K = 9 adds one further class.

**Feature sets (8).**

| Name | Contents |
|---|---|
| `genes_only` | 155 measured genes |
| `spatial_only` | distance to section centre, mean k-NN distance, local density |
| `neighbors_only` | mean expression of 10 nearest neighbours in-section |
| `genes_spatial` | genes + spatial |
| `genes_neighbors` | genes + neighbourhood |
| `all_combined` | genes + spatial + neighbourhood |
| `shape_only` | 12 polygon-morphology descriptors |
| `genes_shape` | genes + morphology |

**Models (5).**

| Model | Type | Role |
|---|---|---|
| Random Forest | Supervised | Tree ensemble, bootstrap + feature subsampling |
| LightGBM | Supervised | Sequential trees, leaf-wise growth |
| Leiden | Unsupervised | Modern graph community detection |
| Louvain | Unsupervised | Moffitt's original method, reproduced |
| Gaussian Mixture | Unsupervised | Non-graph, elliptical-density contrast |

**Evaluation.** Macro-F1 is the headline metric, so a class holding 0.3% of cells counts as much as one holding 37.1%. Unsupervised cluster indices are aligned to labels by Hungarian matching before scoring. ARI and NMI are reported alongside as assignment-free cross-checks.

---

## Pipeline

```
Step 0   Setup and library pinning
Step 1   Control panel — all constants in one place
Step 2   Load cell, boundary, and molecule tables
Step 3   Clean and build the external holdout
Step 4   Preprocess: standardise, PCA (unsupervised models only)
Step 5   Feature engineering: spatial + neighbour-gene blocks
Step 6   Baseline sweep — 5 models x feature sets x K
Step 7   Morphology ablation (boundary subset only)
Step 8-9 Optuna tuning, supervised and unsupervised
Step 10  Tuned vs. baseline comparison
Step 11-12 Parameter transfer across granularities
Step 13  Resolution sweep and clustering regimes
Step 14  Bregma section analysis
Step 15  Reproducibility manifest
```

---

## Reproducibility

- Python 3.12.13, NumPy 1.26.4, pandas 2.2.2, scikit-learn 1.5.1, SciPy 1.13.1, LightGBM 4.5.0, Scanpy 1.10.3, AnnData 0.10.8, python-igraph 0.11.6
- Single global seed of 42, applied to `PYTHONHASHSEED`, the Python `random` module, the NumPy generator, and every estimator, splitter, graph construction and Optuna sampler
- All constants written to `config_v20.json`; tuned hyperparameters to individual JSON files; full Optuna trial histories to CSV
- Pipeline identifier: `thesis_outputs_v20`
- Reproducible to approximately 1e−3 rather than bit-exact, due to multithreaded variation in LightGBM and the neighbour-graph routines

---

## Limitations

- All results rest on a single dataset; the feature-contribution finding may be specific to a tissue where cell classes are not strongly organised by position.
- Boundary polygons were released for one animal only, so the morphology ablation rests on a single biological replicate and is scored on the internal test rather than the external holdout.
- The spatial block is three summary statistics and captures no anatomical coordinates, spatial domains or cell-to-cell adjacency; the neighbourhood block is local expression averaging rather than true neighbourhood composition.
- Hungarian matching maximises overlap on the training rows, making unsupervised macro-F1 an optimistic estimate — which is why ARI and NMI are reported alongside.
- The benchmark operates at 15 broad classes; the ~70 finer neuronal subtypes are not tested.

---

## Citation

If you use this work, please cite the original dataset:

> Moffitt, J. R., Bambah-Mukku, D., Eichhorn, S. W., Vaughn, E., Shekhar, K., Perez, J. D., … Zhuang, X. (2018). Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. *Science*, 362(6416), eaau5324.

---

## Repository Contents

```
notebooks/     Colab pipeline (Steps 0–15)
outputs/       Result tables, figures, config and tuning artefacts
thesis/        LaTeX source and compiled thesis
