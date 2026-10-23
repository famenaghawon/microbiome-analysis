## Compositional transforms

Microbiome data is compositional: abundances are constrained to sum to
a fixed total per sample (due to sequencing depth). We apply **centered
log-ratio (CLR)** transform before any differential abundance test:

$$\text{clr}(x)_i = \log(x_i + 0.5) - \frac{1}{N} \sum_{j=1}^{N} \log(x_j + 0.5)$$

This makes relative differences interpretable and removes the
depth-vs-abundance confounding.

## Differential abundance

Two independent methods:

1. **Wilcoxon rank-sum on CLR abundances.** Non-parametric, robust to
   outliers. BH-FDR correction across genera at alpha = 0.05.
2. **ANCOM-BC-style bias-corrected CLR + Welch t.** Bias correction
   subtracts the median CLR per taxon; effect size is mean CLR
   difference; W-rank combines |effect| and -log10(q).

Both report per-taxon: effect size, raw p, BH-adjusted q, direction,
prevalence. Recomputed for each generation of the dataset.

## Alpha diversity

- **Richness** (observed OTUs): count of taxa with non-zero abundance
- **Shannon**: $- \sum p_i \log(p_i)$
- **Simpson**: $1 - \sum p_i^2$

Computed via scikit-bio `alpha_diversity`.

## Beta diversity

- **Bray-Curtis**: $\frac{\sum |a_i - b_i|}{\sum (a_i + b_i)}$
- **Jaccard**: 1 - |A n B| / |A u B| (presence/absence)

Both produce symmetric N x N matrices with zero diagonal.

## Ordination

**PCoA** (Principal Coordinates Analysis) via `skbio.stats.ordination.pcoa`.
Projects the distance matrix into 2D for visualization.
Reported variance explained: PC1 + PC2 percentages.

## PERMANOVA

Permutation-based multivariate ANOVA on the distance matrix.
Tests the null hypothesis "community composition does not differ by
group", with 999 permutations for the p-value.

## Machine learning

**Leakage-safe classification** via sklearn `Pipeline`:

1. CLR transform (fit-free)
2. `SelectKBest(f_classif, k=50)` — **fit on train fold only**
3. `RandomForestClassifier(500 trees)`

**Critical:** the feature selector is inside the Pipeline, so it sees
only train-fold data during cross-validation. This is enforced by a
**leakage sanity test**: on shuffled labels, balanced accuracy must be
within 0.35-0.65 (chance for binary classification).

5-fold stratified CV. Metrics: balanced accuracy (handles class
imbalance), AUC (for binary), confusion matrix per fold.

Feature importance: fit final model on all data, rank by
`feature_importances_` from the random forest. Reported top-20 genera.
