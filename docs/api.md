# API Reference

## Data loaders

- `microbiome_analysis.data.gmrepo.load_gmrepo(path, min_depth, min_prevalence)`
- `microbiome_analysis.data.loaders.load(source, path, max_samples, seed)`
- `microbiome_analysis.data.synthetic.make_synthetic(seed, effect, ...)`

## Preprocessing

- `microbiome_analysis.preprocess.transforms.clr(counts, pseudocount)`
- `microbiome_analysis.preprocess.transforms.relative(counts)`
- `microbiome_analysis.preprocess.transforms.rarefy(counts, depth, seed)`
- `microbiome_analysis.preprocess.depth.sample_depths(counts)`

## Differential abundance

- `microbiome_analysis.differential.wilcoxon.wilcoxon_da(clr_df, labels, group_a, group_b, alpha)`
- `microbiome_analysis.differential.ancom_bc.ancom_bc_da(counts, labels, ref_group, alpha)`

## Diversity

- `microbiome_analysis.diversity.alpha.alpha_table(counts)`
- `microbiome_analysis.diversity.beta.beta_distances(counts)`
- `microbiome_analysis.diversity.ordination.ordination_coords(distances, n_components)`
- `microbiome_analysis.diversity.permanova.run_permanova(distances, groups, permutations)`

## Machine learning

- `microbiome_analysis.ml.classifier.build_pipeline(k, seed)`
- `microbiome_analysis.ml.classifier.evaluate(counts, labels, k, folds, seed)`
- `microbiome_analysis.ml.classifier.train_final_and_rank_features(counts, labels, k, seed)`
- `microbiome_analysis.ml.classifier.leakage_sanity_check(seed, k)`

## Plotting

- `microbiome_analysis.plotting.alpha_plot(alpha_df, out_path)`
- `microbiome_analysis.plotting.pcoa_plot(coords, explained, labels, out_path, title)`
