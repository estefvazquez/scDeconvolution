# Dataset Human
Data for PDAC was downloaded from GEO (GSE155698). We have used only the PDAC samples.
File markers_for_deconv.csv was generated by performing differential expression analysis on Scanpy. Running the function sc.tl.rank_genes_groups(adata, groupby='New_label', method = 'wilcoxon', corr_method='benjamini-hochberg').
File metadata_complete.csv contains the updated metadata with inferCNV scores (we've used infercnvpy for this. T cells used as reference).
scDeconv_FeatureSelection_ResultsComparison.R contains code for feature selection based on the differentially expressed genes previously identified. Moreover, we show how to compare results between Bayes and Bisque, and correlate the estimated cell fractions to Tumor Purity as a validation exercise.
