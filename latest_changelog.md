# 3.12.0 (2024-05-14)

I'm happy to announce the latest release of RNAlysis, which brings a variety of new features and improvements.

One of the key additions in this version is expanded support for advanced differential expression analysis using DESeq2 and Limma-Voom.
You can now test continuous covariates, perform likelihood ratio tests for factors, interactions, and polynomials, and take advantage of sample quality weights in Limma-Voom analysis.
We've also added several new visualization options, such as a p-value histogram plot and more customization choices for gene expression plots.

This release also includes various usability enhancements and bug fixes to provide a smoother analysis experience. We've clarified error messages, improved parameter annotations, and resolved issues that could lead to crashes or incorrect results in certain scenarios. RNAlysis now integrates with the latest versions of the OrthoInspector and ShortStack APIs as well.

As always, your support and feedback are greatly appreciated.
If you have any questions, encounter issues, or have suggestions for future development, please don't hesitate to reach out.
Happy analysis!
Guy

## Added

- Added support for advanced differential expression analysis with DESeq2/Limma-Voom, including testing continuous covariates, as well as likelihood ratio tests for factors, interactions, and polynomials.
- Added support for sample quality weights in Limma-Voom differential expression analysis.
- Added support for the 'cooksCutoff' parameter in DESq2 differential expression analysis.
- Added support for user-provided scaling factors in DESeq2 differential expression analysis.
- bowtie2 and kallisto paired-end modes now support the new file naming method 'smart', that attempts to automatically determine the common name of each pair of paired-end fastq files.
- Added a p-value histogram plot (DESeqFilter.pval_histogram) that displays the distribution of p-values in a differential expression table.
- Added new visualization options for 'plot expression of specific genes' function (CountFilter.plot_expression).
- Added new function 'Histogram of a column' (Filter.histogram) that generates a histogram of a column in a table.
- Added new function 'Concatenate two tables' (Filter.concatenate) that concatenates two tables along their rows.
- The filtering function 'Filter with a number filter' now supports filtering based on absolute values.

## Changed

- When running differential expression analysis, RNAlysis will automatically ensure that the order of samples in your design matrix matches the order of samples in your count matrix, avoiding erronious results.
- The 'Filter genes with low expression in all columns' function (CountFilter.filter_low_reads) now supports the 'n_samples' parameter, allowing users to filter genes with a minimal expression threshold in a specific number of samples.
- The 'Plot expressino of specific genes' function (CountFilter.plot_expression) now supports the 'jitter', 'group_names' and 'log_scale' parameters, allowing users to further customize the plot.
- The 'Scatter plot - sample VS sample' function (CountFilter.scatter_sample_vs_sample) now always displays highlighted points on top of the plot, making it easier to see which points are highlighted.
- Kallisto quantification (kallisto_quantify_single_end and kallisto_quantify_paired_end) now supports the 'summation_method' parameter, allowing users to choose between 'raw' and 'scaled_tpm' transcript summation methods. The default behavior of the functions did not change (it corresponds to 'scaled_tpm').
- Enrichment bar plots now have optional parameters that control font sizes for titles and labels.
- Moved the enrichment analysis and enrichment graphs actions to the "Enrichment" menu in the graphical interface to make the actions easier to find.
- Improved the clarity of error messages when attempting to read an invalid GTF file.
- RNAlysis now supports the latest version of OrthoInspector API.
- RNAlysis now supports the latest version of Ensembl Ortholog API.
- Improved annotation for the 'metric' parameter of the 'Hierarchical clustergram plot' function (CountFilter.clustergram).
- Improved performance of RNAlysis when generating automatic session reports.
- RNAlysis now offers default values for differential expression tables' column names.
- Functions that average replicates now display clearer group names by default.
- The RNAlysis interface to ShortStack now uses the most recent API (replaced 'knownRNAs' with 'known_miRNAs').
- When running differential expression, RNAlysis session reports will automatically include the auto-generated R script, as well as the sanitized design matrix used.
- Added optional parameters to all differential expression functions, allowing users to return a path to the auto-generated R script and data sanitized design matrix used.

## Fixed

- Fixed bug where some FASTQ/SAM functions could not be added to a FASTQ pipeline.
- Fixed bug where bowtie2 could not be run in/from directories with spaces in their names.
- Fixed bug where RNAlysis would crash when launched without an internet connection.
- Fixed bug that cause ID-mapping functions to raise an error when called from the MacOS stand-alone app (thanks to [Mitchzw](https://github.com/Mitchzw) in [#34](https://github.com/GuyTeichman/RNAlysis/issues/34)).
- Fixed bug that caused R package installations (DESeq2, limma, etc) to fail on some computers (thanks to [Celine-075](https://github.com/Celine-075) in [#35](https://github.com/GuyTeichman/RNAlysis/issues/35)).
- Fixed bug where Limma-Voom differential expression would fit numerical covariates as categorical factors.
- Fixed bug where generating enrichment bar plots with ylim='auto' would cause bars with 100% depletion (log2FC=-inf) to disappear.
- Fixed bug where defining 10 or more sample groups in the graphical interface would cause the groups to be ordered incorrectly in graphs.
- Fixed bug where the 'return_scaling_factors' argument would not return the normalization scaling factors on the graphical interface.
- Fixed various visual issues in the graphical interface
- Fixed bug where the 'filter_by_row_sum' function would raise an error when the table contains non-numerical columns
- Fixed bug where running enrichment on an empty gene set would raise an error.
- Fixed bug where RNAlysis would suggest resuming an auto-report from loaded session even when auto-report is turned off.
- Fixed bug where disabling auto-report in the middle of the session would raise errors when trying to create new graphs.
- Fixed bug where generating multiple gene expression plots (split_plots=True) with auto-generated report would only add the last graph to the session report.
- Fixed bug where the function 'normalize_to_quantile' generated unclear table names.
- Fixed bug where sometimes the first operation performed in a session would not display correctly in the automatic session report.

## New Contributors

- [Mitchzw] in [#34]
- [Celine-075] in [#35]
