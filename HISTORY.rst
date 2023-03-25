=======
History
=======

3.6.2 (2023-03-25)
------------------
This version introduces quality-of-life changes to the graphical interface, as well as bug fixes.

Added
*******
* Sample groupings in functions such as PCA, Pair-plot, etc., can now be imported and exported through the graphical interface.

Fixed
*******
* Fixed bug where the stand-alone Mac version of RNAlysis would sometimes fail to map gene IDs (directly or in enrichment analysis).

3.6.1 (2023-03-22)
------------------
This version introduces minor bug fixes.

Changed
********
* DESeq2 automatic installation should now work more reliably.

Fixed
******
* Fixed bug where PCA plots would not display on the stand-alone app unless another visualization function was applied afterwards.
* Fixed bug where Pipelines that contain specific functions (such as translating gene IDs/filtering biotypes from GTF file) would fail to run through the graphical interface.
* GO Annotations annotated by ComplexPortal are now supported by RNAlysis.

3.6.0 (2023-03-07)
------------------
This version introduces improvements to the usability and clarity of the graphic interface,
new methods for automatic estimation of the number of clusters in a dataset,
and various bug fixes.

Added
******
* Added three new methods for automatic estimation of the number of clusters in a dataset: Callinski-Harabasz, Davies-Bouldin, and Bayesian Information Criterieon.
* Added a 'Close all Figures' actions to the 'View' menu of the *RNAlysis* graphic interface.
* Added an 'interactive' parameter to Volcano Plots (DESeqFilter.volcano_plot) and 'Scatter Sample Vs Sample' (CountFilter.scatter_sample_vs_sample), allowing user to label data points interactively by clicking on them.
* Added more optional plotting parameters to Volcano Plots (DESeqFilter.volcano_plot) and 'Scatter Sample Vs Sample' (CountFilter.scatter_sample_vs_sample).

Changed
********
* Progress bars are now integrated into the main *RNAlysis* window instead of opening as a dialog box.
* Information about running proccesses and functions is now displayed in the main *RNAlysis* window.
* It is now possible to cancel queued jobs through the *RNAlysis* graphic interface.
* When loading multiple data tables at the same time, it is now possible to change the table type of all data tables at once, instead of one-by-one.

Fixed
******
* RNAlysis KEGG enrichment should now match the new KEGG annotation format from March 1st 2023.
* Fixed bug where importing *RNAlysis* would raise ImportError when cutadapt is not installed.
* Fixed bug where the 'Run' button in the Enrichment Analysis window would grey out whenever the enrichment dataset is changed.
* Fixed bug where the *RNAlysis* stand-alone versions were unable to export Figures in specific formats (e.g. PDF, SVG).
* Fixed bug where functions that depend on R scripts (such as DESeq2 and limma) would sometimes fail to run on MacOS (thanks to Matthias Wilm and `sandyl27 <https://github.com/sandyl27>`_ in `#12 <https://github.com/GuyTeichman/RNAlysis/issues/12>`_).
* Fixed bug where running limma-voom with a design matrix whose column names contained spaces or special characterse would raise an error.
* Fixed bug where the 'highlight' parameter of CountFilter.scatter_sample_vs_sample would not work when used through the graphic interface.
* Fixed bug where enrichment analysis would sometimes fail to run when 'exclude_unannotated_genes' is set to False.
* Fixed bug where translate_gene_ids() would fail for RankedSet objects.
* Fixed bug where filtering gene sets by user-defined attributes (FeatureSet.filter_by_attribute()) would occasionally fail to run.

New Contributors
*****************
* `sandyl27`_ in `#12`_

3.5.2 (2023-02-23)
------------------
This version includes bug fixes for a few quality-of-life issues which were introduced in version RNAlysis 3.5.0.

Changed
********
* It is now possible to change the parallel backend of performance-intensive functions such as clustering an enrichment analysis in non-standalone versions of RNAlysis.
* Expanded the Frequently Asked Questions page.
* Added Perl as a dependency for Windows users on the How to Install page.
* Automatic row colours in column-picking tables should no longer mismatch with font colors in a way that obscures visibility.

Fixed
*****
* Fixed bug where occasionally newly-created tabs would open twice instead of once.
* Fixed bug where the 'Add Function' button of the Pipeline window would appear in the wrong location.
* Fixed bug where RNAlysis sessions saved after version 3.5.0 which contain gene sets would raise an error when loaded.
* Fixed typos in the RNAlysis tutorial.

3.5.1 (2023-02-22)
------------------
This version introduces minor bug fixes.

Fixed
*****
* Fixed bug where the *RNAlysis* stand-alone app would sometimes fail to run CutAdapt (thanks to Matthias Wilm).
* Fixed bug where the User Guide action in the graphical interface would point to the Tutorial, and vice versa.
* The X and Y axis labels on volcano plots should now format the 'log' in the label correctly.

3.5.0 (2023-02-08)
------------------
This version introduces new features such as differential expression analysis with the Limma-Voom pipeline,
customizable databases for the quick-search function, basic filtering and procrssing functions for gene sets,
improved progammatic API for FeatureSet and RankedSet objects, and RPKM and TPM read count normalization.
Several changes have been made to improve the user experience, including updated documentation,
improved clarity of function tooltips, clearer output formats, and faster download speeds for tutorial videos.

Added
*******
* Added differential expression analysis with the Limma-Voom pipelines (CountFilter.differential_expression_limma_voom)
* You can now select which databases to display in the right-click quick-search menu, using the settings menu.
* Gene sets now support some basic operations (filtering, gene ID translating, etc) through the graphical interface.
* enrichment.FeatureSet and enrichment.RankedSet now support some filtering operations from the filtering module (such as filtering by user-defined attributes, GO terms, or KEGG pathways).
* Added reads-per-kilobase-million (RPKM) and transcripts-per-million (TPM) normalization methods (CountFilter.normalize_to_rpkm() and CountFilter.normalize_to_tpm()).

Changed
********
* Classes enrichment.FeatureSet and enrichment.RankedSet now inherit from Python base-class set, and can be interacted with like other Python sets. The old API and attributes of these classes were maintained as they were.
* Improved documentation for some functions.
* Function selection tooltips should now display information more clearly.
* Pipelines that contain consecutive clustering/splitting functions will now return their outputs in a clearer format.
* Enrichment bar-plots should now adjust the x-axis limits more tightly to fit the displayed data.
* Improved clarity of automatic titles in enrichment plots.
* Download/update speed of tutorial videos has improved significantly.

Fixed
******
* Fixed bug where Pipelines would not always properly run 'translate_gene_ids'

3.4.2 (2023-02-01)
------------------
This version introduces minor bug fixes.

Fixed
******
* Fixed bug where updating RNAlysis through the graphical interface would not update some of the optional dependencies.
* Fixed typos in the documentation.

3.4.0 (2023-02-01)
------------------
From this release forward, *RNAlysis* is made available as a stand-alone app for Windows and MacOS. You can download these stand-alone versions from the GitHub Releases page.
In addition, new features were added, including new plots, filtering functions, integration of the external tools bowtie2 and featureCounts, and the ability to generate Gene Ontology Graphs and KEGG Pathway Graphs without running enrichment analysis from scratch.

Added
******
* Added a Scree Plot (explained variance per PC) to Principle Component Analysis
* Added CountFilter.split_by_principal_component(), allowing users to filter genes based on their contributions (loadings) to PCA Principal Components.
* Added Filter.drop_columns
* Added support for the Sharpened Cosine distance metric in clustering analyses
* KEGG enrichment can now generate KEGG pathway graphs for pathways that were found to be statistically significant
* Added functions to the enrichment module that can generate KEGG Pathway or Gene Ontology plots based on previously-generated enrichment results
* You can now clear the *RNAlysis* cache from the *RNAlysis* GUI, or through the general module.
* Added bowtie2 alignment to the fastq module.
* Added FeatureCounts feature-counting to the fastq module.
* You can now choose whether or not to discard genes from enrichment analysis if they have no annotations associated with them.
* When right-clicking on a specific cell in a table or line in a gene set view, a context menu will open, allowing you to copy the associated value, or look it up in one of few common biology databases.
* Added sections to the programmatic user guide about the `fastq` module.

Changed
********
* Replaced the 'parallel' parameter in enrichment functions with the 'parallel_backend' parameter, allowing users to choose which parallel backend (if any) will be used in the function.
* Added 'parallel_backend' parameter to all clustering functions under the filtering module.
* When generating Gene Ontology/KEGG Pathway graphs, users can choose whether or not to generate the figure in an additional separate file.
* Updated type annotations of some functions to be more precise and helpful (for example, setting a lower bound on some int function parameters).
* The colorbar ticks in enrichment bar plots now match the axis ticks on the main axis.
* Slight improvements in GUI performance, stability, and looks.
* Slight improvements in performance of enrichment analysis when examining a small number of attributes.
* enrichment.plot_enrichment() was replaced by enrichment.enrichment_bar_plot().
* CountFilter.differential_expression() has new optional parameter `output_folder`, which allows users to save the generated data tables and the R script that generated them into a specified folder.

Fixed
******
* In CountFilter.differential_expression_deseq2(), fixed a bug where design matrix files with non-comma delimiters would cause an error (thanks to `Mintxoklet <https://github.com/Mintxoklet>`_ in `#7 <https://github.com/GuyTeichman/RNAlysis/issues/7>`_)
* Fixed bug where setup.py would install a directory named tests into site-packages folder (thanks to `Bipin Kumar <https://github.com/kbipinkumar>`_ in `#9 <https://github.com/GuyTeichman/RNAlysis/issues/9>`_)
* Fixed bug where the windows of some functions (differential expression, adapter trimming, etc) did not show a link to the function's documentation page.
* Fixed typos in some parts of the *RNAlysis* documentation
* When filtering a table by a single user-defined attribute, the automatic table name will now be more informative about the operation applied.
* Fixed bug where occasionally a Pipeline or Function would generate multiple tables of the same name, but only one of them will appear in the GUI.
* Fixed bug where occasionally significance asterisks on enrichment bar-plots would appear in the wrong location.
* Fixed bug where fastq.create_kallisto_index() (Create Kallisto Index) would not make use of the `make_unique` parameter (thanks to Matthias Wilm)

Removed
********
* Removed the previously-deprecated functions `enrichment.enrich_randomization()` and `enrichment.enrich_hypergeometric()`.



New Contributors
*****************
* `Mintxoklet`_ in `#7`_
* `Bipin Kumar`_ in `#9`_
* Matthias Wilm

3.3.0 (2022-12-02)
------------------
* This version introduced quality-of-life improvements to the graphical user interface.

Added
******
* Added a Frequently Asked Questions page, and linked all *RNAlysis* help material inside the graphical interface Help menu.
* Pipelines can now be edited and deleted through the Pipeline menu of the graphical interface.
* Added Contributing page to the documentation, with basic guidelines on how you can contribute to the *RNAlysis* project!

Changed
*******
* All open tabs are now always visible in the main menu screen. Tab names are now shortened with ellipsis if nessecary.
* The right-click context menu of the main menu tabs now allows users to open a new tab at a specific position, or close a specific tab/tabs to the right/tabs to the left/all other tabs.
* *RNAlysis* documentation is now split into GUI documentation (quick-start video guide, tutorial, GUI user guide), and programmatic documentation (programmatic user guide)
* Improved readability of *RNAlysis* logs
* Pipelines are now exported with additional metadata - the version of *RNAlysis* they were exported from, and the date and time it was exported. This metadata should not affect Pipelines that were created in older versions, and does not affect the way Pipelines are applied to data tables.

Fixed
******
* *RNAlysis* now warns users if they attempt to overwrite an existing Pipeline.
* Fixed an incorrect keyboard shortcut for Export Pipeline action

3.2.2 (2022-11-25)
------------------


Fixed
******
* Fixed bug with DESeq2 automatic installation on Windows computers.

3.2.1 (2022-11-25)
------------------

Changed
*******
* Updated citation information for *RNAlysis*

Fixed
******
* Fixed typos in the *RNAlysis* tutorial

3.2.0 (2022-11-23)
------------------
* This version introduces quality-of-life changes to the graphical user interface, functions for translating gene IDs and running differential expression analysis, and extends RNAlysis to support Python versions 3.9 and 3.10.

Added
******
* Added Filter.translate_gene_ids()
* Added CountFilter.differential_expression_deseq2()
* Added Filter.filter_by_kegg_annotations()
* Added Filter.filter_by_go_annotations()
* Added CountFilter.average_replicate_samples()
* Added fastq module that contains adapter-trimming functions utilizing CutAdapt, and mRNA-sequencing quantification using kallisto.

Changed
*******
* Added additional plotting parameters to visualization functions.
* Improved performance of some aspects of the graphical user interface.
* RNAlysis' basic features are now supported on Python versions 3.9 and 3.10.
* CountFilter.pca() now generates a plot for *every* pair of Principal Components requested by the user.
* CountFilter.split_clicom() now supports clustering each batch of replicates separately, using the 'replicates_grouping' parameter
* Biotype-based filtering and summary can now be done based on GTF annotation files instead of a Biotype Reference Table.
* Filter.biotypes() was refactored into Filter.biotypes_from_ref_table()
* Filter.filter_biotype() was refactored into Filter.filter_biotype_from_ref_table()

Fixed
******
* Users can now queue multiple computationally-intense enrichment/clustering tasks while another task is running.
* Fixed a bug where sometimes some function parameters would disappear from the graphical user interface.
* Fixed a bug where exceptions during computationally-intense tasks would cause *RNAlysis* to crash.
* Auxillary windows are now properly minimized when analysis starts, and restored when analysis ends or encounters an error.

3.1.0 (2022-10-16)
------------------
* This version introduces new count matrix normalization methods, as well as MA plots and minor bug fixes.

Added
******
* Added the visualization function ma_plot() for CountFilter
* Added functions for the normalization functions Relative Log Ratio (RLE), Trimmed Mean of M-values (TMM), Median of Ratios (MRN), Quantile normalization (quantile)

Changed
*******
* CountFilter.normalize_to_rpm() was renamed to CountFilter.normalize_to_rpm_htseqcount(), and was supplemented by the more general function for normalizing to Reads Per Million CountFilter.normalize_to_rpm()

Fixed
******
* Fixed a bug where some elements of the graphical user interface would not display correctly

3.0.1 (2022-10-12)
------------------
* This version fixes a bug with displaying the tutorial videos in the graphical user interface.


3.0.0 (2022-10-10)
------------------
* This version introduces a graphical user interface for RNAlysis, as well as new functions for KEGG Pathways enrichment analysis.


Added
******
* RNAlysis now includes a graphical user interface
* Pipelines can now be imported and exported
* Enrichment and single-set-enrichment for KEGG pathway data

Changed
*******
* Added function FeatureSet.user_defined_enrichment(), which will replace FeatureSet.enrich_hypergeometric() and FeatureSet.enrich_randomization()
* Updated signature of enrichment.venn_diagram
* enrichment.venn_diagram and enrichment.upset_plot can now be generated on a user-supplied FIgure
* Clustering functions now apply a power transform to count data prior to clustering by default
* Non-deprecated enrichment functions no longer filter the background set by biotype by default
* Changed signature of CountFilter.pca, CountFilter.box_plot, CountFilter.enhanced_box_plot, CountFilter.clustergram, and CountFilter.pairplot to ensure consistency among visualization functions.

Fixed
******
* enrichment.venn_diagram can now be plotted with outlines when the circles are unweighted
* Fixed bug in Pipeline.apply_to() where a Filter object would be returned even when the Pipeline was applied inplace


2.1.1 (2022-07-05)
------------------
* This version fixes issues with running GO enrichment that resulted from recent changes to UniProt's API.  Moreover, this version slightly improves the performance of some functions.

Changed
*******
* Fixed issues with running GO enrichment that resulted from changes to UniProt's API.
* Some functions that fetch annotations now cache their results, leading to improved runtimes.
* Updated the documentation of some functions to better reflect their usage and input parameters.

2.1.0 (2022-04-16)
------------------
* This version introduces multiple new features, as well as generally improved graphs and quality-of-life changes.

Added
******
* GO enrichment can now generate Ontology Graphs for the statistically significant GO terms.
* Added CountFilter.split_clicom(), an implementation of the CLICOM ensemble-based clustering method (Mimaroglu and Yagci 2012).
* Added Filter.transform(), a method that can transform your data tables with either predefined or user-defined transformations.

Changed
*******
* CountFilter.pairplot() now uses a logarithmic scale by default.
* Visually improved the graphs generated by many functions, including CountFilter.pairplot() and CountFilter.plot_expression().
* The clusters resulting from all clustering functions are now sorted by size instead of being sorted randomly.

Fixed
******
* Minor bug fixes.


2.0.1 (2022-04-02)
------------------
* This version introduces small bug fixes, as well as a new function in the Filtering module.

Added
******
* Added Filter.majority_vote_intersection(), which returns a set/string of the features that appear in at least (majority_threhold * 100)% of the given Filter objects/sets.

Changed
*******
* When mapping/inferring taxon IDs during GO enrichment analysis, organisms will now be prioritized based on their taxon ID values (numerically lower IDs will be considered to be more relevant).

Fixed
******
* Fixed bug that occured when mapping/inferring taxon IDs during GO enrichment analysis, where integer taxon IDs would be matched by name similarity before trying an exact ID match, leading to spurious matches.
* Fixed bug that occursed when plotting clustering results with style='all' on Python 3.8.

2.0.0 (2021-12-05)
------------------
* This version introduces new method to cluster your read count matrices, including K-Means/Medoids clustering, Hierarchical clustering, and HDBSCAN.
* This version introduces many new ways to perform enrichment analysis and to visualize your results, including highly customizable GO Enrichment, enrichment based on ranked lists of genes, and enrichment for non-categorical attributes.
* This version introduces Pipelines - a quicker and more convenient way to apply a particular analysis pipeline to multiple Filter objects.
* This version improves the performance of many functions in RNAlysis, and in particular the performance of randomization tests.
* This version includes changes to names and signatures of some functions in the module, as elaborated below.


Added
******
* Added class Pipeline to filtering module, which applies a series of filter functions to specified Filter objects.
* Added CountFilter.split_kmeans(), CountFilter.split_kmedoids(), CountFilter.split_hierarchical() and CountFilter.split_hdbscan(), which split your read count matrices into clusters with similar expression patterns.
* Added class RankedSet to enrichment module, which accepts a ranked list of genes/features, and can perform single-list enrichment analysis
* Added RankedSet.single_set_enrichment(), which can perfofm single-list enrichment analysis of user-defined attributes using XL-mHG test (see `Eden et al. (PLoS Comput Biol, 2007) <https://dx.doi.org/10.1371/journal.pcbi.0030039>`_  and `Wagner (PLoS One, 2015) <https://dx.doi.org/10.1371/journal.pone.0143196>`_ ).
* Added FeatureSet.go_enrichment() and RankedSet.single_set_go_enrichment(), which let you compute Gene Ontology enrichment for any organism of your choice, and filter the GO annotations used according to your preferences.
* Added FeatureSet.enrich_hypergeometric(), which can perform enrichment analysis using the Hypergeometric Test.
* Added more visualization functions, such CountFilter.enhanced_box_plot().
* Added FeatureSet.change_set_name(), to give a new 'set_name' to a FeatureSet object.


Changed
*******
* FeatureSet.enrich_randomization_parallel() was deprecated. Instead, you can compute your enrichment analysis with parallel computing by calling FeatureSet.enrich_randomization() with the argument 'parallel_processing=True'. Moreover, parallel session will now start automatically if one was not already active.
* Improved running time of enrich_randomization() about six-fold.
* Filter objects can be created from any delimiter-separated file format (.csv, .tsv, .txt, etc).
* CountFilter.pca() can now be plotted without labeled points.
* Filter.index_string is now sorted by the current order of indices in the Filter object, instead of by alphabetical order.
* CountFilter.violin_plot() now accepts a y_title argument.
* Added more optional arguments to visualization functions such as CountFilter.violin_plot() and CountFilter.clustergram().
* Automatic filenames for Filter objects should now reflect more clearly the operations that were performed.
* The DataFrame returned by enrich_randomization() and enrich_randomization_parallel() now contains the additional column 'data_scale', determined by the new optional argument 'data_scale'.
* The columns 'n obs' and 'n exp' in the DataFrame returned by enrich_randomization() and enrich_randomization_parallel() were renamed to 'obs' and 'exp' respectively.
* FeatureSets no longer support in-place set operations (intersection, union, difference, symmetric difference). Instead, these functions return a new FeatureSet.
* Filter.biotypes_from_ref_table() now accepts the boolean parameter 'long_format' instead of the str parameter 'format'.
* Filter.biotypes_from_ref_table() and FeatureSet.biotypes_from_ref_table() now count features which do not appear in the Biotype Reference Table as '_missing_from_biotype_reference' instead of 'not_in_biotype_reference'.

Fixed
******
* Updated type-hinting of specific functions.
* Filter.biotypes_from_ref_table() and FeatureSet.biotypes_from_ref_table() now support Biotype Reference Tables with different column names.
* Generally improved performance of RNAlysis.
* Fixed bug in Filter.filter_percentile() where the value at the exact percentile speficied (e.g. the median for percentile=0.5) would be removed from the Filter object.
* Fixed bug in enrichment.FeatureSet, where creating a FeatureSet from input string would result in an empty set.
* Various minor bug fixes.





1.3.5 (2020-05-27)
------------------
* This version introduces minor bug fixes and a few more visualization options.

Added
******
* Added Filter.filter_missing_values(), which can remove rows with NaN values in some (or all) columns.
* Added the visualization function CountFilter.box_plot().

Changed
*******
* Updated docstrings and printouts of several functions.
* Slightly improved speed and performance across the board.
* Filter.feature_string() is now sorted alphabetically.
* Enrichment randomization functions in the enrichment module now accept a 'random_seed' argument, to be able to generate consistent results over multiple sessions.
* Enrichment randomization functions can now return the matplotlib Figure object, in addition to the results table.


Fixed
******
* Fixed DepracationWarning on parsing functions from the general module.
* Fixed bug where saving csv files on Linux systems would save the files under the wrong directory.
* Fixed a bug where UTF-8-encoded Reference Tables won't be loaded correctly
* Fixed a bug where enrichment.upsetplot() and enrichment.venn_diagram() would sometimes modify the user dict input 'objs'.
* Fixed a bug in CountFilter.pairplot where log2 would be calculated without a pseudocount added, leading to division by 0.




1.3.4 (2020-04-07)
------------------
* This version fixed a bug that prevented installation of the package.


Changed
*******
* Updated docstrings and printouts of several functions


Fixed
******
* Fixed a bug with installation of the previous version






1.3.3 (2020-03-28)
------------------
* First stable release on PyPI.


Added
******
* Added Filter.sort(), and upgraded the functionality of Filter.filter_top_n().
* Added UpSet plots and Venn diagrams to enrichment module.
* User-defined biotype reference tables can now be used.
* Filter operations now print out the result of the operation.
* Enrichment randomization tests now also support non-WBGene indexing.
* Filter.biotypes_from_ref_table() and FeatureSet.biotypes_from_ref_table() now report genes that don't appear in the biotype reference table.
* Filter.biotypes_from_ref_table() can now give a long-form report with descriptive statistics of all columns, grouped by biotype.
* Added code examples to the user guide and to the docstrings of most functions.


Changed
*******
* Changed argument order and default values in filtering.CountFilter.from_folder().
* Changed default title in scatter_sample_vs_sample().
* Changed default filename in CountFilter.fold_change().
* Settings are now saved in a .yaml format. Reading and writing of settings have been modified.
* Changed argument name 'deseq_highlight' to 'highlight' in scatter_sample_vs_sample(). It can now accept any Filter object.
* Updated documentation and default 'mode' value for FeatureSet.go_enrichment().
* Updated the signature and function of general.load_csv() to be clearer and more predictable.
* Changed argument names in CountFilter.from_folder().
* Modified names and signatures of .csv test files functions to make them more comprehensible.
* Renamed 'Filter.filter_by_ref_table_attr()' to 'Filter.filter_by_attribute()'.
* Renamed 'Filter.split_by_ref_table_attr()' to 'Filter.split_by_attribute()'.
* Renamed 'Filter.norm_reads_with_size_factor()' to 'Filter.normalize_with_scaling_factors()'. It can now use any set of scaling factors to normalize libraries.
* Renamed 'Filter.norm_reads_to_rpm()' to 'Filter.normalize_to_rpm()'.
* Made some functions in the general module hidden.


Fixed
******
* Various bug fixes


Removed
********
* Removed the 'feature_name_to_wbgene' module from RNAlysis.






1.3.2 (2019-12-11)
------------------

* First beta release on PyPI.
