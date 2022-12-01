=======
History
=======
3.3.0 (2022-12-01)
------------------
* This version introduced quality-of-life improvements to the graphical user interface.

Added
******
* Added a Frequently Asked Questions page, and linked all *RNAlysis* help material inside the graphical interface Help menu.
* Pipelines can now be edited and deleted through the Pipeline menu of the graphical interface.

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
