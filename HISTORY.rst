=======
History
=======

1.4.0 (2020-??-??)
------------------
* This version introduces many new ways to perform enrichment analysis and to visualize your results. Moreover, this version introduces Pipelines - a quicker and more convenient way to apply a similar analysis pipeline to multiple Filter objects.
* This version includes changes to names and signatures of some functions in the module, as elaborated below.

Added
******
* Added class Pipeline to filtering module, which applies a series of filter functions to specified Filter objects.
* Added class RankedSet to enrichment module, which accepts a ranked list of genes/features, and can perform single-list enrichment analysis
* Added RankedSet.enrich_single_list(), which can perfofm single-list enrichment analysis of user-defined attributes using XL-mHG test (see `Eden et al. (PLoS Comput Biol, 2007) <https://dx.doi.org/10.1371/journal.pcbi.0030039>`_  and `Wagner (PLoS One, 2015) <https://dx.doi.org/10.1371/journal.pone.0143196>`_ ).
* Added more visualization functions, such CountFilter.enhanced_box_plot().
* Added FeatureSet.change_set_name(), to give a new 'set_name' to a FeatureSet object.


Changed
*******
* CountFilter.pca() can now be plotted without labeled points.
* CountFilter.violin_plot() now accepts a y_title argument.
* Added more optional arguments to visualization functions such as CountFilter.violin_plot() and CountFilter.clustergram().
* Automatic filenames for Filter objects should now reflect more clearly the operations that were performed.
* The DataFrame returned by enrich_randomization() and enrich_randomization_parallel() now contains the additional column 'data_scale', determined by the new optional argument 'data_scale'.
* The columns 'n obs' and 'n exp' in the DataFrame returned by enrich_randomization() and enrich_randomization_parallel() were renamed to 'obs' and 'exp' respectively.
* FeatureSets no longer support in-place set operations (intersection, union, difference, symmetric difference). Instead, these functions return a new FeatureSet.


Fixed
******
* Updated type-hinting of specific functions.
* Filter.biotypes() and FeatureSet.biotypes() now support Biotype Reference Tables with different column names.
* Generally improved performance of RNAlysis.
* Various bug fixes.





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
* Filter.biotypes() and FeatureSet.biotypes() now report genes that don't appear in the biotype reference table.
* Filter.biotypes() can now give a long-form report with descriptive statistics of all columns, grouped by biotype.
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
