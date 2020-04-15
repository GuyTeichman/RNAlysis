=======
History
=======

1.3.5 (2020-??-??)
------------------
* This version introduces Pipelines - a quicker and more convenient way to apply a similar analysis pipeline to multiple Filter objects.

Added
******
* Added class Pipeline to filtering module, which applies a series of filter functions to specified Filter objects.


Changed
*******



Fixed
******
* Updated type-hinting of specific functions.
* Various bug fixes.


Removed
********






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
