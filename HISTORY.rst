=======
History
=======
1.3.3 (2020-02-25)
------------------
* First stable release on PyPI.

Added
++++++
* Added Filter.sort(), and upgraded the functionality of Filter.filter_top_n()
* Added UpSet plots and Venn diagrams to enrichment module
* User-defined biotype reference tables can now be used.
* Filter operations now print out the result of the operation
* Enrichment randomization tests now also support non-WBGene indexing.


Changed
+++++++
* Changed argument order and default values in filtering.CountFilter.from_folder()
* Changed default title in scatter_sample_vs_sample()
* Changed default filename in CountFilter.fold_change()
* Settings are now saved in a .yaml format. Reading and writing of settings have been modified.
* Changed argument name 'deseq_highlight' to 'highlight' in scatter_sample_vs_sample(). It can now accept any Filter object.
* Updated documentation and default 'mode' value for FeatureSet.go_enrichment()

Fixed
++++++
* Various bug fixes

Removed
++++++++
* Removed the 'feature_name_to_wbgene' module from RNAlysis.






1.3.2 (2019-12-11)
------------------

* First beta release on PyPI.
