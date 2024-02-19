=======
History
=======


3.12.0 (2024-??-??)
-------------------

Changed
********
* Enrichment bar plots now have optional parameters that control font sizes for titles and labels.

Fixed
*******
* Fixed bug where some FASTQ/SAM functions could not be added to a FASTQ pipeline.
* Fixed bug where bowtie2 could not be run in/from directories with spaces in their names.
* Fixed bug where RNAlysis would crash when launched without an internet connection.
* Fixed bug that cause ID-mapping functions to raise an error when called from the MacOS stand-alone app (thanks to `Mitchzw <https://github.com/Mitchzw>`_ in `#34 <https://github.com/GuyTeichman/RNAlysis/issues/34>`_).
* Fixed bug that caused R package installations (DESeq2, limma, etc) to fail on some computers (thanks to `Celine-075 <https://github.com/Celine-075>`_ in `#35 <https://github.com/GuyTeichman/RNAlysis/issues/35>`_).
* Fixed bug where generating enrichment bar plots with ylim='auto' would cause bars with 100% depletion (log2FC=-inf) to disappear.

New Contributors
*****************
* `Mitchzw`_ in `#34`_
* `Celine-075`_ in `#35`_

3.11.0 (2024-01-05)
-------------------
This release brings several exciting new features.
Notably, these inclue the ability to run functions from the Picardtools suite, and the ability to automatically save interactive session reports and later resume them from any saved session.
In addition, this release includes several visual upgrades, bug fixes, and quality-of-life improvements.
Happy analysis!

Added
*******
* *RNAlysis* can now run functions from the Picardtools suite, including conversion functions (BAM to SAM, SAM to FASTQ, FASTQ to SAM, etc), quality control (validate SAM), and post-processing functions (remove PCR duplicates, sort SAM, create BAM index).
* *RNAlysis* interactive session reports can now be resumed from any saved session, instead of having to start a new report from scratch. When loading a session created by *RNAlysis* 3.11 and beyond, you will have the option to resume the interactive report from the last saved state.

Changed
********
* CutAdapt adapter trimming functions can now receive an optional "new_filenames" parameter, which allows users to specify the names of the output files.
* Hierarchical clustergram plot (CountFilter.clustergram) now supports the 'colormap' parameter, which allows users to specify a custom color map for the plot.
* Hierarchical clustergram plot (CountFilter.clustergram) now displays continuous values on the color bar, instead of discrete values.
* Generally improved the looks of Hierarchical clustergram plot (CountFilter.clustergram).
* Previously-added functions (such as ortholog/paralog mapping, biotype summary by GTF file, etc) can now be applied to gene sets. Previously, some of these functions could only be applied to data tables.

Fixed
*******
* Function "Summarize feature biotypes (based on a reference table)" (biotypes_from_ref_table) now treats rows with missing values as "_missing_from_biotype_reference" instead of ignoring them entirely.
* Fixed bug where the Ensembl paralog-finding function would appear under the wrong tab in the graphical interface.
* Fixed bug where the description of the MA Plot function and parameters would not display correctly in the graphical interface.

3.10.1 (2023-11-22)
-------------------
Version 3.10.1 introduces several bug fixes, as well as well as support for random effect analysis in Limma-Voom differential expression.

Added
*******
* Limma-Voom differential expression can now fit mixed linear models containing a random effect (e.g. nested design).


Fixed
*******
* Fixed bug where trying to load sessions created with RNAlysis version 3.10.0 would result in an error.
* Fixed bug where using the OrthoInspector ortholog mapping function with database='auto' would sometimes fail to find an appropriate mapping database, even when one exists.
* Fixed bug where kallisto paired-end quantification window would display the 'read1' and 'read2' parameters twice.
* Fixed bug where empty sub-menus would appear under the FASTQ menu. These sub-menus will be implemented in future versions of RNAlysis.

3.10.0 (2023-10-31)
-------------------
I'm thrilled to introduce RNAlysis version 3.10.0.
This version includes features that were requested by users for a while, alongside quality-of-life improvements and bug fixes.
Here is a brief highlight of the most important additions:

**Ortholog Mapping:** *RNAlysis* can now map genes to their closest orthologs in different organisms.
You can map genes to their orthologs using four different databases - Ensembl, Panther, PhylomeDB, and OrthoInspector - extracting both one-to-one and one-to-many ortholog relationships and filtering them based on their reliability.

**Discovering Paralogs:** In the same vein, *RNAlysis* now facilitates the discovery of paralogs within a specific organism, using either the Ensembl or Panther databases.

**New visualization and analysis options for Principal Component Analysis (PCA):** We've introduced new functions and parameters to allow you to get more out of your principal component analysis.

I would also like to extend my personal apology for the delay in bringing you this update.
Due to personal reasons, this release, originally scheduled for the end of August, took longer than expected.
Your patience and support have been invaluable, and I'm eager to share these exciting additions with you.
Thank you for being a part of the RNAlysis community, and stay tuned for more updates in the near future!

Added
*******
* Added new functions to the filtering module that map genes to their closest orthologs in a different organism, using four different databases: Ensembl, Panther, PhylomeDB, and OrthoInspector.
* Added new functions to the filtering module that find paralogs of genes in a given organism, using two different databases: Ensembl and Panther.
* Added new function 'Sort table by contribution to a Principal Component (PCA)' (CountFilter.sort_by_principal_component), which allows sorting of genes in a count matrix by their contribution (gene loadings) to a principal component.
* Added a new parameter called 'legend' to 'Principal Component Analysis (PCA) plot' (CountFilter.pca), which allows users to display a legend on the PCA plot with a name for each sample group/color.

Changed
********
* RNAlysis now issues a warning when users run PCA or PCA-based functions on an unnormalized count matrix.
* The 'seek_fusion_genes' and 'learn_bias' arguments for kallisto quantification (fastq.kallisto_quant_single_end and fastq.kallisto_quant_paired_end), which were depracated in kallisto versions >0.48, are no longer displayed on the graphical interface. Old Pipelines that contain these arguments will still run, but new Pipelines will not contain them.
* Long-running functions now run in the background even when the 'inplace' parameter is set to True, instead of freezing the entire graphical interface.

Fixed
*******
* Fixed bug where functions would sometimes fail to run without displaying an error message.
* Fixed bug where progress bars in the graphical interface would sometimes not disappear after reaching 100% completion.
* RNAlysis should no longer display warning messages about graph layout when graphs are scaled down.
* Fixed bug where the clustergram function (CountFilter.clustergram) would raise an error with specific sets of dependecy versions.
* Loading tables no longer raises a depracation warning when using newer versions of Pandas.

3.9.2 (2023-06-23)
------------------
This patch contains bug fixes and improved functionality for enrichment lollipop plots,
as well as bug fixes for issues with the stand-alone version.

Changed
********
* Single-set enrichment result tables now contain observed/expected values based on the XL-mHG test cutoff.
* When loading a differential expression design matrix, RNAlysis now issues an error if the design matrix column names contain invalid characters.
* Updated the scaling of enrichment lollipop plots to make small 'observed' values easier to discern.

Fixed
*******
* Fixed bug where an error message would sometimes appear after RNAlysis finishes generating an automatic session report on the stand-alone app.
* Fixed bug where enrichment lollipop plots in horizontal mode would display the observed/expected values in reverse order.
* Fixed bug where enrichment lollipop plots and the 'show_exp' parameter would not work on single-set enrichment data.
* Fixed bug where sometimes the tab label of clustering/differential expression output tables would not match the name of the generated table.

3.9.1 (2023-06-19)
------------------

Version 3.9.1 of RNAlysis introduces several improvements and fixes to further improve your analysis experience.
The release includes new optional parameters for single-set enrichment functions, compatibility improvements with newer Python versions,
improved error messaging for R scripts, and adresses minor issues related to enrichment analysis, documentation, plotting parameters, and Pipeline saving.

Added
*******
* Added new optional parameters to single-set enrichment functions, allowing users to determine the top and bottom cutoffs for the XL-mHG test ("X" and "L").

Changed
********
* RNAlysis single-set enrichment analysis using the XL-mHG test now supports Python versions >= 3.8.
* RNAlysis stand-alone app is now built on Python 3.11, improving overall performance.
* Error messages caused by running R tools such as DESeq2, Limma-Voom, and FeatureCounts will now clearly state the reason the script failed, making it easier to understand what went wrong.

Fixed
*******
* Fixed bug where enrichment analysis would raise an error when running enrichment analysis on a gene set with no relevant annotations, or a gene set that does not intersect at all with the background gene set.
* Added missing documentation for plotting parameters in some enrichment functions.
* Depracation Warning should no longer appear when generating a box-plot or enhanced box-plot with scatter=True (CountFilter.box_plot, CountFilter.enhanced_box_plot)
* Fixed bug in featureCounts single-end mode where the 'output_folder' parameter could appear as disabled.
* Fixed bug where RNAlysis would raise an error message after saving a Pipeline, even when the Pipeline was saved successfully.

3.9.0 (2023-06-09)
------------------
Version 3.9.0 of *RNAlysis* introduces several enhancements and fixes to improve your experience.
The release includes additional enrichment plot styles, a new option for PCA plots,
the ability to load and save data tables in Parquet format, and new actions in the Help menu for reporting issues and suggesting improvements.
The update also improves the performance of various functions, ensures consistency in font and theme settings,
and addresses multiple bug fixes, including issues with automatic session reports and visualization functions.

Added
*******
* Added additional parameters to enrichment bar plots (enrichment.enrichment_bar_plot), including a new plot style ('lollipop') and observed/expected labels on the graph.
* Added a new parameter to Principal Component Analysis plots (CountFilter.pca) 'plot_grid', which can enable or disable adding a grid to PCA plots.
* RNAlysis can now load and save data tables in Parquet format (.parquet)
* Added new actions to the Help menu, allowing users to report issues, suggest issues, or open discussions.

Changed
********
* Functions in the FASTQ model are now added to automatic session reports.
* Many of the functions in RNAlysis should now run faster.
* Font type, size, and color for help tooltips should now match the global font settings.
* True/False toggle switches now scale with font size.
* When loading data tables into RNAlysis, you will now see only supported file formats by default.
* Clustering PCA plots are now plotted in proportion to the % variance explained by each PC.
* The legend in clustering PCA plots is now draggable.

Fixed
*******
* Fixed bug where data tables generated through the FASTQ model would not display properly in automatic session reports.
* Fixed bug where graphs generated through the Visualize Gene Sets window would not be added to automatic analysis report.
* When saving a file through the graphical interface, automatically-suggested filenames no longer contain illegal characters.
* Improved clarity of error message when R installation folder is not found.
* Fixed bug where some input parameter widgets in the RNAlysis graphical interface would not display properly.
* RNAlysis now provides a clearer warning message when attempting to run HDBSCAN clustering, if the hdbscan package is not installed.
* Label text in PCA plots and hierarchical clustergrams should no longer be cropped outside of the visible region of the plot.
* Fixed bug where some visualization functions, such as pair-plot (CountFilter.pairplot) would not display properly due to version mismatches between pandas and seaborn.
* Improved clarity of error messages when external apps' (kallisto, bowtie2, etc) installation folders are not found.
* Fixed bug where running the RNAlysis graphical interface on a new computer would sometimes raise an error (thanks to `NeuroRookie <https://github.com/NeuroRookie>`_ in `#25 <https://github.com/GuyTeichman/RNAlysis/issues/25>`_).
* Fixed a bug where the 'min_samples' parameter in HDBSCAN clustering could not be disabled.
* Fixed a bug where applying a function to a gene set with inplace=False would cause the new gene set to be called 'New Table'.
* Fixed a bug where RNAlysis would display the message "Pipeline saved successfully", even when the user cancels the save operation.

New Contributors
*****************
* `NeuroRookie`_ in `#25`_

3.8.0 (2023-05-07)
------------------
Version 3.8.0 of *RNAlysis* comes with several exciting new features, including the ability to generate interactive analysis reports automatically.
This feature allows users to create an interactive graph of all the datasets they loaded into RNAlysis, the functions applied, and the outputs generated during the analysis.
You can read more about this feature in the `Tutorial chapter <https://guyteichman.github.io/RNAlysis/build/tutorial.html#create-analysis-report>`_ and the `User Guide chapter <https://guyteichman.github.io/RNAlysis/build/user_guide_gui.html#rnalysis-interactive-analysis-reports>`_.

The new release also includes Pipelines for FASTQ functions, the ability to export normalization scaling factors, and other changes to improve the software's performance.
RNAlysis now supports Python 3.11, and many functions should now run faster. The software's graphic interface has also improved significantly, and users will now see a clearer error message when attempting to load unsupported file formats.
Lastly, the release also fixes several bugs.

Please note that, since Python 3.7 will be reaching end of life as of June 2023, new versions of *RNAlysis* will no longer support it.

Added
*******
* You can now generate interactive analysis reports automatically using the RNAlysis graphical interface. Read more about this feature `here <https://guyteichman.github.io/RNAlysis/build/user_guide_gui.html>`_.
* Added Pipelines for the FASTQ module (SingleEndPipeline and PairedEndPipeline), allowing you to apply a series of functions (adapter trimming, alignment, quantification, etc) to a batch of sequence/alignment files.
* Added new parameter 'return_scaling_factors' to normalization functions, that allows you to access the scaling factors calculated by RNAlysis.
* Added new parameter 'gzip_output' to CutAdapt adapter trimming (fastq.trim_adapters_single_end and fastq.trim_adapters_paired_end), allowing users to determine whether or not the output files will be gzipped.

Changed
*******
* RNAlysis now supports Python 3.11, and no longer supports Python 3.7.
* Many of the functions in RNAlysis should now run faster.
* The RNAlysis graphic interface should now boot up significantly faster.
* RNAlysis now shows an easier to understand error message when users attempt to load a table in an unsupported format (e.g. Excel files).
* CutAdapt adapter trimming (fastq.trim_adapters_single_end and fastq.trim_adapters_paired_end) now outputs non-gzipped files by default.
* Standardized all plotting functions in the filtering module to return a matplotlib Figure object, which can be further modified by users.

Fixed
*******
* RNAlysis failing to map gene IDs during GO enrichment analysis should no longer raise an error (thanks to `clockgene <https://github.com/clockgene>`_ in `#16 <https://github.com/GuyTeichman/RNAlysis/issues/16>`_).
* Fixed bug where the Command History window would not display history of the current tab immediately after clearing the current session.
* Fixed bug where adapter trimming would fail to run when using CutAdapt version >= 4.4.0.
* Fixed bug where 'Filter rows with duplicate names/IDs' (Filter.filter_duplicate_ids) would raise an error when applied to some tables.

New Contributors
*****************
* `clockgene`_ in `#16`_


3.7.0 (2023-04-07)
------------------
This version introduces small RNA read alignment using ShortStack, new filtering functions, a new optional parameter for Principal Component Analysis, improvements to the graphical interface, and bug fixes.

Added
*******
* Added small RNA read alignment using ShortStack (fastq.shortstack_align_smallrna).
* Added new filtering function 'Filter specific rows by name' (Filter.filter_by_row_name).
* Added new filtering function 'Filter rows with duplicate names/IDs' (Filter.filter_duplicate_ids).
* Function parameters in pop-up windows in the graphical interface can now be imported and exported.
* Added new parameter to Principal Component Analysis (CountFilter.pca) 'proportional_axes', that allows you to make the PCA projection axes proportional to the percentage of variance each PC explains.
* Improved clarity of error messages in the graphical interface.

Changed
*******
* Tables loaded into RNAlysis that use integer-indices will now be converted to use string-indices.
* Refactored CountFilter.from_folder into CountFilter.from_folder_htseqcount, and added a new CountFilter.from_folder method that accepts a folder of count files in any format.
* In the RNAlysis graphical interface, optional parameters that can be disabled will now display the hint "disable this parammeter?" instead of "disable input?".
* Added optional parameter 'ylim' to 'create enrichment bar-plot' function (enrichment.enrichment_bar_plot), allowing users to determine the Y-axis limits of the bar plot.
* Updated function signatures of 'Visualize gene ontology' and 'Visualize KEGG pathway' (enrichment.gene_ontology_graph and enrichment.kegg_pathway_graph) to make more sense.
* Removed parameter 'report_read_assignment_path' from featureCounts feature counting (fastq.featurecounts_single_end and fastq.featurecounts_paired_end).
* The RNAlysis graphical interface should now load more quickly.
* Progress bars in the graphical interface should now reflect elapsed/remaining time for tasks more accurately.

Fixed
*******
* Fixed bug in the function 'Split into Higly and Lowly expressed genes' (Filter.split_by_reads) where the two resulting tables would be named incorrectly (highly-expressed genes would be labeled 'belowXreads' and vice-versa).
* Fixed bug where the 'column' parameter of some functions ('Filter by percentile', 'Split by percentile', 'Filter with a number filter', 'Filter with a text filter') would not automatically detect column names in the graphical interface.
* Fixed bug where the 'numerator' and 'denominator' parameters of of the function 'Calculate fold change' would not automatically detect column names in the graphical interface.
* Fixed bug where tables with integer-indices could not be visualized properly through the graphical interface.
* Fixed bug in the function 'featureCounts single-end' (fastq.featurecounts_single_end) where setting parameter 'report_read_assignment' to any value other than None would raise an error.
* Functions that take column name/names as input (transform, filter_missing_values, filter_percentile, split_percentile) can now be applied to fold change tables (FoldCangeFilter objects).
* Fixed bug where the description for the 'n_bars' parameter of the 'create enrichment bar-plot' function (enrichment.enrichment_bar_plot) would not display properly.
* 'Visualize gene ontology' and 'Visualize KEGG pathway' (enrichment.gene_ontology_graph and enrichment.kegg_pathway_graph) now have proper parameter descriptions.
* Fixed bug where in-place intersection and difference in the filtering module would fail when using recent versions of pandas.
* Fixed bug where graphs generated through the Visualize Gene Sets window would not immediately display when using the RNAlysis stand-alone app.

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
* Changed argument order and default values in filtering.CountFilter.from_folder_htseqcount().
* Changed default title in scatter_sample_vs_sample().
* Changed default filename in CountFilter.fold_change().
* Settings are now saved in a .yaml format. Reading and writing of settings have been modified.
* Changed argument name 'deseq_highlight' to 'highlight' in scatter_sample_vs_sample(). It can now accept any Filter object.
* Updated documentation and default 'mode' value for FeatureSet.go_enrichment().
* Updated the signature and function of general.load_csv() to be clearer and more predictable.
* Changed argument names in CountFilter.from_folder_htseqcount().
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
