############################
User guide
############################



***********************************
*RNAlysis* graphical user interface
***********************************

Open the program
=================
If you installed *RNAlysis* from pypi, you can open the graphical user interface by either executing the command `rnalysis-gui` from your terminal, or by typing the following code into a Python console::

    >>> from *RNAlysis* import gui
    >>> gui.run_gui()

Load a table
============
Choose a file from your computer, and click 'start' to load it into *RNAlysis*.

.. image:: ../../rnalysis/gui/videos/load_table.webp

Examine your table
==================
You will now be able to see an overview of your data, including the table's name, type, and dimensions. Click the 'View full table' button to see your table in its entirety.

.. image:: ../../rnalysis/gui/videos/view_table.webp

Filter your table
=================
Choose a filtering function, set your desired parameters, and click 'Apply' to filter your table. The changes you make will not affect your original file until you save them.

.. image:: ../../rnalysis/gui/videos/filter_table.webp

Undo the operations you applied to your table
=============================================
At any moment, you can use the 'Command history' window to undo or redo an operation you applied to your table.

.. image:: ../../rnalysis/gui/videos/undo_actions.webp

Apply your operations 'in-place' or apply on a new table
========================================================
Instead of applying operations 'in-place', you can choose to apply the operation to a copy of your table in a new tab. The table in the original tab won't be modified.

.. image:: ../../rnalysis/gui/videos/apply_inplace.webp

Save the changes you made to your table
=======================================
To save result of your filtering operations, click the 'Save table' button and choose where to save the modified table.

.. image:: ../../rnalysis/gui/videos/save_table.webp

Work on multiple tables at the same time
========================================
You can work on multiple tables at the same time by opening a new tab and loading another table.

.. image:: ../../rnalysis/gui/videos/new_tab.webp

Different types of tables offer different ways to filter and analyze your data
==============================================================================
When loading a table, you can specify its type. Different types of tables support different types of functions: for example, count matrices support clustering analysis.

.. image:: ../../rnalysis/gui/videos/table_types.webp

Create and save graphs
======================
Some functions can generate graphs of your data. You can resize those graphs, and save them to your computer in multiple file formats.

.. image:: ../../rnalysis/gui/videos/generate_graphs.webp

Sort your tabs and change their icons
=====================================
To help organize your workspace, you can sort tabs by right-clicking a tab and choosing a sorting method. You can also change specific tabs' colors, to help you differentiate them.

.. image:: ../../rnalysis/gui/videos/sort_tabs.webp

Restore tabs you closed
=======================
If you accidentally closed one of your tabs - don't worry! You can restore closed tabs through the 'Edit' menu.

.. image:: ../../rnalysis/gui/videos/restore_tabs.webp

Import lists of genes as Gene Sets
==================================
In addition to tables, *RNAlysis* can also import lists of genes as Gene Sets. We will soon review what we can do with those gene sets.

.. image:: ../../rnalysis/gui/videos/import_gene_sets.webp

Visualize the intersections between your tables and gene sets
=============================================================
In the 'Visualize Gene Sets' window you can create Venn diagrams and UpSet plots that will display the various intersections between your tables and gene sets.

.. image:: ../../rnalysis/gui/videos/visualize_gene_sets.webp

Apply set operations to your tables and gene sets
=================================================
In the 'Set Operations' window you can extract specific subsets from your data. Either use predefined set operations, or click on specific subsets in the preview pane to select them.

.. image:: ../../rnalysis/gui/videos/set_operations.webp

Perform enrichment analysis on your tables and gene sets
========================================================
In the 'Enrichment Analysis' window, you can perform various types of enrichment analysis on the tables and gene sets you filtered.

.. image:: ../../rnalysis/gui/videos/enrichment_analysis.webp

Create Pipelines to streamline your data analysis
=================================================
You can group multiple operations in a specific order and with specific parameters into a Pipeline. Just add those functions to the Pipeline in the order you choose.

.. image:: ../../rnalysis/gui/videos/create_pipeline.webp

Apply Pipelines to one or more of your tables
=============================================
You can apply a Pipeline to a group of tables through the 'Pipelines' menu. Using Pipelines to analyze multiple datasets can make your workflow faster and less error-prone.

.. image:: ../../rnalysis/gui/videos/apply_pipeline.webp

Export and share Pipelines to make your analysis more reproducible
==================================================================
Pipelines you export can be imported from any computer, and can be shared with others to help make your analysis easier to understand and more reproducible.

.. image:: ../../rnalysis/gui/videos/export_pipeline.webp



****************************
*RNAlysis* filtering module
****************************
RNAlysis's filtering module (rnalysis.filtering) is built to allow rapid and easy to understand filtering of various forms of RNA sequencing data. The module also contains specific methods for visualization and clustering of data.

The filtering module is built around :term:`Filter objects`, which are containers for tabular sequencing data. You can use the different types of :term:`Filter objects` to apply filtering operations to various types of tabular data. You will learn more about :term:`Filter objects` in the next section.

Working with Filter objects
============================

All :term:`Filter objects` (:term:`CountFilter`, :term:`DESeqFilter`, :term:`Filter`, :term:`FoldChangeFilter`) work on the same principles,
and share many of the same functions and features. Each of them also has specific filtering, analysis and visualisation functions. In this section we will look into the general usage of :term:`Filter objects`.

Initialize a Filter object
--------------------------

We will start by importing the filtering module::

    >>> from *RNAlysis* import filtering

We can now, for example, create a :term:`DESeqFilter` object from a DESeq2 `csv` output file (see more details about :term:`DESeqFilter` in sections below).
::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")

View a Filter object
--------------------

In order to view a glimpse of the file we imported we can use the 'head' and 'tail' functions.
By default 'head' will show the first 5 rows of the file, and 'tail' will show the last 5 rows,
but you can specify a specific number of lines to show.
::

    >>> d.head()
                   baseMean  log2FoldChange  ...         pvalue           padj
    WBGene00000002  6820.755327        7.567762  ...   0.000000e+00   0.000000e+00
    WBGene00000003  3049.625670        9.138071  ...  4.660000e-302  4.280000e-298
    WBGene00000004  1432.911791        8.111737  ...  6.400000e-237  3.920000e-233
    WBGene00000005  4028.154186        6.534112  ...  1.700000e-228  7.800000e-225
    WBGene00000006  1230.585240        7.157428  ...  2.070000e-216  7.590000e-213
    <BLANKLINE>
    [5 rows x 6 columns]
    >>> d.tail(8)
                   baseMean  log2FoldChange  ...         pvalue           padj
    WBGene00000022   365.813048        6.101303  ...  2.740000e-97  2.400000e-94
    WBGene00000023  3168.566714        3.906719  ...  1.600000e-93  1.340000e-90
    WBGene00000024   221.925724        4.801676  ...  1.230000e-84  9.820000e-82
    WBGene00000025  2236.185837        2.477374  ...  1.910000e-81  1.460000e-78
    WBGene00000026   343.648987       -4.037191  ...  2.320000e-75  1.700000e-72
    WBGene00000027   175.142856        6.352044  ...  1.580000e-74  1.120000e-71
    WBGene00000028   219.163200        3.913657  ...  3.420000e-72  2.320000e-69
    WBGene00000029  1066.242402       -2.811281  ...  1.420000e-70  9.290000e-68
    <BLANKLINE>
    [8 rows x 6 columns]

We can also see the total number of rows and columns by accessing the 'shape' attribute::

    >>> d.shape
    (28, 6)

meaning there are 28 rows and 6 columns in the file.

Filtering operations
--------------------

Now we can start filtering the entries in the file according to parameters of our choosing.
Various filtering operations are applied directly to the :term:`Filter object`. Those operations do not affect the original `csv` file, but its representation within the :term:`Filter object`.
For example, we can the function 'filter_percentile' to remove all rows that are above the specified percentile (in our example, 75% percentile) in the specified column (in our example, 'log2FoldChange')::

    >>> d.filter_percentile(0.75,'log2FoldChange')
    Filtered 7 features, leaving 21 of the original 28 features. Filtered inplace.

If we now look at the shape of d, we will see that 5954 rows have been filtered out of the object, and we remain with 17781 rows.
::

    >>> d.shape
    (21, 6)

By default, filtering operations on :term:`Filter objects` are performed in-place, meaning the original object is modified. However, we can save the results into a new :term:`Filter object` and leave the current object unaffected by passing the argument 'inplace=False' to any filtering function within *RNAlysis*. For example::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")
    >>> d.shape
    (28, 6)
    >>> d_filtered = d.filter_percentile(0.75,'log2FoldChange',inplace=False)
    Filtered 7 features, leaving 21 of the original 28 features. Filtering result saved to new object.
    >>> d_filtered.shape
    (21, 6)
    >>> d.shape
    (28, 6)

In this case, the object 'd' remained unchanged, while 'd_filtered' became a new :term:`Filter object` which contains our filtered results. We can continue applying filters sequentially to the same Filter object, or using 'inplace=False' to create a new object at any point.

Another useful option is to perform an opposite filter. When we specify the parameter 'opposite=True' to any filtering function within *RNAlysis*, the filtering function will be performed in opposite. This means that all of the genomic features that were supposed to be filtered out are kept in the object, and the genomic features that were supposed to be kept in the object are filtered out.
For example, if we now wanted to remove the rows which are below the 25% percentile in the 'log2FoldChange' column, we will use the following code::

    >>> d.filter_percentile(0.25,'log2FoldChange',opposite=True)
    Filtered 7 features, leaving 21 of the original 28 features. Filtered inplace.

Calling this function without the 'opposite' parameter would have removed all values except the bottom 25% of the 'log2FoldChange' column. When specifying 'opposite', we instead throw out the bottom 25% of the 'log2FoldChange' column and keep the rest.

There are many different filtering functions within the filtering module. Some of them are subtype-specific (such as 'filter_low_reads' for :term:`CountFilter` objects and 'filter_significant' for :term:`DESeqFilter` objects), while others can be applied to any :term:`Filter object`. You can read more about the different functions and their usage in the project's documentation.


Performing set operations on multiple Filter objects
----------------------------------------------------

In addition to using regular filters, it is also possible to use set operations such as union, intersection, difference and symmetric difference to combine the results of multiple :term:`Filter objects`. Those set operations can be applied to any Filter object, as well as to python sets. The objects don't have to be of the same subtype - you can, for example, look at the union of a :term:`DESeqFilter` object, an :term:`CountFilter` object and a python set::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")
    >>> counts = filtering.CountFilter('tests/test_files/counted.csv')
    >>> a_set = {'WBGene00000001','WBGene00000002','WBGene00000003'}
    >>> d.difference(counts, a_set)
    {'WBGene00007063', 'WBGene00007064', 'WBGene00007066', 'WBGene00007067', 'WBGene00007069', 'WBGene00007071',
     'WBGene00007074', 'WBGene00007075', 'WBGene00007076', 'WBGene00007077', 'WBGene00007078', 'WBGene00007079',
     'WBGene00014997', 'WBGene00043987', 'WBGene00043988', 'WBGene00043989', 'WBGene00043990', 'WBGene00044022',
     'WBGene00044951', 'WBGene00077502', 'WBGene00077503', 'WBGene00077504'}

When performing set operations, the return type can be either a python set (default) or a string. This means you can use the output of the set operation as an input for yet another set operation. However, since the returned object is a set you cannot use :term:`Filter object` functions such as 'head' and 'save_csv' on it, or apply filters to it directly. Intersection and Difference in particular can be used in-place, which applies the filtering to the first :term:`Filter object`.


Saving Filter results
---------------------

At any point we can save the current result of our filtering to a new `csv` file, by using the 'save_csv' function::

    >>> d.save_csv()

If no filename is specified, the file is given a name automatically based on the filtering operations performed on it, their order and their parameters.
We can view the current automatic filename by looking at the 'fname' attribute::

    >>> d.filter_percentile(0.75,'log2FoldChange')
    Filtered 7 features, leaving 21 of the original 28 features. Filtered inplace.
    >>> d.number_filters('baseMean','greater than',500)
    Filtered 6 features, leaving 15 of the original 21 features. Filtered inplace.
    >>> d.fname
    'D:/myfolder/test_deseq_below0.75baseMeangt500.csv'

Alternatively, you can specify a filename::

    >>> d.save_csv('alternative_filename')

Instead of directly saving the results to a file, you can also get them as a set or string of genomic feature indices::

    >>> print(d.index_set)
    {'WBGene00000005', 'WBGene00000006', 'WBGene00000008', 'WBGene00000009', 'WBGene00000010', 'WBGene00000011',
     'WBGene00000012', 'WBGene00000014', 'WBGene00000015', 'WBGene00000017', 'WBGene00000019', 'WBGene00000021',
     'WBGene00000023', 'WBGene00000025', 'WBGene00000029'}
    >>> print(d.index_string)
    WBGene00000010
    WBGene00000012
    WBGene00000021
    WBGene00000023
    WBGene00000017
    WBGene00000015
    WBGene00000025
    WBGene00000008
    WBGene00000011
    WBGene00000014
    WBGene00000029
    WBGene00000006
    WBGene00000009
    WBGene00000005
    WBGene00000019

Sets of genomic feature indices can be used later for enrichment analysis using the enrichment module (see below).


Using an Attribute Reference Table for filter operations
---------------------------------------------------------

An :term:`Attribute Reference Table` contains various user-defined attributes (such as 'genes expressed in intestine', 'epigenetic genes' or 'genes that have paralogs') and their value for each genomic feature.
You can read more about the :term:`Attribute Reference Table` format and loading an :term:`Attribute Reference Table` in the :ref:`reference-table-ref` section.
Using the function Filter.filter_by_attribute(), you can filter your genomic features by one of the user-defined attributes in the Reference Table::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")
    >>> d.filter_by_attribute('attribute1', ref='tests/test_files/attr_ref_table_for_examples.csv')
    Filtered 27 features, leaving 1 of the original 28 features. Filtered inplace.

Using a Biotype Reference Table for filter operations
--------------------------------------------------------

A :term:`Biotype Reference Table` contains annotations of the biotype of each genomic features ('protein_coding', 'piRNAs', 'lincRNAs', 'pseudogenes', etc).
You can read more about the :term:`Biotype Reference Table` format and loading a :term:`Biotype Reference Table` in the :ref:`reference-table-ref` section.
Using the function Filter.filter_biotype_from_ref_table(), you can filter your genomic features by their annotated biotype in the Biotype Reference Table::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")
    >>> d.filter_biotype_from_ref_table('protein_coding', ref='tests/test_files/biotype_ref_table_for_tests.csv')
    Filtered 2 features, leaving 26 of the original 28 features. Filtered inplace.

You can also view the number of genomic features belonging to each biotype using the function Filter.biotypes_from_ref_table()::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")
    >>> d.biotypes_from_ref_table()
                    gene
    biotype
    protein_coding    26
    pseudogene         1
    unknown            1

Or view more elaborated descriptive statistics for eahc biotype by specifying return_format='long'::

    >>> d.biotypes_from_ref_table(return_format='long', ref='tests/test_files/biotype_ref_table_for_tests.csv')

                   baseMean               ...           padj
                      count         mean  ...            75%            max
    biotype                               ...
    protein_coding     26.0  1823.089609  ...   1.005060e-90   9.290000e-68
    pseudogene          1.0  2688.043701  ...   1.800000e-94   1.800000e-94
    unknown             1.0  2085.995094  ...  3.070000e-152  3.070000e-152
    <BLANKLINE>
    [3 rows x 48 columns]


Filtering DESeq2 output files with filtering.DESeqFilter
=========================================================

:term:`DESeqFilter` objects are built to easily filter differential expression tables, such as those returned by the R package DESeq2.
Like other Filter Objects, filtering operations on :term:`DESeqFilter` are performed in-place by default, meaning the original object is modified.

You can read more about DESeq2 here:
https://doi.org/doi:10.18129/B9.bioc.DESeq2

Loading from a `csv` file
----------------------------

Any `csv` file that contains differential expression analysis data with log2 fold change and adjusted p-values can be used as input for :term:`DESeqFilter`.
By default, *RNAlysis* assumes that log2 fold change values will be specified under a 'log2FoldChange' column, and adjusted p-values will be specified under a 'padj' column (as is the default in differential expression tables generated by DESeq2):

+----------------+----------+----------------+----------+----------+----------+----------+
|                | baseMean | log2FoldChange | lfcSE    | stat     | pvalue   | padj     |
+================+==========+================+==========+==========+==========+==========+
| WBGene00000021 | 2688.044 | 3.329404       | 0.158938 | 20.94783 | 1.96E-97 | 1.80E-94 |
+----------------+----------+----------------+----------+----------+----------+----------+
| WBGene00000022 | 365.813  | 6.101303       | 0.291484 | 20.93189 | 2.74E-97 | 2.40E-94 |
+----------------+----------+----------------+----------+----------+----------+----------+
| WBGene00000023 | 3168.567 | 3.906719       | 0.190439 | 20.51433 | 1.60E-93 | 1.34E-90 |
+----------------+----------+----------------+----------+----------+----------+----------+
| WBGene00000024 | 221.9257 | 4.801676       | 0.246313 | 19.49419 | 1.23E-84 | 9.82E-82 |
+----------------+----------+----------------+----------+----------+----------+----------+
| WBGene00000025 | 2236.186 | 2.477374       | 0.129606 | 19.11463 | 1.91E-81 | 1.46E-78 |
+----------------+----------+----------------+----------+----------+----------+----------+
| WBGene00000026 | 343.649  | -4.03719       | 0.219781 | -18.3691 | 2.32E-75 | 1.70E-72 |
+----------------+----------+----------------+----------+----------+----------+----------+
| WBGene00000027 | 175.1429 | 6.352044       | 0.347777 | 18.26471 | 1.58E-74 | 1.12E-71 |
+----------------+----------+----------------+----------+----------+----------+----------+
| WBGene00000028 | 219.1632 | 3.913657       | 0.217802 | 17.96885 | 3.42E-72 | 2.32E-69 |
+----------------+----------+----------------+----------+----------+----------+----------+

Loading a file that follows this format into a :term:`DESeqFilter` works similarly to other Filter objects::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")


If your differential expression table does not follow this format, you can specify the exact names of the columns in your table that contain log2 fold change values and adjusted p-values::

    >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv", log2fc_col='name of log2 fold change column', padj_col='name of adjusted p-value column')



Unique :term:`DESeqFilter` functions (such as 'filter_significant' and 'filter_abs_log2_fold_change') will not execute properly if the log2 fold change column and adjusted p-value column are not defined correctly.

Filtering operations unique to DESeqFilter
------------------------------------------

There are a few filtering operations unique to DESeqFilter. Those include 'filter_significant', which removes statistically-insignificant rows according to a specified threshold; 'filter_abs_log2_fold_change', removes rows whose absolute value log2 fold change is below the specified threshold; 'filter_fold_change_direction' which removes either up-regulated (positive log2 fold change) or down-regulated (negative log2 fold change) rows; and 'split_fold_change_direction' which returns a :term:`DESeqFilter` object with only up-regulated features and a :term:`DESeqFilter` object with only down-regulated features.


Data visualization and exploratory data analysis with DESeqFilter
------------------------------------------------------------------------
:term:`DESeqFilter` supports methods for visualization and exploratory analysis of differential expression data.


With DESeqFilter.volcano_plot, you can observe the direction, magnitude, and significance of differential expression within your data:

.. figure:: /figures/volcano.png
           :align:   center
           :scale: 70 %

           Example output of DESeqFilter.volcano_plot()


Filtering count matrices with filtering.CountFilter
===============================================================

:term:`CountFilter` objects are capable of visualizing, filtering, normalizing, and clustering of read count matrices (the output of feature-counting software such as HTSeq-count and featureCounts).
Data can be imported into a CountFilter objects either from a `csv` file, or directly from HTSeq-count output files as explained below.

You can read more about HTSeq-count here:
https://htseq.readthedocs.io/en/master/count.html

In principle, any `csv` file where the columns are different conditions/replicates and the rows include reads/normalized reads per genomic feature can be used as input for CountFilter. However, some :term:`CountFilter` functions (such as 'normalize_to_rpm_htseqcount') will only work on HTSeq-count output files, and other unintended interactions may occur.

.. _from-folder-ref:

Generating an CountFilter object from a folder of HTSeq-count output .txt files
---------------------------------------------------------------------------------
HTSeq-count receives as input an aligned SAM/BAM file. The native output of HTSeq-count is a text file with feature indices and read-per-genomic-feature, as well as information about reads that weren't counted for any feature (alignment not unique, low alignment quality, ambiguous, unaligned, aligned to no feature).
An HTSeq-count output file would follow the following format:

+------------------------+-----+
| WBGene00000001         | 376 |
+------------------------+-----+
| WBGene00000002         | 1   |
+------------------------+-----+
| WBGene00000003         | 1   |
+------------------------+-----+
| WBGene00000004         | 18  |
+------------------------+-----+
| WBGene00000005         | 1   |
+------------------------+-----+
| WBGene00000006         | 3   |
+------------------------+-----+
| WBGene00000007         | 6   |
+------------------------+-----+
| WBGene00000008         | 0   |
+------------------------+-----+
| WBGene00000009         | 1   |
+------------------------+-----+
| WBGene00000010         | 177 |
+------------------------+-----+
| __no_feature           | 32  |
+------------------------+-----+
| __ambiguous            | 12  |
+------------------------+-----+
| __too_low_aQual        | 1   |
+------------------------+-----+
| __not_aligned          | 121 |
+------------------------+-----+
| __alignment_not_unique | 100 |
+------------------------+-----+

When running HTSeq-count on multiple SAM files (which could represent different conditions or replicates), the final output would be a directory of .txt files. *RNAlysis* can parse those .txt files into two `csv` tables: in the first each row is a genomic feature and each column is a condition or replicate (a single .txt file), and in the second each row represents a category of reads not mapped to genomic features (alignment not unique, low alignment quality, etc). This is done with the 'from_folder' function::

    >>> counts = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder')

By deault, 'from_folder' does not save the generated tables as `csv` files. However, you can choose to save them by specifying 'save_csv=True', and specifying their filenames in the arguments 'counted_fname' and 'uncounted_fname'::

    >>> counts = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder', save_csv=True, counted_fname='name_for_reads_csv_file', uncounted_fname='name_for_uncounted_reads_csv_file')

It is also possible to automatically normalize the reads in the new :term:`CountFilter` object to reads per million (RPM) using the unmapped reads data by specifying 'norm_to_rpm=True'::

        >>> counts = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder', norm_to_rpm=True)


Loading from a pre-made `csv` file
----------------------------------
If you have previously generated a `csv` file from HTSeq-count output files using *RNAlysis*, or have done so manually, you can directly load this `csv` file into an :term:`CountFilter` object as you would any other Filter object::

    >>> counts = filtering.CountFilter('tests/test_files/counted.csv')

A correct input to a :term:`CountFilter` object would follow the following format:

+----------------+-------+-------+-------+-------+
|                | cond1 | cond2 | cond3 | cond4 |
+================+=======+=======+=======+=======+
| WBGene00007063 | 633   | 451   | 365   | 388   |
+----------------+-------+-------+-------+-------+
| WBGene00007064 | 60    | 57    | 20    | 23    |
+----------------+-------+-------+-------+-------+
| WBGene00044951 | 0     | 0     | 0     | 1     |
+----------------+-------+-------+-------+-------+
| WBGene00007066 | 55    | 266   | 46    | 39    |
+----------------+-------+-------+-------+-------+
| WBGene00007067 | 15    | 13    | 1     | 0     |
+----------------+-------+-------+-------+-------+
| WBGene00007069 | 0     | 2     | 1     | 0     |
+----------------+-------+-------+-------+-------+
| WBGene00077502 | 0     | 0     | 0     | 0     |
+----------------+-------+-------+-------+-------+
| WBGene00077503 | 1     | 4     | 2     | 0     |
+----------------+-------+-------+-------+-------+
| WBGene00077504 | 0     | 0     | 0     | 0     |
+----------------+-------+-------+-------+-------+

Filtering operations unique to CountFilter
--------------------------------------------
There are a few filtering operations unique to CountFilter. Those include 'filter_low_reads', which removes rows that have less than n reads in all columns.

Normalizing reads with CountFilter
------------------------------------
:term:`CountFilter` offers two methods for normalizing reads: normalize with one of the pre-made normalization methods *RNAlysis* supplies, or using user-defined scaling factors. Data normalized in other methods (such as RPKM) can be used as input for CountFilter as well.

*RNAlysis* supplies the following normalization methods:

* Relative Log Expression (RLE - 'normalize_rle'), used by default by R's DESeq2
* Trimmed Mean of M-values (TMM - 'normalize_tmm'), used by default by R's edgeR
* Quantile normalization, a generalization of Upper Quantile normalization (UQ - 'normalize_quantile'), used by default by R's Limma
* Median of Ratios Normalization (MRN - 'normalize_mrn')
* Reads Per Million (RPM - 'normalize_to_rpm')

To normalize a :term:`CountFilter` with one of these functions, simply call the function with your preferred parameters, if there are any. For example::

    >>> counts = filtering.CountFilter('tests/test_files/counted.csv')
    >>> counts.normalize_rle()
    Normalized the values of 22 features. Normalized inplace.

To normalize a :term:`CountFilter` with user-generated scaling factors, we need a `csv` table with the size factor for each sample:

+----------------+----------------+----------------+----------------+
|    sample1     |    sample2     |    sample3     |    sample4     |
+================+================+================+================+
|      0.96      |       1        |      0.78      |      1.23      |
+----------------+----------------+----------------+----------------+

We would then supply the function with the path to the scaling factors file::

    >>> counts = filtering.CountFilter('tests/test_files/counted.csv')
    >>> counts.normalize_with_scaling_factors('scaling_factors.csv')
    Normalized the values of 22 features. Normalized inplace.

The resulting :term:`CountFilter` object will be normalized with the scaling factors (dividing the value of each column by the value of the corresponding scaling factor).


To normalize a :term:`CountFilter` that originated from HTSeq-count to reads per million, we need a `csv` table with the special counters that appear in HTSeq-count output:

+------------------------+---------+---------+---------+---------+
|                        | sample1 | sample2 | sample3 | sample4 |
+========================+=========+=========+=========+=========+
| __ambiguous            | 37      | 12      | 145     | 77      |
+------------------------+---------+---------+---------+---------+
| __no_feature           | 9468    | 11354   | 14009   | 30287   |
+------------------------+---------+---------+---------+---------+
| __alignment_not_unique | 108     | 290     | 557     | 893     |
+------------------------+---------+---------+---------+---------+
| __too_low_aQual        | 0       | 5       | 12      | 9       |
+------------------------+---------+---------+---------+---------+
| __not_aligned          | 109853  | 277653  | 88653   | 96012   |
+------------------------+---------+---------+---------+---------+

Such a `csv` table is generated automatically when you create a :term:`CountFilter` object from a folder of text files (CountFilter.from_folder(), see :ref:`from-folder-ref`).
We would then supply the normalization function with the path to the special counter file::

    >>> counts = CountFilter("tests/test_files/counted.csv")
    >>> counts.normalize_to_rpm_htseqcount("tests/test_files/uncounted.csv")
    Normalized the values of 22 features. Normalized inplace.

The resulting :term:`CountFilter` object will be normalized to RPM with the formula (1,000,000 * reads in cell) / (sum of aligned reads + __no_feature + __ambiguous + __alignment_no_unique)


Data clustering with CountFilter
----------------------------------
*RNAlysis* supports a wide variety of clustering methods, which can group genomic features into clusters according to their similarity across different samples.

When clustering genomic features in a :term:`CountFilter` object, the called function returns a tuple of :term:`CountFilter` objects, with each object corresponding to one cluster of genomic features.

Expression plots of the resulting clusters can be generated in one of multiple styles:

 .. figure:: /figures/kmeans_all.png
           :align:   center
           :scale: 40 %

           Example expression plot of clustering results with plot_style='all'

 .. figure:: /figures/kmeans_std_area.png
           :align:   center
           :scale: 40 %

           Example expression plot of clustering results with plot_style='std_area'

 .. figure:: /figures/kmeans_std_bar.png
           :align:   center
           :scale: 40 %

           Example expression plot of clustering results with plot_style='std_bar'

 .. figure:: /figures/clustering_PCA_clicom.png
           :align:   center
           :scale: 40 %

           Example PCA plot of clustering results

The expression plots can also by split into separate graphs, one for each discovered cluster.

All clustering methods in *RNAlysis* which require you to specify the expected number of clusters (such as K in K-Means clustering) allow multiple ways of specifying the number of clusters you want to find.
You can specify a single value::

    >>> counts = CountFilter("tests/test_files/counted.csv")
    >>> five_clusters = counts.split_kmeans(n_clusters=5)
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 7 features, leaving 15 of the original 22 features. Filtering result saved to new object.
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 21 features, leaving 1 of the original 22 features. Filtering result saved to new object.
    >>> print(len(five_clusters))
    5

You can specify a list of values to be used, and each value will be calculated and returned separately::

    >>> counts = CountFilter("tests/test_files/counted.csv")
    >>> five_clusters, two_clusters = counts.split_kmeans(n_clusters=[5,2])
    Filtered 21 features, leaving 1 of the original 22 features. Filtering result saved to new object.
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 7 features, leaving 15 of the original 22 features. Filtering result saved to new object.
    Filtered 4 features, leaving 18 of the original 22 features. Filtering result saved to new object.
    Filtered 18 features, leaving 4 of the original 22 features. Filtering result saved to new object.
    >>> print(len(five_clusters))
    5
    >>> print(len(two_clusters))
    2

Finally, you can use a model selection method to estimate the number of clusters in your dataset. *RNAlysis* supports both the Silhouette method and the Gap Statistic method::

    >>> counts = CountFilter("tests/test_files/counted.csv")
    >>> silhouette_clusters = counts.split_kmeans(n_clusters='silhouette')
    Estimating the optimal number of clusters using the Silhouette method in range 2:5...
    Using the Silhouette method, 4 was chosen as the best number of clusters (k).
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 6 features, leaving 16 of the original 22 features. Filtering result saved to new object.
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.
    >>> print(len(silhouette_clusters))
    4
    >>> gap_stat_clusters = counts.split_kmeans(n_clusters='gap')
    Estimating the optimal number of clusters using the Gap Statistic method in range 2:5...
    Using the Gap Statistic method, 2 was chosen as the best number of clusters (K).
    Filtered 4 features, leaving 18 of the original 22 features. Filtering result saved to new object.
    Filtered 18 features, leaving 4 of the original 22 features. Filtering result saved to new object.
    >>> print(len(gap_stat_clusters))
    2

To help in evaluating the result of these model selection methods, *RNAlysis* will also plot a summary of their outcome:

.. image:: /figures/ gap_statistic.png
           :width: 60 %
.. image:: /figures/ silhouette.png
           :width: 30 %

|

K-Means clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^
K-means is a clustering method which partitions all of the data points into K clusters by minimizing the squared eucliean distance between points within each cluster.

The algorithm is initiated by picking a random starting point, and therefore the exact clustering results can change between runs.

The main advantage of K-means clustering is its simplicity - it contains one main tuning parameter (*K*, the expected number of clusters in the data).

.. image:: /figures/kmeans_all.png
  :width: 450
  :alt: K-means clustering output figure

|

K-Medoids clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The K-medoids method is very similar to K-means. The main difference between the two is the way they define clusters and the distances between them:
K-medoids picks one data point as the 'center' (medoid) of each cluster.
In addition, K-medoids attempts to minimize the sum of dissimilarities within each cluster, instead of minimizing squared euclidean distance.

Due to these differences, the K-medoids algorithm supports the use of distance metrics other than eucliean distance through the `metric` parameter.

K-medoids clustering in *RNAlysis* supports the following distance metrics:

* eucliidean
* cosine
* pearson
* spearman
* manhattan (cityblock)
* l1
* l2
* jackknife (see `Heyer, Kruglyak and Yooseph 1999 <https://doi.org/10.1101%2Fgr.9.11.1106>`_)
* YS1 (see `Son and Baek 2007 <https://doi.org/10.1016/j.patrec.2007.09.015>`_)
* YR1 (see `Son and Baek 2007 <https://doi.org/10.1016/j.patrec.2007.09.015>`_)
* hammming
* all other pairwise distance metrics supported by scikit-learn

.. image:: /figures/kmedoids_all.png
  :width: 450
  :alt: K-medoids clustering output figure

|

Hierarchical clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Hierarchical clustering (or agglomerative clustering) is clustering method which aims to build a hierarchy of clusters.

In agglomerative hierarchical clustering, each data points starts in its own clusters.
The clustering algorithm then uses a distance metric (a measure of distance between pairs of data points)
and a linkage criterion
(determines the distance between sets of data points as a function of the pairwise distances between observations)
to group merge data points into clusters, and then further group those clusters into larger clusters based on their similarity.
Eventually, all of the observations are connected into a hierarchical tree.

We can decide to 'cut' the tree at any height in order to generate the final clustering solution.
This can be done by either specifying the estimated number of clusters like in K-means,
or by specifiying a distance threshold above which clusters will not be merged.

Hierarchical clustering in *RNAlysis* supports the following distance metrics:

* euclidean
* cosine
* pearson
* spearman
* manhattan (cityblock)
* l1
* l2
* jackknife (see `Heyer, Kruglyak and Yooseph 1999 <https://doi.org/10.1101%2Fgr.9.11.1106>`_)
* YS1 (see `Son and Baek 2007 <https://doi.org/10.1016/j.patrec.2007.09.015>`_)
* YR1 (see `Son and Baek 2007 <https://doi.org/10.1016/j.patrec.2007.09.015>`_)


.. image:: /figures/hierarchical_all.png
  :width: 450
  :alt: Hierarchical clustering output figure

|

HDBSCAN clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
HDBSCAN stands for Hierarchical Density-Based Spatial Clustering of Applications with Noise (see https://link.springer.com/chapter/10.1007/978-3-642-37456-2_14 ).
HDBSCAN offers multiple advantages over more traditional clustering methods:

1. HSBSCAN makes relatively few assumptions about the data - it assumes that the data contains noise, as well as some real clusters which we hope to discover.
2. Unlike most other clustering methods, HDBSCAN does not force every data point to belong to a cluster. Instead, it can classify data points as outliers, excluding them from the final clustering solution.
3. HDBSCAN does not require you to guess the number of clusters in the data. The main tuning parameter in HDBSCAN is *minimum cluster size* (`min_cluster_size`), which determines the smallest "interesting" cluster size we expect to find in the data.

HDBSCAN supports additional tuning parameters, which you can read more about in the `HDBSCAN documentation <https://hdbscan.readthedocs.io/en/latest/parameter_selection.html>`_:

HDBSCAN in *RNAlysis* supports the following distance metrics:

* eucliidean
* cosine
* pearson
* spearman
* manhattan (cityblock)
* l1
* l2
* jackknife (see `Heyer, Kruglyak and Yooseph 1999 <https://doi.org/10.1101%2Fgr.9.11.1106>`_)
* YS1 (see `Son and Baek 2007 <https://doi.org/10.1016/j.patrec.2007.09.015>`_)
* YR1 (see `Son and Baek 2007 <https://doi.org/10.1016/j.patrec.2007.09.015>`_)
* hammming
* all other pairwise distance metrics elaborated in the `HDBSCAN documentation <https://hdbscan.readthedocs.io/en/latest/basic_hdbscan.html?#what-about-different-metrics>`_.

.. image:: /figures/hdbscan_all.png
  :width: 450
  :alt: HDBSCAN output figure

|

CLICOM clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CLICOM is an ensemble-based clustering algorithm (see https://doi.org/10.1016/j.eswa.2011.08.059 ).
The CLICOM algorithm incorporates the results of multiple clustering solutions, which can come from different clustering algorithms with differing clustering parameters, and uses these clustering solutions to create a combined clustering solution.
CLICOM offers multiple advantages over more traditional clustering methods:

1. The ensemble clustering approach allows you to combine the results of multiple clustering algorithms with multiple tuning parameters, potentially making up for the weaknesses of each individual clustering method, and only taking into account patterns that robustly appear in many clustering solutions.
2. Unlike most other clustering methods, CLICOM does not have to force every data point to belong to a cluster. Instead, it can classify data points as outliers, excluding them from the final clustering solution.
3. CLICOM does not require you to guess the final number of clusters in the data. The main tuning parameter in HDBSCAN is the *evidence threshold* (`evidence_threshold`).

*RNAlysis* offers a modified implementation of CLICOM. This implementation of CLICOM supports a few tuning parameters, in addition to the clustering solutions themselves:
Moreover, ths modified version of the algorithm can cluster each batch of biological/technical replicates in your data separately, which can reduce the influence of batch effect on clustering results, and increases the accuracy and robustness of your clustering results.

* `evidence_threshold`: a higher evidence threshold leads to fewer, large clusters, with fewer features being classified as outliers.
* `cluster_unclustered_features`: if True, CLICOM will force every feature to belong to a discovered cluster. Otherwise, features can be classified as noise and remain unclustered.
* `min_cluster_size`: determines the minimal size of a cluster you would consider meaningful. Clusters smaller than this would be classified as noise and filtered out of the final result, or merged into other clusters (depending on the value of `cluster_unclustered_features`).
* `replicates_grouping`: allows you to group samples into technical/biological batches. The algorithm will then cluster each batch of samples separately, and use the CLICOM algorithm to find an ensemble clustering result from all of the separate clustering results.


.. image:: /figures/clicom_all.png
  :width: 450
  :alt: CLICOM output figure

|

Data visualization and exploratory data analysis with CountFilter
------------------------------------------------------------------------
:term:`CountFilter` includes multiple methods for visualization and exploratory analysis of count data.


With CountFilter.pairplot, you can get a quick overview of the distribution of counts within each sample, and the correlation between different samples:

.. figure:: /figures/pairplot.png
           :align:   center
           :scale: 40 %

           Example output of CountFilter.pairplot()

With CountFilter.pca, you can perform a principal component analysis and look for strong patterns in your dataset:

 .. figure:: /figures/pca.png
           :align:   center
           :scale: 40 %

           Example plot of CountFilter.pca()

With CountFilter.plot_expression, you can examine the average expression of specific genomic features under the specific conditions:

 .. figure:: /figures/plot_expression.png
           :align:   center
           :scale: 60 %

           Example plot of CountFilter.plot_expression()

With CountFilter.clustergram, you can cluster your samples according to specified distance and linkage metrics:

 .. figure:: /figures/clustergram.png
           :align:   center
           :scale: 40 %

           Example plot of CountFilter.clustergram()

Filtering fold-change data of features using filtering.FoldChangeFilter
=======================================================================

:term:`FoldChangeFilter` objects can perform filtering operations and randomization tests on fold change values between two conditions.

A :term:`FoldChangeFilter` object can be calculated from a :term:`CountFilter` object (you can read more about it in the :ref:`fold-change-from-count-ref`), or imported from a `csv` file like other :term:`Filter objects`.

.. warning:: by default, :term:`FoldChangeFilter` assumes that fold change is calculated as (numerator_reads+1)/(denominator_reads+1), and does not support 0 and inf values. If you load a `csv` file which contains 0 and/or inf values into a :term:`FoldChangeFilter` object, unintended results and interactions may occur.

Unlike other Filter object, the underlying data structure storing the values is a pandas Series and not a pandas DataFrame, and lacks the Columns attribute.

Loading fold change data from a `csv` file
------------------------------------------------

Like with other objects from the Filter family, you can simply load a pre-existing or pre-calculated `csv` file into a :term:`FoldChangeFilter` object. However, in addition to the file path you will also have to enter the name of the numerator condition and the name of the denominator condition::

    >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv','name of numerator condition', 'name of denominator condition')

The names of the conditions are saved in the object attributes 'numerator' and 'denominator'::

    >>> f.numerator
    'name of numerator condition'
    >>> f.denominator
    'name of denominator condition'

.. warning:: by default, :term:`FoldChangeFilter` assumes that fold change is calculated as (mean_numerator_reads+1)/(mean_denominator_reads+1), and does not support 0 and inf values. If you load a `csv` file which contains 0 and/or inf values into a :term:`FoldChangeFilter` object, unintended results and interactions may occur.

.. _fold-change-from-count-ref:

Generating fold change data from an existing CountFilter object
-----------------------------------------------------------------

Alternatively, you can generate a :term:`FoldChangeFilter` object from count data in a :term:`CountFilter` object. We will start by loading a :term:`CountFilter` object::

    >>> counts = filtering.CountFilter('tests/test_files/counted_fold_change.csv')

The :term:`CountFilter` has the following columns::

    >>> counts.columns
    ['cond1_rep1', 'cond1_rep2', 'cond2_rep1', 'cond2_rep2', 'cond3_rep1', 'cond3_rep2']

We will now calculate the fold change between the mean of condition1 and condition2. Fold change is calculated as (mean_numerator_reads+1)/(mean_denominator_reads+1). We will need to specify the numerator columns, the denominator columns, and the names of the numerator and denominator. Specifying names is optional - if no names are specified, they will be generator automatically from columns used as numerator and denominator. Since we have multiple replicates of each condition, we will specify all of them in a list::

    >>> f = counts.fold_change(['cond1_rep1','cond1_rep2'],['cond2_rep1','cond2_rep2'])

In this example we did not specify names for the numerator and denominator, and therefore they were generated automatically::

    >>> f.numerator
    "Mean of ['cond1_rep1', 'cond1_rep2']"
    >>> f.denominator
    "Mean of ['cond2_rep1', 'cond2_rep2']"

We now have a :term:`FoldChangeFilter` object that we can perform further filtering operations on.

Performing randomization tests on a FoldChangeFilter object
------------------------------------------------------------

You can perform a randomization test to examine whether the fold change of a group of specific genomic features (for example, genes with a specific biological function) is significantly different than the fold change of a background set of genomic features.
To perform a randomization test you need two :term:`FoldChangeFilter` objects: one which contains the fold change values of all background genes, and another which contains the fold change values of your specific group of interest. For example::

    >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv' , 'numerator' , 'denominator')
    >>> f_background = f.filter_biotype_from_ref_table('protein_coding', ref='tests/test_files/biotype_ref_table_for_tests.csv', inplace=False) #keep only protein-coding genes as reference
    Filtered 9 features, leaving 13 of the original 22 features. Filtering result saved to new object.
    >>> f_test = f_background.filter_by_attribute('attribute1', ref='tests/test_files/attr_ref_table_for_examples.csv', inplace=False)
    Filtered 6 features, leaving 7 of the original 13 features. Filtering result saved to new object.
    >>> rand_test_res = f_test.randomization_test(f_background)
    Calculating...
       group size  observed fold change  ...      pval  significant
    0           7              2.806873  ...  0.360264        False
    <BLANKLINE>
    [1 rows x 5 columns]

The output table would look like this:

+------------+----------------------+----------------------+--------+-------------+
| group size | observed fold change | expected fold change | pval   | significant |
+============+======================+======================+========+=============+
|   7        |       2.806873       |  2.51828             |0.36026 | False       |
+------------+----------------------+----------------------+--------+-------------+

Sequentially applying filtering operations using Pipelines
============================================================
:term:`Pipeline` objects allow you to group together multiple operations from the *filtering* module (such as filtering, splitting, normalizing, plotting or describing your data), and apply this group of operations to :term:`Filter objects` of your choice in a specific and consistent order.
Pipelines make your workflow easier to read and understand, help you avoid repetitive code, and makes your analyses more reproducible and less error-prone.

Creating a new Pipeline
________________________
To create a new empty :term:`Pipeline`, simply create a new Pipeline object::

    >>> from *RNAlysis* import filtering
    >>> pipe = Pipeline()

Because every :term:`Filter object` has its own unique functions, a particular Pipeline can only contain functions of a specific Filter object type, and can only be applied to objects of that type.
By default, a new Pipeline's `filter_type` is :term:`Filter`, and can only contain general functions from the *filtering* module that can apply to any Filter object.
If we wanted, for example, to create a Pipeline for DESeqFilter objects, we would have to specify the parameter `filter_type`::

    >>> from *RNAlysis* import filtering
    >>> deseq_pipe = filtering.Pipeline('deseqfilter')

One we have an empty :term:`Pipeline`, we can start adding functions to it.
We can do that either via the function's name::

    >>> from *RNAlysis* import filtering
    >>> pipe = filtering.Pipeline('DESeqFilter')
    >>> pipe.add_function('filter_significant')
    Added function 'DESeqFilter.filter_significant()' to the pipeline.

or via the function itself::

    >>> from *RNAlysis* import filtering
    >>> pipe = filtering.Pipeline('DESeqFilter')
    >>> pipe.add_function(filtering.DESeqFilter.filter_significant)
    Added function 'DESeqFilter.filter_significant()' to the pipeline.

We can also specify the function's arguments. We can specify both non-keyworded and keyworded arguments, just as we would if we called the function normally::

    >>> from *RNAlysis* import filtering
    >>> pipe = filtering.Pipeline()
    >>> pipe.add_function(filtering.Filter.filter_biotype_from_ref_table, biotype='protein_coding')
    Added function 'Filter.filter_biotype_from_ref_table(biotype='protein_coding')' to the pipeline.
    >>> pipe.add_function('number_filters', 'column1', 'gt', value=5, opposite=True)
    Added function 'Filter.number_filters('column1', 'gt', value=5, opposite=True)' to the pipeline.

We can also view the functions currently in the Pipeline object, their arguments, and their order::

    >>> print(pipe)
    Pipeline for Filter objects:
        Filter.filter_biotype_from_ref_table(biotype='protein_coding')
        Filter.number_filters('column1', 'gt', value=5, opposite=True)
    >>> print(repr(pipe))
    Pipeline('Filter'): Filter.filter_biotype_from_ref_table(biotype='protein_coding')-->Filter.number_filters('column1', 'gt', value=5, opposite=True)


We can also remove functions from the Pipeline::

    >>> pipe.remove_last_function()
    Removed function number_filters with parameters ['column1', 'gt', value=5, opposite=True] from the pipeline.

Now that we have a Pipeline with multiple functions, we can apply it to our Filter objects.

Applying Pipelines to Filter objects
_____________________________________
Just like with other functions in the *filtering* module, the functions in a :term:`Pipeline` can be applied either inplace or returned as a new object.
You can determine that via the `inplace` argument of the function `Pipeline.apply_to()`::

    >>> from *RNAlysis* import filtering
    >>> # create the pipeline
    >>> pipe = filtering.Pipeline('DESeqFilter')
    >>> pipe.add_function(filtering.DESeqFilter.filter_missing_values)
    Added function 'DESeqFilter.filter_missing_values()' to the pipeline.
    >>> pipe.add_function(filtering.DESeqFilter.filter_top_n, by='padj', n=3)
    Added function 'DESeqFilter.filter_top_n(by='padj', n=3)' to the pipeline.
    >>> pipe.add_function('sort', by='baseMean')
    Added function 'DESeqFilter.sort(by='baseMean')' to the pipeline.
    >>> # load the Filter object
    >>> d = filtering.DESeqFilter('tests/test_files/test_files/test_deseq_with_nan.csv')
    >>> # apply the Pipeline not-inplace
    >>> d_filtered = pipe.apply_to(d, inplace=False)
    Filtered 3 features, leaving 25 of the original 28 features. Filtering result saved to new object.
    Filtered 22 features, leaving 3 of the original 25 features. Filtering result saved to new object.
    Sorted 3 features. Sorting result saved to a new object.
    >>> # apply the Pipeline inplace
    >>> pipe.apply_to(d)
    Filtered 3 features, leaving 25 of the original 28 features. Filtered inplace.
    Filtered 22 features, leaving 3 of the original 25 features. Filtered inplace.
    Sorted 3 features. Sorted inplace.

Note that only functions that can be applied inplace (such as filtering/normalizing) will be applied inplace.
If our pipeline contained other types of functions, they will not be applied inplace, and will instead be returned at the end of the Pipeline.

If we apply a Pipeline with functions that return additional outputs (such as Figures, DataFrames, etc), they will be returned in a dictionary alongside the Filter object::

    >>> from *RNAlysis* import filtering
    >>> # create the pipeline
    >>> pipe = filtering.Pipeline('DESeqFilter')
    >>> pipe.add_function('biotypes_from_ref_table', ref='tests/test_files/test_files/biotype_ref_table_for_tests.csv')
    Added function 'DESeqFilter.biotypes_from_ref_table(ref='tests/test_files/test_files/biotype_ref_table_for_tests.csv')' to the pipeline.
    >>> pipe.add_function('filter_biotype_from_ref_table', 'protein_coding', ref='tests/test_files/test_files/biotype_ref_table_for_tests.csv')
    Added function 'DESeqFilter.filter_biotype_from_ref_table('protein_coding', ref='tests/test_files/test_files/biotype_ref_table_for_tests.csv')' to the pipeline.
    >>> pipe.add_function('biotypes_from_ref_table', ref='tests/test_files/test_files/biotype_ref_table_for_tests.csv')
    Added function 'DESeqFilter.biotypes_from_ref_table(ref='tests/test_files/test_files/biotype_ref_table_for_tests.csv')' to the pipeline.
    >>> # load the Filter object
    >>> d = filtering.DESeqFilter('tests/test_files/test_files/test_deseq_with_nan.csv')
    >>> # apply the Pipeline not-inplace
    >>> d_filtered, output_dict = pipe.apply_to(d, inplace=False)
    Biotype Reference Table used: tests/test_files/test_files/biotype_ref_table_for_tests.csv
    Biotype Reference Table used: tests/test_files/test_files/biotype_ref_table_for_tests.csv
    Filtered 2 features, leaving 26 of the original 28 features. Filtering result saved to new object.
    Biotype Reference Table used: tests/test_files/test_files/biotype_ref_table_for_tests.csv
    >>> print(output_dict['biotypes_1'])
                    gene
    biotype
    protein_coding    26
    pseudogene         1
    unknown            1
    >>> print(output_dict['biotypes_2'])
                    gene
    biotype
    protein_coding    26
    >>> # apply the Pipeline inplace
    >>> output_dict_inplace = pipe.apply_to(d)
    Biotype Reference Table used: tests/test_files/test_files/biotype_ref_table_for_tests.csv
    Biotype Reference Table used: tests/test_files/test_files/biotype_ref_table_for_tests.csv
    Filtered 2 features, leaving 26 of the original 28 features. Filtered inplace.
    Biotype Reference Table used: tests/test_files/test_files/biotype_ref_table_for_tests.csv

When an output dictionary is returned, the keys in the dictionary will be the name of the function appended to the number of call made to this function in the Pipeline (in the example above, the first call to 'biotypes_from_ref_table' is under the key 'biotypes_1', and the second call to 'biotypes_from_ref_table' is under the key 'biotypes_2'); and the values in the dictionary will be the returned values from those functions.
We can apply the same Pipeline to as many Filter objects as we want, as long as the type of the Filter object matches the Pipeline's `filter_type`.

****************************
*RNAlysis* enrichment module
****************************
RNAlysis's enrichment module (rnalysis.enrichment) can be used to perform various enrichment analyses including Gene Ontology (GO) enrichment and enrichment for user-defined attributes. The module also includes basic set operations (union, intersection, difference, symmetric difference) between different sets of genomic features.


Working with FeatureSet objects
=========================================
The enrichment module is built around :term:`FeatureSet` objects. A Featureset is a container for a set of gene/genomic feature IDs, and the set's name (for example, 'genes that are upregulated under hyperosmotic conditions'). All further anslyses of the set of features is done through the :term:`FeatureSet` object.


Initialize an FeatureSet object
------------------------------------------
We will start by importing the enrichment module::

    >>> from *RNAlysis* import enrichment

A :term:`FeatureSet` object can now be initialized by one of three methods.
The first method is to specify an existing Filter object::

    >>> my_filter_obj = filtering.CountFilter('tests/test_files/counted.csv') # create a Filter object
    >>> my_set = enrichment.FeatureSet(my_filter_obj, 'a name for my set')

The second method is to directly specify a python set of genomic feature indices, or a python set generated from an existing :term:`Filter object` (see above for more information about :term:`Filter objects` and the filtering module) using the function 'index_set'::

    >>> my_python_set = {'WBGene00000001','WBGene0245200',' WBGene00402029'}
    >>> my_set = enrichment.FeatureSet(my_python_set, 'a name for my set')
    # alternatively, using 'index_set' on an existing Filter object:
    >>> my_other_set = enrichment.FeatureSet(my_filter_obj.index_set,' a name for my set')

The third method is not to specify a gene set at all::

    >>> en = enrichment.FeatureSet(set_name = 'a name for my set')

At this point, you will be prompted to enter a string of feature indices seperated by newline. They will be automatically paresd into a python set.

FeatureSet objects have two attributes: gene_set, a python set containing genomic feature indices; and set_name, a string that describes the feature set (optional).

GO Enrichment
---------------
Using the *enrichment* module, you can perform enrichment analysis for Gene Ontology terms (GO enrichment).
You can read more about Gene Ontology on the `Gene Ontology Consortium website <http://geneontology.org/docs/ontology-documentation/?>`_.

To perform GO Enrichment analysis, we will start by creating an FeatureSet object::

    >>> counts = filtering.CountFilter('path_to_my_file.csv')
    >>> en = enrichment.FeatureSet(counts.index_set, 'my set')

Define the correct *organism* and *gene ID type* for your dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Since GO annotations refer to specific gene products, which can differ between different species, *RNAlysis* needs to know which organism your dataset refers to.
The organism can be specified as either the organism's name, or the organism's *NCBI Taxon ID* (for example: 6239 for *Caenorhabditis elegans*).

It is recommended to manually determine your organism's *NCBI Taxon ID* to avoid mischaracterization of annotations.
However, if you are not sure, *RNAlysis* will attempt to automatically determine the correct `organism` by default, based on the gene IDs in your FeatureSet.


Furthermore, since different annotations use different gene ID types to annotate the same gene products (such as UniProtKB ID, Entrez Gene ID, or Wormbase WBGene), *RNAlysis* can translate gene IDs from one gene ID type to another.
In order to do that, you need to specify which gene ID type your dataset uses.

Define the background set
^^^^^^^^^^^^^^^^^^^^^^^^^^
In enrichment analysis, we test whether our set of genomic features is enriched/depleted for a certain *GO Term*, in comparison to a more generalized set of genomic features that we determined as 'background'.
This could be the set of all protein-coding genes, the set of all genomic features that show expression above a certain threshold, or any other set of background genes which you deem appropriate. Importantly, the background set must contain all of the genes in the enrichment set.

Enrichment analysis is usually performed on protein-coding genes. Therefore, by default, *RNAlysis* uses all of the protein-coding genes that have at least one GO Annotation as a background set.
If you don't want to use the default setting, there are two methods of defining the background set:

The first method is to specify a biotype (such as 'protein_coding', 'miRNA' or 'all') under the parameter 'biotype'::

    >>> en.go_enrichment(biotype='all')

In this example, instead of using all of the protein-coding genes that have GO Annotations as background, we use every genomic feature with GO Annotations as background.
When specifying a biotype, the Biotype Reference Table that you specified is used to determine the biotype of each genomic feature.

The second method of defining the background set is to define a specific set of genomic features to be used as background::

    >>> my_background_set = {'feature1','feature2','feature3'}
    >>> en.go_enrichment(background_genes=my_background_set)

In this example, our background set consists of *feature1*, *feature2* and *feature3*.

It is not possible to specify both a biotype and a specific background set.

If some of the features in the background set or the enrichment set do no appear in the Reference Table, they will be ignored when calculating enrichment.

Choose the statistical test (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Significance testing for GO enrichment analysis can be done using either the Hypergeometric Test, Fisher's Exact Test, or a randomization test.

The hypergeometric test is defined as: Given *M* genes in the background set, *n* genes in the test set, with *N* genes from the background set belonging to a specific attribute ('success') and *X* genes from the test set belonging to that attribute.
If we were to randomly draw *n* genes from the background set (without replacement), what is the probability of drawing *X* or more (in case of enrichment)/*X* or less (in case of depletion) genes belonging to the given attribute?

The Fisher's Exact test is similar in principle to the hypergeometric test, but is two-tailed by default, as opposed to the hypergeometric test which examines enrichment and depletion separately.

The randomization test is defined as: Given *M* genes in the background set, *n* genes in the test set, with *N* genes from the background set belonging to a specific attribute and *X* genes from the test set belonging to that attribute.
We performs the number of randomizations specified by the user (10,000 by default).
In each randomization we randomly draw a set of *n* genes from the background set (without replacement), and marks the randomization as a 'success' if the number of genes in the random set belonging to the attribute is >= *X* (in case of enrichment) or <= *X* (in case of depletion).
The p-values are calculated as *(number of sucesses + 1)/(number of repetitions + 1)*.
This is a positive-bias estimator of the exact p-value, which avoids exactly-zero p-values.
You can read more about the topic in the following publication: https://www.ncbi.nlm.nih.gov/pubmed/21044043

If you don't specify which statistical test you want to use, the Fisher's Exact Test will be used by default.

To choose the statistical test you want to use, utilize the `statistical_test` parameter, which accepts either 'fisher', 'hypergeometric', or 'randomization'.
If you choose to use a randomization test, you can specify the number of randomization repititions to run using the `randomization_reps` parameter, and set the random seed using the `random_seed` parameter.

Filter GO Terms by *GO aspects* (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Gene Ontology considers three discrete aspects by which gene products can be described:

1. Biological process - the general 'biological objective' the gene product contributes to (e.g. 'locomotion', 'cell-cell signaling by wnt')
2. Molecular function - the molecular process or activity carried out by the gene product (e.g. 'antioxidant activity', 'ATP-dependent protein folding chaperone')
3. Cellular component - the location of the gene product when it carries out its action (e.g. 'P granule', 'mitochondrion')

Every GO term is exclusively associated with one of these *GO aspects*.
If you are interested in testing enrichment only for GO terms associated with a subset of these *GO aspects* you can specify which *GO aspects* to use through the `aspects` parameter.

If you don't specify *GO aspects* to be included, *RNAlysis* will test enrichment for GO Terms from all *GO aspects* by default.

Filter GO Annotations by Evidence Codes (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Every GO annotations is supported by an evidence code, which specifies what kind of evidence supports this annotation.
Evidence codes are divided into six categories:

1. experimental (there is evidence from an experiment directly supporting this annotation)
2. phylogenetic (annotations are derived from a phylogenetic model)
3. computational (annotations are based on in-silico analysis of gene sequence or other computational analysis)
4. author (annotations are based on the statement of the author in the cited reference)
5. curator (annotations are based on a curator's judgement)
6. electronic (annotations are based on homology and/or sequence information, and were not manually reviewed)

Each evidence category contains multiple evidence codes, each with its own definition.

You can choose to include only annotations with specific evidence codes, or to exclude annotations with specific annotation codes, using the `evidence_types` and `excluded_evidence_types` parameters.
You can specify either specific evidence codes (e.g. 'IEA', 'IKR'), evidence categories ('experimental', 'electronic'), or any combination of those.

If you don't specify evidence types to be included/excluded, *RNAlysis* will use annotations with all evidence codes by default.

You can read more about GO evidence codes here:
http://geneontology.org/docs/guide-go-evidence-codes/

Filter GO Annotations by database (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GO annotations are curated by different databases, such as UniProt, WormBase, or The Arabidopsis Information Resource.
You can choose to include only annotations from specific databases, or to exclude annotations from specific databases, using the `databases` and `excluded_databases` parameters.

If you don't specify databases to be included/excluded, *RNAlysis* will use annotations from all databases by default.

Filter GO Annotations by Qualifiers (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some GO annotations are modified by qualifiers. Each qualifier has a specific meaning within Gene Ontology:

1. the *NOT* qualifier - an explicit annotation that this particular gene product has been experimentally demonstrated *not* to be associated with the particular GO term.
Annotations with the *NOT* qualifier are usually ignored during enrichment analysis.
2. the *contributes_to* qualifier - indicates that this gene product facilitates, but does not directly carry out a function of a protein complex.
3. the *colocalizes_with* qualifier - indicates that this gene product associates with an organelle or complex.

You can choose to include only annotations with specific qualifiers, or to exclude annotations with a specific qualifier, using the `qualifiers` and `excluded_qualifiers` parameters.

If you don't specify qualifiers to be included/excluded, *RNAlysis* will ignore annotations with the *NOT* qualifier by default, and use annotations with any other qualifiers (or no qualifiers at all).

You can read more about GO qualifiers here:
http://geneontology.org/docs/go-annotations/

Choose annotation propagation method (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Gene Ontology terms have a somewhat hierarchical relationship that is defined as a directed a-cyclic graph (DAG). This means that each GO term is a node in the graph, and that each node has defined parents that are less specific than itself, going up to the top of the graph.

For example:

        .. figure:: /figures/http://geneontology.org/assets/hexose-biosynthetic-process.png
           :align:   center
           :scale: 35 %

           taken from the Gene Ontology Consortium site

We can see that in this example the GO term 'hexose biosynthetic process' has two parents, one of which is the less specific term 'hexose metabolic process', and these relationships go all the way up to the GO term 'metabolic process'.

Due to the relationships defined between GO terms, when a gene is annotated with a specific GO term, it makes sense that all of the less-specific parents of this GO term will also apply to said gene.
Therefore, when performing GO enrichment, we would usually 'propagate' every GO annotation to all of the GO terms upstream to it, all the way to the top of the GO graph.

Unfortunately, propagating GO annotations comes with some issues:
the defined relationship between GO terms introduces dependencies between neighboring GO terms, leading to correlation between enrichment results of different GO terms, and under-estimation of the False Discovery Rate of our analysis.
Moreover, since more specific GO terms by definition have less annotations than their less-specific parents, the most statistically significant enrichment results usually belong to the least-specific GO terms, which are not very biologically relevant.

To deal with this problem, several alternative propagation methods were developed to help de-correlate the GO graph structure and increase the specificity of our results without compromising accuracy.
You can read more about some suggested methods in the following publication:
https://pubmed.ncbi.nlm.nih.gov/16606683/

*RNAlysis* implements three of these propagation methods: *elim*, *weight*, and *all.m*.
You can decide which propagation method to use by specifying the `propagation_method` parameter: 'no' for no propagation of GO annotations, 'classic' for classic propagation of GO annotations, and 'elim'/'weight'/'all.m' for propagation using the *elim*/*weight*/*all.m* propagation algorithms.

If you don't specify which propagation method to use in enrichment analysis, the *elim* method will be used by default.

Choose plotting parameters (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After *RNAlysis* is done calculating the results of your enrichment analysis, it will automatically plot a summary of the enrichment results.
*RNAlysis* plots the results as a bar plot, with the Y axis showing log2 fold enrichment, and asterisks indicating whether this enrichment is statistically significant after correcting for multiple comparisons.

You can determine the orientation of the bar plot (horizontal or vertical) using the `plot_horizontal` parameter:

        .. figure:: /figures/plot_enrichment_results_go.png
           :align:   center
           :scale: 40 %

           `plot_horizontal`=True


        .. figure:: /figures/plot_enrichment_results_go_vertical.png
           :align:   center
           :scale: 40 %

           `plot_horizontal`=False


If you want to further customize this plot, you can request *RNAlysis* to return a Matplotlib Figure object of the barplot, by using the `return_fig` parameter.

If you don't specify plotting parameters, *RNAlysis* will generate a horizontal bar plot by default, and will not return a Matplotlib Figure object of the bar plot.


In addition, *RNAlysis* can generate an ontology graph, depicting all of the statistically significant GO terms and their hierarchical relationships:

        .. figure:: /figures/ontology_graph.png
           :align:   center
           :scale: 40 %

           `plot_ontology_graph`=True

If you don't want to generate this graph, you can set the parameter `plot_ontology_graph` to False.

Moreover, you can determine the file format of the generated graph (.pdf, .png, .svg, etc'), by setting the `ontology_graph_format` parameter.

If you don't specify plotting parameters, *RNAlysis* will generate an ontology graph by default in a PDF format.

Enrichment analysis output
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Running enrichment analysis will calculate enrichment for each of the GO terms, and return a pandas DataFrame in the following format:

+-------------+------------------+--------------+-----+-------+----------------------+----------+----------+-------------+
|             |       name       |    samples   | obs |   exp | log2_fold_enrichment |   pval   |   padj   | significant |
+=============+==================+==============+=====+=======+======================+==========+==========+=============+
|  GO:0001556 | oocyte maturation|    1327      | 451 | 319.52| 0.49722119558        | 0.0000999| 0.0000999| True        |
+-------------+------------------+--------------+-----+-------+----------------------+----------+----------+-------------+
|  GO:0043186 |     P granule    |    1327      | 89  | 244.87| -1.46013879322       | 0.0000999| 0.0000999| True        |
+-------------+------------------+--------------+-----+-------+----------------------+----------+----------+-------------+

'samples' is the number of features that were used in the enrichment set. 'obs' is the observed number of features positive for the attribute in the enrichment set.
'exp' is the expected number of features positive for the attribute in the background set. 'log2_fold_enrichment' is log2 of the fold change 'obs'/'exp'.

KEGG Pathways enrichment
_________________________
Using the *enrichment* module, you can perform enrichment analysis for KEGG pathways.
You can read more about KEGG pathways on the `KEGG website <https://www.genome.jp/kegg/pathway.html>`_.


To perform KEGG Enrichment analysis, we will start by creating an FeatureSet object::

    >>> counts = filtering.CountFilter('path_to_my_file.csv')
    >>> en = enrichment.FeatureSet(counts.index_set, 'my set')

Define the correct *organism* and *gene ID type* for your dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Since KEGG annotations refer to specific gene products, which can differ between different species, *RNAlysis* needs to know which organism your dataset refers to.
The organism can be specified as either the organism's name, or the organism's *NCBI Taxon ID* (for example: 6239 for *Caenorhabditis elegans*).

It is recommended to manually determine your organism's *NCBI Taxon ID* to avoid mischaracterization of annotations.
However, if you are not sure, *RNAlysis* will attempt to automatically determine the correct `organism` by default, based on the gene IDs in your FeatureSet.


Furthermore, since different annotations use different gene ID types to annotate the same gene products (such as UniProtKB ID, Entrez Gene ID, or Wormbase WBGene), *RNAlysis* can translate gene IDs from one gene ID type to another.
In order to do that, you need to specify which gene ID type your dataset uses.

Define the background set
^^^^^^^^^^^^^^^^^^^^^^^^^^
In enrichment analysis, we test whether our set of genomic features is enriched/depleted for a certain KEGG pathway, in comparison to a more generalized set of genomic features that we determined as 'background'.
This could be the set of all protein-coding genes, the set of all genomic features that show expression above a certain threshold, or any other set of background genes which you deem appropriate. Importantly, the background set must contain all of the genes in the enrichment set.

Enrichment analysis is usually performed on protein-coding genes. Therefore, by default, *RNAlysis* uses all of the protein-coding genes that have at least one KEGG annotation as a background set.
If you don't want to use the default setting, there are two methods of defining the background set:

The first method is to specify a biotype (such as 'protein_coding', 'miRNA' or 'all') under the parameter 'biotype'::

    >>> en.kegg_enrichment(biotype='all')

In this example, instead of using all of the protein-coding genes that have GO Annotations as background, we use every genomic feature with GO Annotations as background.
When specifying a biotype, the Biotype Reference Table that you specified is used to determine the biotype of each genomic feature.

The second method of defining the background set is to define a specific set of genomic features to be used as background::

    >>> my_background_set = {'feature1','feature2','feature3'}
    >>> en.kegg_enrichment(background_genes=my_background_set)

In this example, our background set consists of *feature1*, *feature2* and *feature3*.

It is not possible to specify both a biotype and a specific background set.

If some of the features in the background set or the enrichment set do no appear in the Reference Table, they will be ignored when calculating enrichment.

Choose the statistical test (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Significance testing for KEGG enrichment analysis can be done using either the Hypergeometric Test, Fisher's Exact Test, or a randomization test.

The hypergeometric test is defined as: Given *M* genes in the background set, *n* genes in the test set, with *N* genes from the background set belonging to a specific attribute ('success') and *X* genes from the test set belonging to that attribute.
If we were to randomly draw *n* genes from the background set (without replacement), what is the probability of drawing *X* or more (in case of enrichment)/*X* or less (in case of depletion) genes belonging to the given attribute?

The Fisher's Exact test is similar in principle to the hypergeometric test, but is two-tailed by default, as opposed to the hypergeometric test which examines enrichment and depletion separately.

The randomization test is defined as: Given *M* genes in the background set, *n* genes in the test set, with *N* genes from the background set belonging to a specific attribute and *X* genes from the test set belonging to that attribute.
We performs the number of randomizations specified by the user (10,000 by default).
In each randomization we randomly draw a set of *n* genes from the background set (without replacement), and marks the randomization as a 'success' if the number of genes in the random set belonging to the attribute is >= *X* (in case of enrichment) or <= *X* (in case of depletion).
The p-values are calculated as *(number of sucesses + 1)/(number of repetitions + 1)*.
This is a positive-bias estimator of the exact p-value, which avoids exactly-zero p-values.
You can read more about the topic in the following publication: https://www.ncbi.nlm.nih.gov/pubmed/21044043

If you don't specify which statistical test you want to use, the Fisher's Exact Test will be used by default.

To choose the statistical test you want to use, utilize the `statistical_test` parameter, which accepts either 'fisher', 'hypergeometric', or 'randomization'.
If you choose to use a randomization test, you can specify the number of randomization repititions to run using the `randomization_reps` parameter, and set the random seed using the `random_seed` parameter.


Choose plotting parameters (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After *RNAlysis* is done calculating the results of your enrichment analysis, it will automatically plot a summary of the enrichment results.
*RNAlysis* plots the results as a bar plot, with the Y axis showing log2 fold enrichment, and asterisks indicating whether this enrichment is statistically significant after correcting for multiple comparisons.

You can determine the orientation of the bar plot (horizontal or vertical) using the `plot_horizontal` parameter:

        .. figure:: /figures/plot_enrichment_results.png
           :align:   center
           :scale: 40 %

           `plot_horizontal`=True


        .. figure:: /figures/plot_enrichment_results_vertical.png
           :align:   center
           :scale: 40 %

           `plot_horizontal`=False


If you want to further customize this plot, you can request *RNAlysis* to return a Matplotlib Figure object of the barplot, by using the `return_fig` parameter.

If you don't specify plotting parameters, *RNAlysis* will generate a horizontal bar plot by default, and will not return a Matplotlib Figure object of the bar plot.

Enrichment analysis output
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Running enrichment analysis will calculate enrichment for each of the KEGG pathways, and return a pandas DataFrame in the following format:

+-----------+-----------------------------------------------------------------+--------------+-----+-------+----------------------+----------+----------+-------------+
|   KEGG ID |                              name                               |    samples   | obs |   exp | log2_fold_enrichment |   pval   |   padj   | significant |
+===========+=================================================================+==============+=====+=======+======================+==========+==========+=============+
|  cel00010 | Glycolysis / Gluconeogenesis - Caenorhabditis elegans (nematode)|    1327      | 451 | 319.52| 0.49722119558        | 0.0000999| 0.0000999| True        |
+-----------+-----------------------------------------------------------------+--------------+-----+-------+----------------------+----------+----------+-------------+
|  cel00030 |  Pentose phosphate pathway - Caenorhabditis elegans (nematode)  |    1327      | 89  | 244.87| -1.46013879322       | 0.0000999| 0.0000999| True        |
+-----------+-----------------------------------------------------------------+--------------+-----+-------+----------------------+----------+----------+-------------+

'samples' is the number of features that were used in the enrichment set. 'obs' is the observed number of features positive for the attribute in the enrichment set.
'exp' is the expected number of features positive for the attribute in the background set. 'log2_fold_enrichment' is log2 of the fold change 'obs'/'exp'.


Enrichment analysis for user-defined attributes
--------------------------------------------------
Using the *enrichment* module, you can perform enrichment analysis for user-defined attributes (such as 'genes expressed in intestine', 'epigenetic genes', 'genes that have paralogs'). The enrichment analysis can be performed using either the hypergeometric test or a randomization test.

Enrichment analysis for user-defined attributes is performed using either FeatureSet.enrich_hypergeometric, FeatureSet.enrich_randomization or FeatureSet.enrich_randomization_parallel. We will start by creating an FeatureSet object::

    >>> counts = filtering.CountFilter('path_to_my_file.csv')
    >>> en = enrichment.FeatureSet(counts.index_set, 'my set')

Choose which user-defined attributes to calculate enrichment for
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Our attributes should be defined in a Reference Table `csv` file. You can read more about Reference Tables and their format in the section :ref:`reference-table-ref`.
Once we have a Reference Table, we can perform enrichment analysis for those attributes using the function FeatureSet.enrich_randomization.
If your Reference Tables are set to be the default Reference Tables (as explained in :ref:`reference-table-ref`) you do not need to specify them when calling enrich_randomization. Otherwise, you need to specify your Reference Tables' path.
The names of the attributes you want to calculate enrichment for can be specified as a list of names (for example, ['attribute1', 'attribute2']).

Define the background set
^^^^^^^^^^^^^^^^^^^^^^^^^^
In enrichment analysis, we test whether our set of genomic features is enriched/depleted for a certain attribute, in comparison to a more generalized set of genomic features that we determined as 'background'.
This could be the set of all protein-coding genes, the set of all genomic features that show expression above a certain threshold, or any other set of background genes which you deem appropriate. Importantly, the background set must contain all of the genes in the enrichment set.

Enrichment analysis is usually performed on protein-coding genes. Therefore, by default, *RNAlysis* uses all of the protein-coding genes that appear in the Attribute Reference Table as a background set.
If you don't want to use the default setting, there are two methods of defining the background set:

The first method is to specify a biotype (such as 'protein_coding', 'miRNA' or 'all') under the parameter 'biotype'::

    >>> result = en.enrich_randomization(['attribute1','attribute2'], biotype='all')

In this example, instead of using all of the protein-coding genes in the Attribute Reference Table as background, we use all of the genomic features in the Attribute Reference Table as background.
When specifying a biotype, the Biotype Reference Table that you specified is used to determine the biotype of each genomic feature.

The second method of defining the background set is to define a specific set of genomic features to be used as background::

    >>> my_background_set = {'feature1','feature2','feature3'}
    >>> result = en.enrich_hypergeometric(['attribute1','attribute2'], background_genes=my_background_set)

In this example, our background set consists of feature1, feature2 and feature3.

It is not possible to specify both a biotype and a specific background set.

If some of the features in the background set or the enrichment set do no appear in the Reference Table, they will be ignored when calculating enrichment.

Choose the statistical test (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Significance testing for enrichment analysis can be done using either the Hypergeometric Test, Fisher's Exact Test, or a randomization test.

The hypergeometric test is defined as: Given *M* genes in the background set, *n* genes in the test set, with *N* genes from the background set belonging to a specific attribute ('success') and *X* genes from the test set belonging to that attribute.
If we were to randomly draw *n* genes from the background set (without replacement), what is the probability of drawing *X* or more (in case of enrichment)/*X* or less (in case of depletion) genes belonging to the given attribute?

The Fisher's Exact test is similar in principle to the hypergeometric test, but is two-tailed by default, as opposed to the hypergeometric test which examines enrichment and depletion separately.

The randomization test is defined as: Given *M* genes in the background set, *n* genes in the test set, with *N* genes from the background set belonging to a specific attribute and *X* genes from the test set belonging to that attribute.
We performs the number of randomizations specified by the user (10,000 by default).
In each randomization we randomly draw a set of *n* genes from the background set (without replacement), and marks the randomization as a 'success' if the number of genes in the random set belonging to the attribute is >= *X* (in case of enrichment) or <= *X* (in case of depletion).
The p-values are calculated as *(number of sucesses + 1)/(number of repetitions + 1)*.
This is a positive-bias estimator of the exact p-value, which avoids exactly-zero p-values.
You can read more about the topic in the following publication: https://www.ncbi.nlm.nih.gov/pubmed/21044043

If you don't specify which statistical test you want to use, the Fisher's Exact Test will be used by default.

To choose the statistical test you want to use, utilize the `statistical_test` parameter, which accepts either 'fisher', 'hypergeometric', or 'randomization'.
If you choose to use a randomization test, you can specify the number of randomization repititions to run using the `randomization_reps` parameter, and set the random seed using the `random_seed` parameter.

Choose plotting parameters (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After performing enrichment analysis, *RNAlysis* will automatically plot a summary of your enrichment results as a bar plot of log-transformed enrichment scores.
You can determine the orientation of the bar plot (horizontal or vertical) using the `plot_horizontal` parameter:


        .. figure:: /figures/plot_enrichment_results.png
           :align:   center
           :scale: 40 %

           `plot_horizontal`=True


        .. figure:: /figures/plot_enrichment_results_vertical.png
           :align:   center
           :scale: 40 %

           `plot_horizontal`=False


If you want to further customize your plot, you can retreive the matplotlib Figure object of your plot using the `return_fig` parameter.
When it is set as 'True', *RNAlysis* will return the Figure object it generated in addition to the results table.

Enrichment analysis output
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Running enrichment analysis will calculate enrichment for each of the specified attributes, and return a pandas DataFrame in the following format:

+----------------+--------------+-----+-------+----------------------+----------+----------+-------------+
|                |    samples   | obs |   exp | log2_fold_enrichment |   pval   |   padj   | significant |
+================+==============+=====+=======+======================+==========+==========+=============+
|     attribute1 |    1327      | 451 | 319.52| 0.49722119558        | 0.0000999| 0.0000999| True        |
+----------------+--------------+-----+-------+----------------------+----------+----------+-------------+
|     attribute2 |    1327      | 89  | 244.87| -1.46013879322       | 0.0000999| 0.0000999| True        |
+----------------+--------------+-----+-------+----------------------+----------+----------+-------------+

'samples' is the number of features that were used in the enrichment set. 'obs' is the observed number of features positive for the attribute in the enrichment set.
'exp' is the expected number of features positive for the attribute in the background set. 'log2_fold_enrichment' is log2 of the fold change 'obs'/'exp'.

Performing enrichment analysis for non-categorical user-defined attributes
---------------------------------------------------------------------------
Instead of peforming enrichment analysis for categorical attributes ("genes which are expressed exclusively in neurons", "genes enriched in males", "epigenetic gene-products", etc), you can test whether your FeatureSet is enriched for a non-categorical attribute ("number of paralogs", "gene length", or any other numeric attribute) using the function `FeatureSet.non_categorical_enrichment()`.

Choose which user-defined non-categorical attributes to calculate enrichment for
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The attributes should be defined in an Attribute Reference Table `csv` file. You can read more about Reference Tables and their format in the section :ref:`reference-table-ref`.
Once we have an Attirubte Reference Table, we can perform enrichment analysis for those non-categorical attributes using the function `FeatureSet.non_categorical_enrichment`.
If your Reference Tables are set to be the default Reference Tables (as explained in :ref:`reference-table-ref`) you do not need to specify them when calling non_categorical_enrichment. Otherwise, you need to specify your Reference Tables' path.
The names of the attributes you want to calculate enrichment for can be specified as a list of names (for example, ['attribute1', 'attribute2']).

Note that the variables you use for non-categorical enrichment analysis must be non-categorical, and must be defined for every genomic feature in the background and test sets (meaning, no NaN values).

Define the background set
^^^^^^^^^^^^^^^^^^^^^^^^^^
In enrichment analysis, we test whether our set of genomic features is enriched/depleted for a certain attribute, in comparison to a more generalized set of genomic features that we determined as 'background'.
This could be the set of all protein-coding genes, the set of all genomic features that show expression above a certain threshold, or any other set of background genes which you deem appropriate. Importantly, the background set must contain all of the genes in the enrichment set.

Enrichment analysis is usually performed on protein-coding genes. Therefore, by default, *RNAlysis* uses all of the protein-coding genes that appear in the Attribute Reference Table as a background set.
If you don't want to use the default setting, there are two methods of defining the background set:

The first method is to specify a biotype (such as 'protein_coding', 'miRNA' or 'all') under the parameter 'biotype'::

    >>> result = en.non_categorical_enrichment(['attribute1','attribute2'], biotype='all')

In this example, instead of using all of the protein-coding genes in the Attribute Reference Table as background, we use all of the genomic features in the Attribute Reference Table as background.
When specifying a biotype, the Biotype Reference Table that you specified is used to determine the biotype of each genomic feature.

The second method of defining the background set is to define a specific set of genomic features to be used as background::

    >>> my_background_set = {'feature1','feature2','feature3'}
    >>> result = en.non_categorical_enrichment(['attribute1','attribute2'], background_genes=my_background_set)

In this example, our background set consists of feature1, feature2 and feature3.

It is not possible to specify both a biotype and a specific background set.

If some of the features in the background set or the enrichment set do no appear in the Reference Table, they will be ignored when calculating enrichment.

Choose the statistical test (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When calculating enrichment for a non-categorical attribute, you can use either a parametric or non-parametric statistical test.

If the `parametric_test` parameter is True, *RNAlysis* will calculate enrichment using the *one-sample Student's T-test*, comparing the value of the non-categorical attribute in your enrichment set to the mean value of the non-categorical attribute in the background set.

If the `parametric_test` parameter is False, *RNAlysis* will calculate enrichment using the non-parametric *one-sample Sign test*, comparing the value of the non-categorical attribute in your enrichment set to the median value of the non-categorical attribute in the background set.

The parametric T-test assumes that the values of your non-categorical attribute distribute normally. Therefore, if you are not sure your data meets this assumption, it is recommended to use the less-powerful non-parametric test.

If you don't specify the statistical test, *RNAlysis* will automatically use a non-parametric test.

Choose plotting parameters (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The results of non-categorical enrichment analysis will be plotted as a histogram, comparing the distribution of the non-categorical attribute within the enrichment set compared to the background set.
If you used a parametric test, the mean of each set will be marked on the histogram. Alternatively, if you used a non-parametric test, the median of each set will be marked on the graph.
The enrichment results of each attribute will be plotted in a different Figure.

To better fit the visualization of the results to your needs, you can specify two parameters:

`plot_log_scale` can accept either `True` or `False`, and will determine whether the X-axis of the histogram will be plotting on a logarithmic or linear scale.

`plot_style' can accept either 'overlap' or 'interleaved', and will draw the histogram in one of the following two styles:

        .. figure:: /figures/hist_overlap.png
           :align:   center
           :scale: 40 %

           `plot_style`='overlap'

        .. figure:: /figures/hist_interleaved.png
           :align:   center
           :scale: 40 %

           `plot_style`='interleaved'

If you don't specify plotting parameters, *RNAlysis* will plot the histogram on a logarithmic scale using the 'overlap' style by default.

If you want to further customize your plot, you can retreive the matplotlib Figure object of your plot using the `return_fig` parameter.
When it is set as 'True', *RNAlysis* will return the Figure object it generated in addition to the results table.

Non-Categorical Enrichment analysis output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Running non-categorical enrichment analysis will calculate enrichment for each of the specified attributes, and return a pandas DataFrame in the following format:

+----------------+--------------+-------+--------+----------+----------+-------------+
|                |    samples   |  obs  |  exp   |   pval   |   padj   | significant |
+================+==============+=======+========+==========+==========+=============+
|     attribute1 |    1327      | 451   | 319.52 | 0.0000999| 0.0000999| True        |
+----------------+--------------+-------+--------+----------+----------+-------------+
|     attribute2 |    1327      | 89    | 244.87 | 0.0000999| 0.0000999| True        |
+----------------+--------------+-------+--------+----------+----------+-------------+

'samples' is the number of features that were used in the enrichment set. 'obs' is the observed mean/median (depending if the statistical test was parametric or not) value of the non-categorical attribute in the enrichment set.
'exp' is the expected mean/median value of the non-categorical attribute in the background set.


Performing set operations and visualisation on multiple FeatureSet objects
-------------------------------------------------------------------------------

Similarly to Filter objects, it is possible to use set operations such as union, intersection, difference and symmetric difference to combine the feature sets of multiple FeatureSet objects. Those set operations can be applied to both FeatureSet objects and python sets. The objects don't have to be of the same subtype - you can, for example, look at the union of an FeatureSet object and a python set::

    >>> en = enrichment.FeatureSet({'WBGene00003002','WBGene00004201','WBGene00300139'})

    >>> s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
    >>> union_result = en.union(s)

When performing set operations, the return type will always be a python set. This means you can use the output of the set operation as an input for yet another set operation, or as input to a new FeatureSet object.

In addition, the enrichment module includes functions for visualisation of sets and overlap, such as enrichment.venn_diagram() and enrichment.upset_plot().
Both functions receive a similar input: a dictionary whose keys are the names of the sets, and values are either FeatureSet objects, Filter objects, sets of genomic feature names, or the name of an attribute from the :term:`Attribute Reference Table` (you can read more about attributes in :ref:`reference-table-ref`).
Venn diagrams are limited to 2-3 sets:

       .. figure:: /figures/venn.png
           :align:   center
           :scale: 70 %

           Example plot of venn_diagram()

While UpSet plots can include any number of sets:

        .. figure:: /figures/upsetplot.png
           :align:   center
           :scale: 70 %

           Example plot of upset_plot()


Saving indices from FeatureSet to a .txt file
--------------------------------------------------------

It is possible to save the feature indices from an FeatureSet object to a .txt file, for use in online enrichment tools or simply to share the list of genomic features. This is done with the 'save_txt' function::

    >>> en.save_txt('D:\path\filename')

The feature indices will be saved to the text file in the specified path, separated by newline ('\n').

Working with RankedSet objects
=========================================
The :term:`RankedSet` class represents a set of genomic features which are ranked by an inherent order. For example, genes could be ranked based on their mean expression level, their fold change between a pair of conditions, etc.
:term:`RankedSet` objects behave similarly to :term:`FeatureSet` objects, but in addition they support *single-set enrichment analysis* - enrichment analysis without a background set.

The statistical test used to test significance in *single-set enrichment analysis* is the *generalized minimum hypergeometric test (XL-mHG)*, which measures the enrichment of the genomic features at the top and bottom of your ranked gene set compared to the rest of the gene set.
You can read more about the *minimum hypergeometric test* and its generalized version the *XL-mHG test* in the following publications:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030039
https://arxiv.org/abs/1507.07905

Initialize an RankedSet object
------------------------------------------
We will start by importing the enrichment module::

    >>> from *RNAlysis* import enrichment

:term:`RankedSet` objects can be initialized by one of two methods.
The first method is to specify an existing Filter object::

    >>> my_filter_obj = filtering.CountFilter('tests/test_files/counted.csv') # create a Filter object
    >>> my_ranked_set = enrichment.RankedSet(my_filter_obj, 'a name for my set')

When using this method, the ranking of the genomic features will be determined by their current order within the :term:`Filter` object. You can use the 'sort' method of your Filter object to modify the current order of genomic features within the Filter object.

The second method is to directly specify a list, tuple, or numpy array of genomic feature names::

    >>> my_list = ['WBGene00000001','WBGene0245200',' WBGene00402029']
    >>> my_ranked_set = enrichment.RankedSet(my_list, 'a name for my set')


RankedSet objects have three attributes: ranked_genes, a numpy array containing the ranked genomic features; gene_set, a python set containing the genomic features; and set_name, a string that describes the feature set (optional).


Performing single-set GO Enrichment analysis without a background set
-----------------------------------------------------------------------
Single-set GO Enrichment works much the same as normal Go Enrichment analysis, with two key differences:

First, when performing single-set GO Enrichment you do not define a background set to comapre your enrichment gene set against. Instead, you supply a :term:`RankedSet` that defines a meaningful ranking for your gene set.
Second, you cannot specify which statistical test to use, since the *XL-mHG* test has to be used.

An example for running single-set GO Enrichment would look like so::

    >>> from *RNAlysis* import enrichment
    >>> ranked_set = enrichment.RankedSet(['WBGene00000019', 'WBGene00000106', 'WBGene00000041', 'WBGene00000105'])
    >>> go_en_result = ranked_set.single_set_go_enrichment(gene_id_type='WormBase')


Performing single-set KEGG Enrichment analysis without a background set
-------------------------------------------------------------------------
Single-set KEGG Enrichment works much the same as normal KEGG Enrichment analysis, with two key differences:

First, when performing single-set KEGG Enrichment you do not define a background set to comapre your enrichment gene set against. Instead, you supply a :term:`RankedSet` that defines a meaningful ranking for your gene set.
Second, you cannot specify which statistical test to use, since the *XL-mHG* test has to be used.

An example for running single-set KEGG Enrichment would look like so::

    >>> from *RNAlysis* import enrichment
    >>> ranked_set = enrichment.RankedSet(['WBGene00000019', 'WBGene00000106', 'WBGene00000041', 'WBGene00000105'])
    >>> go_en_result = ranked_set.single_set_kegg_enrichment(gene_id_type='WormBase')

Performing single-set enrichment analysis for user-defined attributes without a background set
------------------------------------------------------------------------------------------------
Single-set enrichment analysis for user-defined attributes works much the same as normal enrichment analysis, with two key differences:

First, when performing single-set enrichment you do not define a background set to comapre your enrichment gene set against. Instead, you supply a :term:`RankedSet` that defines a meaningful ranking for your gene set.
Second, you cannot specify which statistical test to use, since the *XL-mHG* test has to be used.

An example for running single-set enrichment analysis would look like so::

    >>> from *RNAlysis* import enrichment
    >>> ranked_set = enrichment.RankedSet(['WBGene00000019', 'WBGene00000106', 'WBGene00000041', 'WBGene00000105'])
    >>> en_result = ranked_set.single_set_enrichment(['attribute1', 'attribute3'], attr_ref_path='tests/test_files/attr_ref_table_for_examples.csv')



Visualizing sets, intersections, and enrichment
================================================

Plotting results of enrichment analysis
-----------------------------------------
If you want to plot existing enrichment results using *RNAlysis*, you can use the `enrichment.plot_enrichment_results()` function. It employs a similar API to the enrichment functions in the enrichment modules, but accepts pre-calculated DataFrames.

Plotting Venn Diagrams and UpSet Plots
---------------------------------------

To visualize set intersection using a Venn Diagram or UpSet plot, you first need to define the sets you want to visualize.
Those sets can be represented using a python dictionary, where each key represents the name of the set, and the corresponding value represents the set's content (the genomic features that are part of the set).

There are three ways of specifying a set's content:

1. A python set containing the names of the genomic features
2. A `FeatureSet` or `RankedSet` object containing the names of the genomic features
3. A column name from an Attribute Reference Table

For example::

    >>> from *RNAlysis* import enrichment
    >>> sets_to_visualize = {'First set':{'gene1', 'gene2', 'gene3'}, 'Second set':enrichment.FeatureSet({'gene2','gene3','gene4'}), 'Third set':'attribute1'}

After defining the dictionary of sets, we can use it to plot set intersection via the `enrichment.venn_diagram()` or `enrichment.upset_plot()` functions.
Venn diagrams in *RNAlysis* are plotted with relatively accurate proportions and overlaps, and therefore only support up to 3 sets.

       .. figure:: /figures/venn.png
           :align:   center
           :scale: 70 %

           Example plot of venn_diagram()

UpSet plots support the visualization of much larger groups of sets. You can read more about UpSet plots here: https://upset.app/

        .. figure:: /figures/upsetplot.png
           :align:   center
           :scale: 70 %

           Example plot of upset_plot()

****************************
*RNAlysis* general module
****************************
RNAlysis's general module (rnalysis.general) contains general functions that can be useful during analysis of RNA sequencing data, including regular expression parsers and setting the Reference Table path.

.. _reference-table-ref:

Set and load a Reference Table
===============================

What is an Attribute Reference Table?
----------------------------------------
You can perform enrichment analysis or filtering operations based on user-defined attributes (such as 'genes expressed in intestine', 'epigenetic genes', 'genes that have paralogs').
User-defined attributes should be defined in an :term:`Attribute Reference Table` `csv` file. The format of the :term:`Attribute Reference Table` is one row for each gene/genomic feature, and one column for each attribute. Features that are negative for the attribute (for example, genes that have no paralogs under the attribute 'genes that have paralogs') should have the value NaN specified for the attribute, and features that are positive for the attribute (for example, genes that have paralogs under the attribute 'genes that have paralogs') should have any value other than NaN. The value could be either a boolean value (in our example, 'True' or '1' for genes that have paralogs), a number (in our example, the number of paralogs the gene has or the genomic distance to the nearest paralog), or any other value which is not NaN. See example for an :term:`Attribute Reference Table` below:

+----------------+--------------+-------------+-------------+
| feature_indices| attribute1   | attribute2  | attribute3  |
+================+==============+=============+=============+
| WBGene0000001  |      1       |     NaN     |     13.7    |
+----------------+--------------+-------------+-------------+
| WBGene0000002  |     NaN      |      1      |     241     |
+----------------+--------------+-------------+-------------+
| WBGene0000003  |     NaN      |      1      |     3.6     |
+----------------+--------------+-------------+-------------+
| WBGene0000004  |      1       |      1      |     NaN     |
+----------------+--------------+-------------+-------------+
| WBGene0000005  |      1       |     NaN     |     21.5    |
+----------------+--------------+-------------+-------------+

What is a Biotype Reference Table?
---------------------------------------
You can perform filtering operations or generate background-sets for enrichment analysis based on user-annotated biotypes (such as 'protein_coding', 'pseudogene', 'piRNA', etc).
User-annotated biotypes should be defined in a :term:`Biotype Reference Table` `csv` file. The format of the :term:`Biotype Reference Table` is one row for each gene/genomic feature, and a column titled 'biotype' (case insensitive). See example for a Biotype Reference Table below:

+----------------+----------------+
| feature_indices|    biotype     |
+================+================+
| WBGene0000001  | protein_coding |
+----------------+----------------+
| WBGene0000002  | protein_coding |
+----------------+----------------+
| WBGene0000003  |   pseudogene   |
+----------------+----------------+
| WBGene0000004  |     piRNA      |
+----------------+----------------+
| WBGene0000005  |    lincRNA     |
+----------------+----------------+



Set a Reference Table as default
----------------------------------
Once we have an Attribute and/or Biotype Reference Table, we can set it to be the default reference table for all future uses of *RNAlysis*::

    >>> from *RNAlysis* import general
    >>> path="the_new_attribute_reference_table_path"
    >>> general.set_attr_ref_table_path(path)
    Attribute Reference Table path set as: the_new_attribute_reference_table_path

    >>> path="the_new_biotype_reference_table_path"
    >>> general.set_biotype_ref_table_path(path)
    Attribute Reference Table path set as: the_new_biotype_reference_table_path

This will create a file called 'settings.yaml', which will store the full paths of your reference tables.
Whenever *RNAlysis* needs to use an Attribute/Biotype Reference Table and no other path is specified, *RNAlysis* will automatically use the path saved in the settings file.
The saved path can be changed any time using the general.set_attr_ref_table_path() and general.set_biotype_ref_table_path() functions.

Load the default Attribute Reference Table path
-------------------------------------------------
You can read the saved path from the settings file using the general.read_attr_ref_table_path() and general.read_biotype_ref_table_path() functions::

    >>> from *RNAlysis* import general
    >>> attr_table_path = general.read_attr_ref_table_path()
    Attribute Reference Table used: the_attribute_reference_table_path_that_was_saved_in_the_settings_file

    >>> biotype_table_path = general.read_biotype_ref_table_path()
    Biotype Reference Table used: the_biotype_reference_table_path_that_was_saved_in_the_settings_file

If an :term:`Attribute Reference Table` path was not previously defined, you will be requested to define it when you run this function.

Parse *C. elegans* gene names, WBGene indices and sequence names using regular expressions
===========================================================================================

The general module includes functions which can parse *C. elegans* gene names (like *daf-2* or *lin-15B*), WBGene indices (like WBGene00023495) and sequence names (like Y55D5A.5 or T23G5.6).
For example, we could extract all WBGene indices from the following string::

    >>> from *RNAlysis* import general
    >>> my_string='''WBGene00000001 and WBGene00000002WBGene00000003

            WBGene00000004g
            '''
    >>> general.parse_wbgene_string(my_string)
    {'WBGene00000001','WBGene000000002','WBGene00000003','WBGene00000004'}

