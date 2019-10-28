############################
RNAlysis user guide
############################


****************************
RNAlysis filtering module
****************************
RNAlysis's filtering module (rnalysis.filtering) is built to allow rapid and easy to understand filtering of various forms of RNA sequencing data. The module also contains specific methods for visualization and clustering of data.

The filtering module is built around Filter objects, which are containers for tabular sequencing data. You can use the different types of Filter objects to apply filtering operations to various types of tabular data. You will learn more about Filter objects in the next section.

Working with Filter objects
============================

All Filter objects (DESeqFilter, HTCountFilter, FoldChangeFilter) work on the same principles,
and share many of the same functions and features. Each of them also has specific filtering, analysis and visualisation functions. In this section we will look into the general usage of Filter objects.

Initialize a Filter object
--------------------------

We will start by importing the filtering module::

    from rnalysis import filtering

We can now, for example, create a DESeqFilter object from a DESeq2 csv output file (see more details about DESeqFilter in sections below).
::

    d = filtering.DESeqFilter('D:\myfolder\my_deseq2_output.csv')

View a Filter object
--------------------

In order to view a glimpse of the file we imported we can use the 'head' and 'tail' functions.
By default 'head' will show the first 5 rows of the file, and 'tail' will show the last 5 rows,
but you can specify a specific number of lines to show.
::

    d.head()
    d.tail(8)

will give the following output:
#insert output here

We can also see the total number of rows an columns by accessing the 'shape' attribute::

    d.shape

which will return the following output:
#insert output here
meaning there are {} rows and {} columns in the file.

Filtering operations
--------------------

Now we can start filtering the entries in the file according to parameters of our choosing.
Various filtering operations are applied directly to the Filter object. Those operations do not affect the original csv file, but its representation within the Filter object.
For example, we can the function 'filter_percentile' to remove all rows that are above the specified percentile (in our example, 75% percentile) in the specified column (in our example, 'log2FoldChange')::

    d.filter_percentile(0.75,'log2FoldChange')

If we now look at the shape of d, we will see that {} rows have been filtered out of the object, and we remain with {} rows.

By default, filtering operations on Filter objects are performed in-place, meaning the original object is modified. However, we can save the results into a new Filter object and leave the current object unaffected by passing the argument 'inplace=False' to any filtering function within RNAlysis. For example::

    d_filtered = d.filter_percentile(0.75,'log2FoldChange',inplace=False)

In this case, the object 'd' will remain unchanged, while 'd_filtered' will be a new Filter object which contains our filtered results. We can continue applying filters sequentially to the same Filter object, or using 'inplace=False' to create a new object at any point.

Another useful option is to perform an opposite filter. When we specify the parameter 'opposite=True' to any filtering function within RNAlysis, the filtering function will be performed in opposite. This means that all of the genomic features that were supposed to be filtered out are kept in the object, and the genomic features that were supposed to be kept in the object are filtered out.
For example, if we now wanted to remove the rows which are below the 25% percentile in the 'log2FoldChange' column, we will use the following code::

    d.filter_percentile(0.25,'log2FoldChange',opposite=True)

Calling this function without the 'opposite' parameter would have removed all values except the bottom 25% of the 'log2FoldChange' column. When specifying 'opposite', we instead throw out the bottom 25% of the 'log2FoldChange' column and keep the rest.

There are many different filtering functions within the filtering module. Some of them are subtype-specific (such as 'filter_low_reads' for HTCountFilter objects and 'filter_significant' for DESeqFilter objects), while others can be applied to any Filter object. You can read more about the different functions and their usage in the project's documentation.


Performing set operations on multiple Filter objects
----------------------------------------------------

In addition to using regular filters, it is also possible to use set operations such as union, intersection, difference and symmetric difference to combine the results of multiple Filter objects. Those set operations can be applied to any Filter object, as well as to python sets. The objects don't have to be of the same subtype - you can, for example, look at the union of a DESeqFilter object, an HTCountFilter object and a python set::

    d = filtering.DESeqFilter('deseqfile.csv')
    h = filtering.HTCountFilter('htseq_count_file.csv')
    s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
    union_result = d.union(h, s)

When performing set operations, the return type will always be either a python set (default) or a string. This means you can use the output of the set operation as an input for yet another set operation. However, since the returned object is a set you cannot use Filter object functions such as 'head' and 'save_csv' on it, or apply filters to it directly.


Saving Filter results
---------------------

At any point we can save the current result of our filtering to a new csv file, by using the 'save_csv' function::

    d.save_csv()

If no filename is specified, the file is given a name automatically based on the filtering operations performed on it, their order and their parameters.
We can view the current automatic filename by looking at the 'fname' attribute::

    d.fname

In this example, the automatic filename is::

    'D:/myfolder/my_deseq2_output_below0.75percentile_below0.25percentileopposite.csv'
Alternatively, you can specify a filename::

    d.save_csv('alt_filename')

Instead of directly saving the results to a file, you can also get them as a set or string of genomic feature indices::

    set_output = d.features_set()
    str_output = d.features_string()

Sets of genomic feature indices can be used later for enrichment analysis using the enrichment module (see below).


Using a reference table for filter operations
---------------------------------------------

#set reference path
#read reference path
#what does reference table look like
#what do reference table attributes look like
#biotype table


Filtering DESeq2 output files with filtering.DESeqFilter
=========================================================

DESeqFilter objects are built to easily filter the output of R's DESeq2 package. This package is meant to analyze differential expression of genomic features in sequencing data. You can read more about it here: {}
Like other filter objects, filtering operations on DESeqFilter are performed in-place by default,meaning the original object is modified.

In principle, any .csv file that contains differential expression analysis data with log2 fold change and adjusted p values can be used as input for DESeqFilter.
However, some DESeqFilter functions (such as 'filter_significant' and 'filter_abs_log2_fold_change') may only work on DESeq2 output files, and other unintended interactions may occur.

Loading from a .csv file
------------------------
Loading a file into a DESeqFilter works as explained above for any Filter object::

    d = filtering.DESeqFilter('my_file.csv')

Filtering operations unique to DESeqFilter
------------------------------------------

There are a few filtering operations unique to DESeqFilter. Those include 'filter_significant', which removes statistically-insignificant rows according to a specified threshold; 'filter_abs_log2_fold_change', removes rows whose absolute value log2 fold change is below the specified threshold; 'filter_fold_change_direction' which removes either up-regulated (positive log2 fold change) or down-regulated (negative log2 fold change) rows; and 'split_fold_change_direction' which returns a DESeqFilter object with only up-regulated features and a DESeqFilter object with only down-regulated features.

The unique DESeqFilter filter operations expect specific column names (the column names automatically generated by DESeq2), and will not work with other column names:
'log2FoldChange','pval','padj'.


Filtering HTSeq-count output files with filtering.HTCountFilter
===============================================================

You can read more about HTSeq-count here:
https://htseq.readthedocs.io/en/release_0.11.1/count.html

In principle, any .csv file where the columns are different conditions/replicates and the rows include reads/normalized reads per genomic feature can be used as input for HTCountFilter. However, some HTCountFilter functions (such as 'norm_reads_to_rpm') will only work on HTSeq-count output files, and other unintended interactions may occur.

Generating an HTCountFilter object from a folder of HTSeq-count output .txt files
---------------------------------------------------------------------------------
HTSeq-count receives as input an aligned SAM/BAM file. The native output of HTSeq-count is a text file with feature indices and read-per-genomic-feature, as well as information about reads that weren't counted for any feature (alignment not unique, low alignment quality, ambiguous, unaligned, aligned to no feature). When running HTSeq-count on multiple SAM files (which could represent different conditions or replicates), the final output would be a directory of .txt files. RNAlysis can parse those .txt files into two .csv tables: in the first each row is a genomic feature and each column is a condition or replicate (a single .txt file), and in the second each row represents a category of reads not mapped to genomic features (alignment not unique, low alignment quality, etc). This is done with the 'from_folder' function::

    h = filtering.HTCountFilter.from_folder('my_folder_path', save_reads_fname='name_for_reads_csv_file', save_not_counted_fname='name_for_unmapped_reads_csv_file')

By deault, 'from_folder' saves the generated tables as .csv files. However, you can avoid that by specifying 'save_csv=False'.
It is also possible to automatically normalize the reads in the new HTCountFilter object to reads per million (RPM) using the unmapped reads data by specifying 'norm_to_rpm=True'.


Loading from a pre-made .csv file
----------------------------------
If you have previously generated a .csv file from HTSeq-count output files using RNAlysis, or have done so manually, you can directly load this .csv file into an HTCountFilter object as you would any other Filter object::

    h = filtering.HTCountFilter('my_csv_file.csv')


Filtering operations unique to HTCountFilter
--------------------------------------------
There are a few filtering operations unique to HTCountFilter. Those include 'filter_low_reads', which removes rows that have less than n reads in all columns.

Normalizing reads with HTCountFilter
------------------------------------
HTCountFilter offers two methods for normalizing reads: reads per million (RPM) and DESeq2's size factors. Data normalized in other methods (such as RPKM) can be used as input for HTCountFilter, but it cannot perform such normalization methods on its own.
#normalize to rpm
#normalize with size factors

Data visualization and clustering analysis with HTCountFilter
-------------------------------------------------------------
HTCountFilter includes multiple methods for visualization and clustering of count data.


With HTCountFilter.pairplot, you can get a quick overview of the distribution of counts within each sample, and the correlation between different samples:

.. figure::  pairplot.png
           :align:   center
           :scale: 40 %

           Example output of HTCount.pairplot()


With HTCount.clustergram, you can cluster your samples according to specified distance and linkage metrics:

 .. figure::  clustergram.png
           :align:   center
           :scale: 40 %

           Example plot of HTCount.clustergram()

With HTCount.pca, you can perform a principal component analysis and look for strong patterns in your dataset:

 .. figure::  pca.png
           :align:   center
           :scale: 40 %

           Example plot of HTCount.pca()

Filtering fold-change data of features using filtering.FoldChangeFilter
=======================================================================

#fold change is calculated (1+reads)/(1+reads)
#0 and inf are not supported, and may lead to unintented results

Loading fold change data from a .csv file
-----------------------------------------



Generating fold change data from an existing HTCountFilter object
-----------------------------------------------------------------




****************************
RNAlysis enrichment module
****************************
RNAlysis's enrichment module (rnalysis.enrichment) can be used to perform various enrichment analyses including gene ontology (GO) enrichment and enrichment for user-defined attributes. The module also includes basic set operations (union, intersection, difference, symmetric difference) between different sets of genomic features.


Working with EnrichmentProcessing objects
=========================================

The enrichment module is built around EnrichmentProcessing objects, which are a container for a set of genomic features and their name (for example, 'genes that are upregulated under hyperosmotic conditions'). All further anslyses of the set of features is done through the EnrichmentProcessing object.


Initialize an EnrichmentProcessing object
------------------------------------------
We will start by importing the enrichment module::

    from rnalysis import enrichment

An EnrichmentProcessing object can now be initialized by one of three methods.
The first method is to directly specify a python set of genomic feature indices::

    myset = {'WBGene00000001','WBGene0245200',' WBGene00402029'}
    en = enrichment.EnrichmentProcessing(myset, 'a name for my set')

The second method is to extract a python set of genomic feature indices from an existing Filter object (see above for more information about Filter objects and the filtering module) using the function 'features_set'::

    h = filtering.HTCountFilter('path_to_my_file.csv')
    en = enrichment.EnrichmentProcessing(h.features_set(), 'a name for my set')

The third method is not to specify a gene set at all::

    en = enrichment.EnrichmentProcessing(set_name = 'a name for my set')

At this point, you will be prompted to enter a string of feature indices seperated by newline. They will be automatically paresd into a python set.

EnrichmentProcessing objects have two attributes: gene_set, a python set containing genomic feature indices; and set_name, a string that describes the feature set (optional).


Randomization test enrichment analysis for user-defined attributes
-------------------------------------------------------------------
Using the enrichment module, you can perform enrichment analysis for user-defined attributes (such as 'genes expressed in intestine', 'epigenetic genes', 'genes that have paralogs'). The enrichment analysis is performed using a randomization test.

Enrichment analysis is performed using either EnrichmentProcessing.enrich_randomization or EnrichmentProcessing.enrich_randomization_parallel. We will start by creating an EnrichmentProcessing object::

    h = filtering.HTCountFilter('path_to_my_file.csv')
    en = enrichment.EnrichmentProcessing(h.features_set(), 'my set')

Our attributes should be defined in a Reference Table csv file. You can read more about Reference Tables and their format in the section :ref:`reference-table-ref`.
Once we have a Reference Table, we can perform enrichment analysis for those attributes using the function EnrichmentProcessing.enrich_randomization.
If our Reference Table is set to be the default Reference Table (as explained in :ref:`reference-table-ref`) we do not need to specify it when calling enrich_randomization. Otherwise, we need to specify our Reference Table's path.
The names of the attributes we want to calculate enrichment for can be specified as a list of names (for example, ['attribute1', 'attribute2']).

Next, we need to determine the set of genes to be used as background. Enrichment analysis is usually performed on protein-coding genes. Therefore, by default, enrich_randomization uses all of the protein-coding genes that appear in the Reference Table as a background set.
There are two methods of changing the default background set:

The first method is to specify a biotype (such as 'protein_coding', 'miRNA' or 'all') under the parameter 'biotype'::

    en.enrich_randomization(['attribute1','attribute2'], biotype='all')

In this example, instead of using all of the protein-coding genes in the Reference Table as background, we use all of the genomic features in the Reference Table as background.

The second method of changing the background set is to define a specific set of genomic features to be used as background::

    my_background_set = {'feature1','feature2','feature3'}
    en.enrich_randomization(['attribute1','attribute2'], background_genes=my_background_set)

In this example, our background set consists of feature1, feature2 and feature3.

If some of the features in the background set or the enrichment set do no appear in the Reference Table, they will be ignored when calculating enrichment.

Calling enrich_randomization will perform a randomization test for each of the specified attributes, and return a pandas DataFrame with the following format:

+----------------+--------------+-------+-------+----------------------+----------+----------+-------------+
|     name       |    samples   | n obs | n exp | log2_fold_enrichment |   pval   |   padj   | significant |
+================+==============+=======+=======+======================+==========+==========+=============+
|     attribute1 |    1327      | 451   | 319.52| 0.49722119558        | 0.0000999| 0.0000999| True        |
+----------------+--------------+-------+-------+----------------------+----------+----------+-------------+
|     attribute2 |    1327      | 89    | 244.87| -1.46013879322       | 0.0000999| 0.0000999| True        |
+----------------+--------------+-------+-------+----------------------+----------+----------+-------------+

'samples' is the number of features that were used in the enrichment set. 'n obs' is the observed number of features positive for the attribute in the enrichment set.
'n exp' is the expected number of features positive for the attribute in the enrichment set. 'log2_fold_enrichment' is log2 of the fold change 'n obs'/'n exp'.
enrich_randomization performs the number of randomizations specified by the user (10,000 by default), and marks each randomization as either a success or a failure.
The p values specified in 'pval' are calculated as (sucesses+1)/(repetitions+1). This is a positive-bias estimator of the exact p-value, which avoids exactly-zero p-values. You can read more about the topic in the following publication: https://www.ncbi.nlm.nih.gov/pubmed/21044043

If we want to perform the enrichment analysis in parallel and save time, we could use the enrich_randomization_parallel function instead of enrich_randomization.
To use it, you must first start a parallel session::

    from rnalysis import enrichment, general
    general.start_parallel_session()

To read more about parallel sessions, visit the :ref:`parallel-ref` section.
Afterwards, enrich_randomization_parallel is used exactly like enrich_randomization.

Performing set operations on multiple EnrichmentProcessing objects
-------------------------------------------------------------------

Similarly to Filter objects, it is possible to use set operations such as union, intersection, difference and symmetric difference to combine the feature sets of multiple EnrichmentProcessing objects. Those set operations can be applied to both EnrichmentProcessing objects and python sets. The objects don't have to be of the same subtype - you can, for example, look at the union of an EnrichmentProcessing object and a python set::

    en = enrichment.EnrichmentProcessing({'WBGene00003002','WBGene00004201','WBGene00300139'})

    s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
    union_result = en.union(s)

When performing set operations, the return type will always be a python set. This means you can use the output of the set operation as an input for yet another set operation, or as input to a new EnrichmentProcessing object.


Saving indices from EnrichmentProcessing to a .txt file
--------------------------------------------------------

It is possible to save the feature indices from an EnrichmentProcessing object to a .txt file, for use in online enrichment tools or simply to share the list of genomic features. This is done with the 'save_txt' function::

    en.save_txt('D:\path\filename')


The feature indices will be saved to the text file in the specified path, separated by newline ('\n').


****************************
RNAlysis general module
****************************
RNAlysis's general module (rnalysis.general) contains general functions
can be used to perform various enrichment analyses including gene ontology (GO) enrichment and enrichment for user-defined attributes. The module also includes basic set operations (union, intersection, difference, symmetric difference) between different sets of genomic features.

.. _parallel-ref:

Start and stop a parallel processing session
==============================================

Parallel processing in RNAlysis is performed using the ipyparallel package. You can read more about it here: https://ipyparallel.readthedocs.io/en/latest/
To use parallel processing features, you must first start an ipyparallel ipcluster. This is done using the general.start_parallel_session() function::

    from rnalysis import general
    general.start_parallel_session()

Your python console will then become unavailable for 30 seconds while the ipcluster is being started.
By default, the parallel session will use all available processors on the machine to perform parallel processing. You can specify the exact number of processors you want to use in the current session.

start_parallel_session() will automatically close the previous parallel session, start a new session, and sleep for 30 seconds while the ipcluster is being started. You can perform the same operations manually in order to skip the sleep period::

    from rnalysis import general
    general.start_ipcluster()
    #perform parallel processing here
    general.stop_ipcluster()


.. _reference-table-ref:

Set and load a Reference Table
===============================

What is a Reference Table?
----------------------------
You can perform enrichment analysis or filtering operations based on user-defined attributes (such as 'genes expressed in intestine', 'epigenetic genes', 'genes that have paralogs').
User-defined attributes should be defined in a Reference Table csv file. The format of the reference table is one row for each gene/genomic feature, and one column for each attribute. Features that are negative for the attribute (for example, genes that have no paralogs under the attribute 'genes that have paralogs') should have the value NaN specified for the attribute, and features that are positive for the attribute (for example, genes that have paralogs under the attribute 'genes that have paralogs') should have any value other than NaN. The value could be either a boolean value (in our example, 'True' or '1' for genes that have paralogs), a number (in our example, the number of paralogs the gene has or the genomic distance to the nearest paralog), or any other value which is not NaN. See example for a Reference Table below:

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


Set a Reference Table as default
---------------------------------
Once we have a Reference Table, we can set it to be the default Reference Table for all future uses of RNAlysis::

    from rnalysis import general
    general.set_reference_table_path('path/to/my/reference/table.csv')

This will create a file called 'settings.ini', which will store the full path of your Reference Table file.
Whenever RNAlysis needs to use a Reference Table and no other path is specified, RNAlysis will automatically use the path saved in the settings file.
The saved path can be changed any time using the general.set_reference_table_path() function.

Load the default Reference Table path
--------------------------------------
You can load the saved path from the settings file using the read_reference_table_path function::

    from rnalysis import general
    general.read_reference_table_path()

If a Reference Table path was not previously defined, you will be requested to define it when you run this function.

Parse *C. elegans* gene names, WBGene indices and sequence names using regular expressions
===========================================================================================








