############################
RNAlysis user guide
############################


****************************
RNAlysis filtering module
****************************
RNAlysis's filtering module (rnalysis.filtering) is build to allow rapid and easy to understand filtering of various forms of RNA sequencing data. The module also contains specific methods for visualization and clustering of data. 

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
#link

In principle, any .csv file where the columns are different conditions/replicates and the rows include reads/normalized reads per genomic feature can be used as input for HTCountFilter. However, some HTCountFilter functions (such as 'norm_reads_to_rpm') will only work on HTSeq-count output files, and other unintended interactions may occur. 

Generating an HTCountFilter object from a folder of HTSeq-count output .txt files
---------------------------------------------------------------------------------
HTSeq-count receives as input an aligned SAM file. The native output of HTSeq-count is a text file with feature indices and read-per-genomic-feature, as well as information about reads that weren't counted for any feature (alignment not unique, low alignment quality, ambiguous, unaligned, aligned to no feature). When running HTSeq-count on multiple SAM files (which could represent different conditions or replicates), the final output would be a directory of .txt files. RNAlysis can parse those .txt files into two .csv tables: in the first each row is a genomic feature and each column is a condition or replicate (a single .txt file), and in the second each row represents a category of reads not mapped to genomic features (alignment not unique, low alignment quality, etc). This is done with the 'from_folder' function::

	h = filtering.HTCountFilter.from_folder('my_folder_path', save_reads_fname='name_for_reads_csv_file', save_not_counted_fname='name_for_unmapped_reads_csv_file')

By deault, 'from_folder' saves the generated tables as .csv files. However, you can avoid that by specifying 'save_csv=False'. 
It is also possible to automatically normalize the reads in the new HTCountFilter object to reads per million (RPM) using the unmapped reads data by specifying 'norm_to_rpm=True'. 


Loading from a pre-made .csv file
----------------------------------
If you have previously generated a .csv file from HTSeq-count output files using RNAlysis, or have done so manually, you can directly load this .csv file into an HTCountFilter object as you would any other Filter object::

	h = filtering.HTCountFilter('my_csv_file.csv')


Filtering operations unique to HTCountFilter
--------------------------------------------
#filter low reads

Normalizing reads with HTCountFilter
------------------------------------
HTCountFilter offers two methods for normalizing reads: reads per million (RPM) and DESeq2's size factors. Data normalized in other methods (such as RPKM) can be used as input for HTCountFilter, but it cannot perform such normalization methods on its own. 
#normalize to rpm
#normalize with size factors

Data visualization and clustering analysis with HTCountFilter
-------------------------------------------------------------
#pairplot
#clustermap
#pca

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
