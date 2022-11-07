####################################
A-to-Z tutorial - example analysis
####################################


***********************************
Intro and recommendations
***********************************

In this A-to-Z tutorial, we are going to analyze two datasets step-by-step, with each analysis covering different aspects of using *RNAlysis*.

These analyses both focus on RNA sequencing data from nematodes, but the same principles can be applied to data from any organism.

****************************************
Analysis #1 - time series RNA sequencing
****************************************

The dataset
=================
For the first analysis, we will use a publicly available dataset from `WormBase <https://downloads.wormbase.org/species/c_elegans/annotation/RNASeq_controls_FPKM/c_elegans.PRJNA13758.current.RNASeq_controls_FPKM.dat>`_. This dataset describes the mean expression level of each gene over the developmental stages of *C. elegans* nematodes, averaged from RNA-sequencing short-read experiments that have been identified in the literature as being 'controls'.

The original table linked above contains the mean expression, median expression, and number of samples for each developmental stage. Since we are only interested in the means, before starting the analysis, I copied the Mean columns to a new table, and gave the columns more readable names.
You can find the table I created `here <https://raw.githubusercontent.com/GuyTeichman/RNAlysis/master/tests/test_files/elegans_developmental_stages.tsv>`_.


Let's start by loading the dataset - in the main window, click on the "Load" button and choose the table's file from your computer.
We will then use the drop-down menu to change our table type from "Other" to "Count matrix". This will allow us to later on use analysis methods that are applicable only to count matrix-style datasets.
.. image:: /tutorial_screenshots/01a01_load_table.png
  :width: 600
  :alt: Loading a table into *RNAlysis* - choose the table type

Once we picked a table type, a new option will appear, allowing us to specify whether our table was pre-normalized. We will set this option to "True", since our dataset was already normalized to Fragments Per Kilobase Million (FPKM).
This is not a mandatory step, but if we don't do it, *RNAlysis* will warn us whenever we run analyses that expect normalized values.

.. image:: /tutorial_screenshots/01a02_load_table.png
  :width: 600
  :alt: Loading a table into *RNAlysis* - set the table as pre-normalized

Finally, we can click the "start" button to finish loading our dataset into *RNAlysis*.
The window will now display a preview of our table, as well as a a short summary of our table's content (table name, table type, number of rows and columns).

.. image:: /tutorial_screenshots/01b01_view_table.png
  :width: 600
  :alt: Preview the loaded table


If we want to see the entire table, we can click on the "View full table" button, and the entire table will appear in a new window:

.. image:: /tutorial_screenshots/01b02_view_table.png
  :width: 600
  :alt: View the table

At any point during the analysis, we can click the "Save table" button to save our table into a new file.

Data preprocessing and exploratory data analysis
=================================================

We can now begin exploring our data! Let's start with a pre-processing step - removing lowly-expressed genes.

Filter out lowly-expressed genes
---------------------------------

We want to filter out the genes that have not been expressed or that have low expression across all samples.
Lowly-expressed genes can negatively affect our analysis downstream, since the % error in the expression of these genes is relatively high, and these genes are likely to add noise rather than useful signal to our analysis.
We are going to filter our table, so that we keep only genes with 50 or more normalized reads in at least 1 experimental condition.

To apply a filtering function to our table, click on the "Filter" tab, and select a function from the drop-down menu that opens:

.. image:: /tutorial_screenshots/01c01_filter_low_reads.png
  :width: 600
  :alt: Choose a filtering function

We are going to select "Filter genes with low expression in all columns". This function will filter our any gene whose expression level is lower than X in every single column. This means we only keep genes that are reasonably expressed in at least one experimental condition.
If you are not sure what a function does, you can click on the blue question mark button next to the function's name to read a short description, or go to the function's help page by clicking on the blue link at the bottom of the main window.

Once we choose a function, the window will expand, and the filtering function's parameters will be displayed below. We can modify these parameters to choose exactly how to filter our table.

In this example, we are going to set the parameter `threshold` to 50 FPKM. This means that any genes which have less than 50 FPKM in **all** of the table's columns will be filtered out.
If you are not sure what a parameter does, you can hover your cursor over its name, or click on the blue question mark button next to the parameter's name.

Every filtering function in *RNAlysis* supports two additional parameters: `opposite` and `inplace`.

`opposite`, as the name indicates, allows you to apply the **inverse** of your filtering operation to your table. For example, in our case, instead of filtering out lowly-expressed genes, we will filter out everything **but** the lowly-expressed genes.

`inplace` allows us to choose whether we want to apply the filtering operation in the same window (default), or to keep an unfiltered version of the table and apply the filtering operation in a new tab.
This can be conveniet if we want to try multiple filtering approaches in parallel, or if we need to use the unfiltered table later down the line.
Regardless of what we choose, the original file we loaded will not be changed by our filtering unless we expelicitly save our filtering results, and any filtering operation we apply can be undone with the click of a button.

Once we are happy with the parameters we chose, we can scroll down and click on the "Apply" button to apply the filtering function to our table.

.. image:: /tutorial_screenshots/01c02_filter_low_reads.png
  :width: 600
  :alt: Set filtering parameters

After clicking "Apply", we can see that our action has been added to the *command history* pane on the right of the window. We can undo and redo any operation we apply inplace by clicking on a specific row in this pane.
Moreover, we can see that a short summary of our filtering action has printed to the log box at the bottom, and that the name and size of the table have been updated.

.. image:: /tutorial_screenshots/01c03_filter_low_reads.png
  :width: 600
  :alt: Filtering results

Examine variance in our data with Principal Component Analysis
---------------------------------------------------------------

The first analysis we will perform is Principal Component Analysis, or PCA.
Principal component analysis is a dimensionality-reduction method. Meaning, it can reduce the dimensionality of large data sets (e.g. the expression of thousands of genes), by transforming the expression data of these genes into a smaller dataset that still contains most of the information from the original dataset.

Since our dataset contains thousands of genes, it can be difficult to see how the different conditions differ in the expression of those genes. To remedy that, PCA analysis can "rearrange" the axes of our data such that we can summarize most of the variance of our dataset in very few dimensions (~2-10 dimensions).

You can read more about PCA `here <https://builtin.com/data-science/step-step-explanation-principal-component-analysis>`_.

To run a PCA, click on the "Visualize" tab, and select "Principal Component Analysis" from the drop-down menu.

.. image:: /tutorial_screenshots/01d01_pca.png
  :width: 600
  :alt: Choose a visualization function function

The 'samples' parameter allows you to choose which samples will be analyzed with PCA, and also lets you group these samples into sub-groups (for example, group replicate data by experimental condition), so that each sub-group will be drawn with a different color on the final graph.
In our case we only have one column per condition, and we want to examine them all, so we don't need to change this parameter.

By default, *RNAlysis* will apply a power-transform (Box-Cox) to the data before standardazing it and running PCA. This is the case for many functions in *RNAlysis*, since applying power transform minimizes undesirable characteristics of counts data, such as skeweness, mean-variance dependence, and extreme values.
However, this feature can always be disabled with the `power_transform` parameter.

Let's apply the analysis and look at the output:

.. image:: /tutorial_screenshots/01d02_pca.png
  :width: 600
  :alt: Principal Component Analysis

*RNAlysis* generated three graphs, depicting all of the pair-wise combinations between Principal Components #1-3.
We can visualize less or more principal components by changing the value of the `n_components` parameter.

Usually the most relevant graph is the one depicting the first and second PCs, since they explain the largest amount of variance in the data.

In our case, we can see that PCs 1-2 together explain ~75% of the variance in our dataset. Interestingly, the PCA shows a semi-circular pattern, ordered by the developmental stage of the worms.
My hypothesis would be that PC1 arranges the samples by their relative "germline" content - embryos are mostly gonads, adult nematodes contain a rather large quantity of germ cells, L4 larvae is the developmental stage where germline proliferation begins, and during the L1-L3 stages the relative germline content of the worms is relatively minimal.
PC2 appears to arrange the samples by their developmental stage, with embryos appearing at the top of the graph and adults at the bottom.

Examine similarity between developmental stages
-------------------------------------------------

Let's now examine the distribution of gene expression across developmental stages with a `Pair-plot<https://pythonbasics.org/seaborn-pairplot/>`_.
Pair-plots displays the pairwise relationships between samples, or experimental conditions, in our dataset, and also display a histogram of gene expression values within each sample/condition.

#TODO: stopped here


.. image:: /tutorial_screenshots/01e02_pairplot.png
  :width: 600
  :alt: Pair-plot

Compare the expression of specific genes over the developmental stages
-----------------------------------------------------------------------
We already have a hypothesis about the expression of some genes over the developmental stages.
Let's go to the "Visualize" tab to plot the expression of these genes over the developmental stages:

.. image:: /tutorial_screenshots/01f01_plot_expression.png
  :width: 600
  :alt: Choose genes to plot


.. image:: /tutorial_screenshots/01f02_plot_expression.png
  :width: 600
  :alt: Choose genes to plot - part 2


.. image:: /tutorial_screenshots/01f03_plot_expression.png
  :width: 600
  :alt: Expression plot

Clustering analysis
====================
We are interested in how different groups of genes change in expression level over the life cycle of the worm. We can use clustering analysis to group the genes in this dataset by their expression pattern over the developmental stages.
There is an abundance of approaches when it comes to clustering analysis of gene expression. To illustrate this point, we will cluster our data using three different types of clustering algorithms, arranged from the simplest to the most complex.

The simple approach - distance-based clustering with K-Medoids
--------------------------------------------------------------

K-Medoids clustering is a distance-based clustering method, where the algorithm attempts to divide all of our genes into K clusters (K is the number of clusters we are looking for), with each cluster being centered around one gene (Medoid).
K-Medoids clustering takes in two major parameters - the number of clusters we expect (K), and the distance metric by which we want to measure the similarity of expression between genes.

Specifying the exact number of clusters we expect is a bit challenging for us, since we aren't actually sure how many biologically-meaningful clusters are there in our data.
Moreover, this number could depend on how fine-grained we want our analysis to be - we could reasonably divide our genes into a small number of more generalized clusters (such as "genes expressed more in the start of development" vs "genes expressed more near the end of development"), or we could further divide our genes into smaller groups based on their exact temporal expression pattern.

Fortunately, some methods were developed to suggest a good number of clusters for our dataset (a "good number of clusters" meaning that the genes in each clusters are most similar to each other and most different than genes in other clusters). Two of these methods, known as the Silhouette Method and the Gap Statistic, are available in *RNAlysis*.
We will use the Gap Statistic method to determine some good options for the number of clusters in our analysis.

To start, let's click on the Clustering tab and choose K-Medoids clustering from the drop-down menu.
We can then set the value of the parameter `n_clusters` to 'gap', to indicate we want to use the Gap Statistic to determine the number of clusters in this analysis:

.. image:: /tutorial_screenshots/01g01_kmedoids.png
  :width: 600
  :alt: K-Medoids clustering setup - choose the number of clusters using the Gap Statistic

Next, we can set the distance metric. Different distance metrics can be more or less effective on specific types of data.
We will use a distance metric called YR1, that was developed especially for time-series gene expression data. You can read more about it in `Son and Baek 2007 <https://doi.org/10.1016/j.patrec.2007.09.015>`_:

.. image:: /tutorial_screenshots/01g02_kmedoids.png
  :width: 600
  :alt: K-Medoids clustering setup - choose the distance metric

We can now scroll all the way down, click the "Apply" button, and wait for the analysis to finish:

.. image:: /tutorial_screenshots/01g03_kmedoids.png
  :width: 600
  :alt: K-Medoids clustering - loading screen

One the clustering anslysis is finished, a few figures will open up. Let's examine them one by one.
The first figure will show us the results of the Gap Statistic algorithm. The graph on the left will show us, for each value of `n_clusters` tested, the natural logarithm (ln) of within-cluster dispersion.
We expect this value to go down the more clusters there are (the more clusters there are, the fewer genes will be in each cluster - therefore the genes within each cluster will be more similar to each other), and therefore we show both the actual dispersion values for the clustering solutions we calculated (the blue line), and also dispersion values for clustering solutions on random data, drawn from the same distribution (the orange line).
On the graph to the right we can see the natural logarithm of the ratio between the observed and expected dispersion - this is the Gap Statistic.
We are looking for local 'peaks' in the graph -  values of `n_clusters` that have a larger Gap Statistic than their neighbors.
*RNAlysis* automatically picks the lowest value of `n_clusters` that fits this criterion, but also suggests other good values based on the results.

.. image:: /tutorial_screenshots/01g04_kmedoids.png
  :width: 600
  :alt: Gap Statistic results

In our case, *RNAlysis* recommended 3 clusters as the optimal number of clusters. This value might not be granular enough for the kind of analysis we want to run. Therefore, we will re-run the K-Medoids algorithm with the same parameters, but set the value of `n_clusters` to one of the next good values discovered in the Gap Statistic analysis - 11 clusters.


Next, we can look at the rest of the output of the K-Medoids clustering algorithm for `n_clusters`=11.
The first graph will show us the distributions of gene expression in each discovered cluster. Note that the expression values are regularized and power-transformed, since we are interested in grouping the different genes by their relative pattern of expression, and not by their absolute expression levels (highly/lowly-expressed genes).
The clusters are sorted by their size, from the biggest to the smallest cluster.
This type of graph can help us see the general expression pattern that characterizes each cluster. Moreover, it can help point out how internally similar our clusters are - indicating the quality of our clustering result.

.. image:: /tutorial_screenshots/01g05_kmedoids.png
  :width: 600
  :alt: K-Medoids clustering results

In this case, we can see that while some clusters seem very internally consistent, quite a few clusters seem to contain a significant number of 'outlier' genes.

*RNAlysis* also generates a Principal Component Analysis graph of our gene expression data, marking the genes in each cluster with a different color.
This is another useful way to look at our clustering results - we would hope to see that the first two principal components explain a large degree of the variance in gene expression, and the genes in the same clusters will be grouped together in the graph.

.. image:: /tutorial_screenshots/01g06_kmedoids.png
  :width: 600
  :alt: K-Medoids clustering results - principal component analysis

Finally, the following window will open, prompting us to choose which output clusters we want to keep, and giving us the option to give these clusters a new name:

.. image:: /tutorial_screenshots/01g07_kmedoids.png
  :width: 600
  :alt: K-Medoids clustering results - choose which clusters to keep

For now, we will choose to keep none of the clusters, so that we can try out other clustering approaches. Therefore, we click either OK or Cancel without selecting any clusters.

Fine-tuning our approach - density-based clustering with HDBSCAN
-----------------------------------------------------------------

The next clustering approach we will use, HDBSCAN, belongs to a different category of clustering algorithms - density-based clustering.
HDBSCAN stands for Hierarchical Density-Based Spatial Clustering of Applications with Noise (see ` the publication <https://link.springer.com/chapter/10.1007/978-3-642-37456-2_14>`_ for more details).
HDBSCAN offers multiple advantages over more traditional clustering methods:

1. HSBSCAN makes relatively few assumptions about the data - it assumes that the data contains noise, as well as some real clusters which we hope to discover.
2. Unlike most other clustering methods, HDBSCAN does not "force" every gene to belong to a cluster. Instead, it can classify genes as outliers, excluding them from the final clustering solution.
3. HDBSCAN does not require you to guess the number of clusters in the data. The main tuning parameter in HDBSCAN is *minimum cluster size* (`min_cluster_size`), which determines the smallest "interesting" cluster size we expect to find in the data.

To run HDBSCAN, we need to pick a value for `min_cluster_size`.
Lower values of `min_cluster_size` will return a larger number of small clusters, revealing more fine-grained patterns in our gene expression data.
Higher values of `min_cluster_size` will return a smaller number of large clusters, revealing the most obvious and significant patterns in the data.
For our example, let's pick a value of 75:


.. image:: /tutorial_screenshots/01g11_hdbscan.png
  :width: 600
  :alt: HDBSCAN clustering setup

We will, once again, use YR1 as the distance metric.

If we look at the clustering results, we can see that HDBSCAN ended up generating a much larger number of clusters than the previous method, and they look fairly internally consistent.
.. image:: /tutorial_screenshots/01g12_hdbscan.png
  :width: 600
  :alt: HDBSCAN clustering results


.. image:: /tutorial_screenshots/01g13_hdbscan.png
  :width: 600
  :alt: HDBSCAN clustering results - principal component analysis

Let's now move on to the final clustering approach - ensemble-based clustering.

The complex approach - ensemble-based clustering with CLICOM
--------------------------------------------------------------

The last clustering approach we will use, CLICOM, is an emsemble-based clustering algorithm.
CLICOM (see ` the publication <https://doi.org/10.1016/j.eswa.2011.08.059>`_ ) incorporates the results of multiple clustering solutions, which can come from different clustering algorithms with differing clustering parameters, and uses these clustering solutions to create a combined "concensus" clustering solution.
CLICOM offers multiple advantages over more traditional clustering methods:

1. The ensemble clustering approach allows you to combine the results of multiple clustering algorithms with multiple tuning parameters, potentially making up for the weaknesses of each individual clustering method, and only taking into account patterns that robustly appear in many clustering solutions.
2. CLICOM does not require you to guess the final number of clusters in the data. The main tuning parameter in HDBSCAN is the *evidence threshold* (`evidence_threshold`).

*RNAlysis* offers a modified implementation of CLICOM. The modified version of the algorithm can, like the HDBSCAN algorithm, classify genes as outliers, excluding them from the final clustering solution.

This modified version of CLICOM supports a few tuning parameters, in addition to the clustering solutions themselves:

* `evidence_threshold`: how many clustering solutions (fraction between 0-1) have to agree about  two genes being clustered together in order for them to appear together in the final solution? A lower evidence threshold leads to fewer, large clusters, with fewer features being classified as outliers.
* `cluster_unclustered_features`: if True, CLICOM will force every gene to belong to a discovered cluster. Otherwise, genes can be classified as noise and remain unclustered.
* `min_cluster_size`: determines the minimal size of a cluster you would consider meaningful. Clusters smaller than this would be classified as noise and filtered out of the final result, or merged into other clusters (depending on the value of `cluster_unclustered_features`).


To start the analysis, let's choose the CLICOM algorithm from the Clustering drop-down menu. A new window will open:

.. image:: /tutorial_screenshots/01g21_clicom.png
  :width: 600
  :alt: CLICOM clustering setup

On the left half of the window we can set the tuning parameters of the CLICOM algorithm. For our example, let's set the evidence threshold to 0.5, and the minimum cluster size to 75.

On the right half of the window we can add new clustering setups to our run of CLICOM. These clustering setups can be any of the clustering algorithms available in *RNAlysis*, and we can add as many as we want - including multiple clustering setups using the same algorithm.
To add a setup, pick an algorithm from the drop-down menu, set it's parameters, and click the "Add Setup" button:

.. image:: /tutorial_screenshots/01g22_clicom.png
  :width: 600
  :alt: CLICOM clustering - add clustering setups

The setups you added will appear under "added setups" on the right, and you can delete a setup from that list if you want:

.. image:: /tutorial_screenshots/01g23_clicom.png
  :width: 600
  :alt: CLICOM clustering - examine added setups

Let's add the two clustering setups we used earlier, plus a few more:

.. image:: /tutorial_screenshots/01g24_clicom.png
  :width: 600
  :alt: CLICOM clustering - multiple clustering setups

Once we are happy with the clustering solutions and tuning parameters, we can click on the "Start CLICOM" button, and see progress reports in the output box on the main window of *RNAlysis*.

.. image:: /tutorial_screenshots/01g25_clicom.png
  :width: 600
  :alt: CLICOM clustering loading screen

Let's look at the final result:

.. image:: /tutorial_screenshots/01g26_clicom.png
  :width: 600
  :alt: CLICOM clustering results

.. image:: /tutorial_screenshots/01g27_clicom.png
  :width: 600
  :alt: CLICOM clustering results - principal component analysis

Let's choose a few good-looking clusters to keep, and give them a name that indicates their expression pattern:

.. image:: /tutorial_screenshots/01g28_clicom.png
  :width: 600
  :alt: The clusters we chose to keep

For this tutorial, I chose to keep clusters #1 ("down over development"), #2 ("L4 peak"), and #9 ("Down from L1 to adult").

Enrichment analysis
====================

Now that we have extracted a few clusters of interest, we can try to characterize their biological significance using enrichment analysis. We will look at the most commonly-used enrichment analysis - Gene Ontology enrichment.

Running enrichment analysis
----------------------------

To open the Enrichent Analysis window, open the 'Gene sets' menu and click "Enrichment Analysis":

.. image:: /tutorial_screenshots/01h01_go_enrichment.png
  :width: 600
  :alt: Pick 'Enrichment analysis' from the 'Gene sets' menu

For basic enrichment analysis, we first need to choose our *enrichment set* (the gene set we are interested in - for example, "Down from L1 to adult", cluster #9 we found earlier), and our background set (the reference genes we will be comparing our enrichment results to - for example, the genes in our original filtered table).
We can choose our sets from the two drop-down menus in the Enrichment window:

.. image:: /tutorial_screenshots/01h02_go_enrichment.png
  :width: 600
  :alt: Enrichment analysis setup - choose the enrichment and background sets

Next, we can choose the dataset we want to draw annotations from. In our case, we will choose Gene Ontology (GO).
After picking the dataset, the window expanded to show us all of the parameters we can modify for our analysis:

.. image:: /tutorial_screenshots/01h03_go_enrichment.png
  :width: 600
  :alt: Enrichment analysis setup - choose the analysis type, organism, and gene ID type

Let's set the statistical test to 'hypergeometric', the organism to 'Caenorhabditis elegans' (matching our gene expression data), and the Gene ID Type to "WormBase" (matching the gene ID type of our original data table).

We will leave the rest of the settings on the default values, but keep in mind that you can customize the analysis to a significant degree: using different statistical tests (including a statistical test that doesn't require a background gene set), using only specific types of GO annotations, propagating the annotations differently, etc'. You can read more about these options in the complete user guide.

Scroll to the bottom of thw window and click on the "run" button to run the analysis:

.. image:: /tutorial_screenshots/01h04_go_enrichment.png
  :width: 600
  :alt: Enrichment analysis loading screen

The enrichment window is going to minimize to allow you to read the log box on the main *RNAlysis* window, but you can enlarge it back if you want to look at the analysis parameters, or start a new analysis with the same parameters.

Once the analysis is done, we will be able to observe our results in multiple formats.

The first is a tabular format, showing all of the statistically significant GO terms we found (or the all of the tests GO terms, if we set the `return_nonsignificant` parameter to True).
The GO terms will be sorted by specificity, the most specific GO terms appearing at the top of the table.
The table also includes the statistics for each GO term (number of genes in the cluster, number of genes matching the GO term, expected number of genes to match the GO term, log2 fold change, p-value, and adjusted p-value).

.. image:: /tutorial_screenshots/01h06_go_enrichment.png
  :width: 600
  :alt: Enrichment analysis results - tabular format

The second output format is an Ontology Graph, depicting the statistically-significant GO terms in each GO aspect (biological process, cellular component, and molecular function), as well as their ancestors in the ontology graph.
The color of the terms on the graph indicates their log2 fold change, and the depth in the tree indicates the specificity of the term, with more specific GO terms being at the bottom.

.. image:: /tutorial_screenshots/01h05_go_enrichment.png
  :width: 800
  :alt: Enrichment analysis results - ontology graph

The final output format is a bar plot depicting the log2 fold change values, as well as significance, of the 10 most specific GO terms that were found to be statistically significant in our analysis.

.. image:: /tutorial_screenshots/01h07_go_enrichment.png
  :width: 600
  :alt: Enrichment analysis results - bar plot

*******************************************
Analysis #2 - differential expression data
*******************************************

The dataset
=================
We will start by analyzing XYZ. This dataset describes...

Let's start by loading the dataset:

#TODO: PtrSc

Data filtering and visualization with Pipelines
=================================================

#TODO

Apply the Pipeline to our datasets
-----------------------------------

Visualizing and extracting gene set interesctions
===================================================

Create a Venn diagram
-----------------------

Extract the subsets we are interested in
-----------------------------------------

Enrichment analysis
====================

Define a background set
-------------------------

Define our custom enrichment attributes
----------------------------------------

Running enrichment analysis
----------------------------
