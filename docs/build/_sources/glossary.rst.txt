.. _glossary:

Glossary
===================================

.. glossary::
   :sorted:

   Filter object
   Filter objects
      An object that stores your tabular data and its filename.
      You can apply filtering/normalizing functions to your data using the Filter object.
      Each type of Filter object supports different data structures and offers different utilities.
      For more details, see :term:`Filter`, :term:'CountFilter', :term:`DESeqFilter' and :term:`FoldChangeFilter'.

   Filter
      The most basic type of :term:`Filter object`.
      Can support any form of tabular data, and contains only the basic filtering functions.

   CountFilter
      A :term:`Filter object` that supports count matrices, where each row represents a gene or genomic feature,
      and every column represents a condition or replicate.

   DESeqFilter
      A :term:`Filter object` that supports differential expression data, formatted like the output files of DESeq2.
      See the user guide for more details.

   FoldChangeFilter
      A :term:`Filter object` that supports fold change values, where each row represents a gene or genomic feature,
      and the **only** other column contains fold-change values between two conditions.
      Can be loaded from a `csv` file, or generated from a :term:`CountFilter`.

   Pipeline
      An object that stores a series of functions and their corresponding arguments.
      You can then use the Pipeline object to apply this series of functions, with the same order and same arguments, to any number of :term:`Filter objects`.
      Pipelines can contain filtering functions, normalization functions, splitting functions and visualization functions.

   Attribute Reference Table
      A `csv` table that contains user-defined information ('attributes', such as 'genes expressed in intestine', 'epigenetic genes' or 'genes that have paralogs')
      abut genes/genomic features.
      Every row in the file is a gene/genomic feature, and every column is an attribute.
      You can define such table to be your default Attribute Reference Table using the function `general.set_attr_ref_table_path()`.
      The attributes in the table can be used both to filter your data with functions like `Filter.filter_by_attribute',
      and to perform enrichment analysis with various functions from the `enrichment` module.

   Biotype Reference Table
      A `csv` table that contains information about the biotype (protein coding, pseudogene, lncRNA, etc...) of genomic features.
      You can define such table to be your default Biotype Reference Table using the function `general.set_biotype_ref_table_path()`.
      Various functions in RNAlysis, such `Filter.biotypes_from_ref_table()` and `Filter.filter_biotype_from_ref_table()` will use the information in the Biotype Reference Table
      to filter data based on biotype, or display information about the biotype of the genomic features in your Filter objects.

   FeatureSet
      A container for a set of gene/feature IDs. A FeatureSet can optionally be named.
      Using the Enrichmemt module, you can run various enrichment analyses on FeatureSet objects.

   RankedSet
      A subtype of :term:`FeatureSet` that, instead of storing an unsorted set of gene/feature IDs,
      also stores their order. You can run enrichment analyses on RankedSet objects as you would on a FeatureSet object.
      However, you can also perform single-list enrichment (enrichment analysis without a background set, based on the ranking/order of the gene IDs)
      on RankedSet objects.
