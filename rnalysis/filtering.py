"""
This module can import HTCount and DESeq .csv files, perform various filtering operations on them,  \
perform set operations (union, intersection, etc), run basic exploratory analyses and plots (such as PCA, clustergram, \
violin plots, scatter, etc), save the filtered files to the computer, and return set of features that appear in an \
imported DESeq file. \
When a filtered/modified DESeq/HTCount is saved back to disk, its new name will include by default \
 all of the operations performed on it, in the order they were performed, to allow easy traceback of analyses.

"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
from rnalysis import general, __gene_names_and_biotype__
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt


class Filter:
    """
    A parent class for DESeqFilter and HTCountFilter.

    """

    def __init__(self, fname: str):
        if isinstance(fname, tuple):
            assert isinstance(fname[1], pd.DataFrame) and isinstance(fname[0], (str, Path))
            self.fname = fname[0]
            self.df = fname[1]
        else:
            assert isinstance(fname, (str, Path))
            self.fname = Path(fname)
            self.df = general.load_csv(fname, 0)
        if self.df.index.has_duplicates:
            warnings.warn("Warning: this Filter object contains multiple rows with the same WBGene index.")
        self.shape = self.df.shape
        self.columns = tuple(self.df.columns)

    def __str__(self):
        return f"{type(self).__name__} of file {self.fname}: \n{self.df.__str__()}"

    def __copy__(self):
        return type(self)((self.fname, self.df.copy(deep=True)))

    @staticmethod
    def _get_ref_path(ref):
        if ref == 'predefined':
            return general.read_reference_table_path()
        else:
            return ref

    def _inplace(self, new_df: pd.DataFrame, opposite: bool, inplace: bool, suffix: str):
        """
        Executes the user's choice whether to filter in-place or create a new instance of the Filter object.

        :param new_df: the post-filtering DataFrame
        :param opposite: boolean. Determines whether to return the filtration ,or its opposite.
        :param inplace: boolean. Determines whether to filter in-place or not.
        :param suffix: The suffix to be added to the filename
        :return:
        If inplace is False, returns a new instance of the Filter object.
        """
        assert isinstance(inplace, bool), "'inplace' must be True or False!"
        assert isinstance(opposite, bool), "'opposite' must be True or False!"
        if opposite:
            new_df = self.df.loc[self.df.index.difference(new_df.index)]
            suffix += 'opposite'

        new_fname = Path(f"{str(self.fname.parent)}\\{self.fname.stem}{suffix}{self.fname.suffix}")

        if inplace:
            self.df, self.fname = new_df, new_fname
            self.shape = self.df.shape
        else:
            tmp_df, tmp_fname = self.df, self.fname
            self.df, self.fname = new_df, new_fname
            new_obj = self.__copy__()
            self.df, self.fname = tmp_df, tmp_fname
            return new_obj

    def save_csv(self, alt_filename=None):
        """
        Saves the current filtered data to a .csv file.

        :param alt_filename: If None, file name will be generated automatically according to the filtering methods used. \
        if string, will be the name of the saved file. Example input: 'myfilename'
        """
        if alt_filename is None:
            alt_filename = self.fname
        else:
            alt_filename = f"{str(self.fname.parent)}\\{alt_filename}{self.fname.suffix}"
        general.save_to_csv(self.df, alt_filename, "")

    @staticmethod
    def _color_gen():
        """
        A generator that randomizes RGB values.

        :return:
        a numpy.ndarray of size (3,) containing three random values, each between 0 and 1.
        """
        while True:
            yield np.random.random(3)

    @staticmethod
    def _from_string(msg: str = '', delimiter: str = '\n'):
        """
        Takes a manual string input from the user, and then splits it using a delimiter into a list of values. \
        Called when an EnrichmentProcessing instance is created without input, \
        or when EnrichmentProcessing.enrich_big_table is called without input.

        :param msg: a promprt to be printed to the user
        :param delimiter: the delimiter used to separate the values. Default is '\n'
        :return:
        A list of the comma-seperated values the user inserted.
        """
        string = input(msg)
        split = string.split(sep=delimiter)
        if split[-1] == '':
            split = split[:-1]
        return split

    def head(self, n=5):
        """
        Return the first n rows of the DataFrame. See pandas.DataFrame.head documentation.

        :type n: int, default 5
        :param n: Number of rows to show.
        :return:
        returns the first n rows of the Filter object.
        """
        return self.df.head()

    def filter_percentile(self, percentile: float, column: str, opposite: bool = False, inplace: bool = True):
        """
        Removes all entries above the specified percentile in the specified column. \
        For example, if the column were 'pvalue' and the percentile was 0.5, then all features whose pvalue is above \
        the median pvalue will be filtered out.

        :type percentile: float between 0 and 1
        :param percentile: The percentile that all features above it will be filtered out.
        :type column: str
        :param column: Name of the DataFrame column according to which the filtering will be performed.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return:
        If inplace is False, returns a new and filtered instance of the Filter object.
        """
        assert isinstance(percentile, float), "percentile must be a float between 0 and 1!"
        assert isinstance(column, str) and column in self.df, "Invalid column name!"
        suffix = f'_below{percentile}percentile'
        new_df = self.df[self.df[column] < self.df[column].quantile(percentile)]
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_by_percentile(self, percentile: float, column: str):
        """
        Splits the Filter object into two Filter objects: \
        above and below the specified percentile in the spcfieid column.

        :type percentile: float between 0 and 1
        :param percentile: The percentile that all features above it will be filtered out.
        :type column: str
        :param column: Name of the DataFrame column according to which the filtering will be performed.
        :return:
        a list of two Filter objects: the first contains all of the features below the specified percentile, \
        and the second contains all of the features above and equal to the specified percentile.
        """
        return [self.filter_percentile(percentile=percentile, column=column, opposite=False, inplace=False),
                self.filter_percentile(percentile=percentile, column=column, opposite=True, inplace=False)]

    def filter_biotype(self, biotype='protein_coding',
                       ref: str = __gene_names_and_biotype__, opposite: bool = False, inplace: bool = True):
        """
        Filters out all features that do not match the indicated biotype. \
        Legal inputs: 'protein_coding','pseudogene','piRNA','miRNA','ncRNA','lincRNA','rRNA','snRNA','snoRNA'.

        :type biotype: str or list
        :param biotype: the biotypes which will not be filtered out.
        :param ref: Name of the reference file used to determine biotype. Default is the BigTable.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of Filter object.
        """
        legal_inputs = ('protein_coding', 'pseudogene', 'piRNA', 'miRNA', 'ncRNA', 'lincRNA', 'rRNA', 'snRNA', 'snoRNA')
        assert isinstance(biotype, (str, list)), "biotype must be a string or a list!"
        if isinstance(biotype, str):
            biotype = [biotype]
        for bio in biotype:
            assert bio in legal_inputs, f"biotype {bio} is not a legal string!"
        ref_df = general.load_csv(ref, 0)
        suffix = f"_{'_'.join(biotype)}"

        mask = pd.Series(np.zeros_like(ref_df['bioType'], dtype=bool), index=ref_df['bioType'].index, name='bioType')
        for bio in biotype:
            mask = mask | (ref_df['bioType'] == bio)

        gene_names = ref_df[mask].index.intersection(self.df.index)
        new_df = self.df.loc[gene_names]
        return self._inplace(new_df, opposite, inplace, suffix)

    # TODO: add 'remove unindexed rows' to here!

    def filter_by_bigtable_group(self, attributes: list = None, mode='union', inclusion: bool = True,
                                 ref: str = 'predefined',
                                 opposite: bool = False, inplace: bool = True):
        """
        Filters features by inclusion in or exclusion from a Big Table attribute, or multiple Big Table attributes. \
        When multiple attributes are given, filtering can be done in 'union' mode \
        (where features that belong to at least one attribute are not filtered out), or in 'intersection' mode \
        (where only features that belong to ALL attributes are not filtered out).

        :param attributes: list of Big Table attributes to filter by.
        :type mode: 'union' or 'intersection'.
        :param mode: 'union' or 'intersection'. If 'union', filters out every feature that does not match at least one \
        of the indicated Big Table attributes. If 'intersection', \
        filters out every feature that does not match all of the indicated Big Table attributes.
        :param inclusion: boolean. If True (default), features will be filtered by inclusion in the specified criteria \
        (meaning features that  belong to the attributes are kept \
        and features that don't belong to the attributes are filtered out). \
        If False, features will be filtered out by exclusion from the specified criteria \
        (meaning features that DON'T belong to the attributes \
        are kept, and features that do belong to the attributes are filtered out).
        :param ref: filename/path of the Big Table to be used as reference.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return:
        If 'inplace' is False, returns a new and filtered instance of the Filter object.
        """
        ref = self._get_ref_path(ref)
        if attributes is None:
            attributes = self._from_string(
                "Please insert attributes separated by newline "
                "(for example: \n'epigenetic_related_genes\nnrde-3 targets\nALG-3/4 class small RNAs')")
        elif isinstance(attributes, str):
            attributes = [attributes]
        else:
            assert isinstance(attributes, (list, tuple, set))
        assert isinstance(mode, str), "'mode' must be a string!"
        big_table = general.load_csv(ref, 0)
        sep_idx = [big_table[big_table[attr].notnull()].index for attr in attributes]

        if mode == 'intersection':
            suffix = '_bigtableintersection'
            indices = self.df.index
            for idx in sep_idx:
                indices = indices.intersection(idx)
        elif mode == 'union':
            suffix = '_bigtableUnion'
            indices = pd.Index([])
            for idx in sep_idx:
                indices = indices.union(idx)
            indices = indices.intersection(self.df.index)
        else:
            raise ValueError(f"Illegal input {mode}: mode must be either 'union' or 'intersection'")
        new_df = self.df.loc[set(indices)]
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_by_bigtable_group(self, attributes: tuple = None,
                                ref: str = 'predefined'):
        """
        Splits the Filter object into multiple Filter objects, \
        each corresponding to one of the specified Big Table attributes. \
        Each object contains only features that match its indicated Big Table attribute.

        :param attributes: list of Big Table attributes to filter by.
        :param ref: filename/path of the Big Table to be used as reference.
        :return:
        A list of Filter objects, each containing only features that match one Big Table attribute; the Filter objects \
        appear in the list in the same order the Big Table attributes were given in.
        """
        assert isinstance(attributes, (tuple, list, set))
        ref = self._get_ref_path(ref)
        return [self.filter_by_bigtable_group(attributes=[att], mode='union', ref=ref, inplace=False) for att in
                attributes]

    def describe(self, percentiles: list = [0.01, 0.25, 0.5, 0.75, 0.99]):
        """
        Generate descriptive statistics that summarize the central tendency, dispersion and shape \
        of the dataset’s distribution, excluding NaN values. \
        For more information see the documentation of pandas.DataFrame.describe.

        :type percentiles: list-like of numbers, optional
        :param percentiles: The percentiles to include in the output. \
        All should fall between 0 and 1. \
        The default is [.25, .5, .75], which returns the 25th, 50th, and 75th percentiles.
        :return:
        Summary statistics of the dataset.
        :rtype: Series or DataFrame
        """
        return self.df.describe(percentiles=percentiles)

    def features_set(self) -> set:
        """
        Returns all of the features in the current DataFrame (which were not removed by previously used filter methods) \
        as a set. \
        Warning: if any duplicate features exist in the filter object (same WBGene appears more than once), \
        the corresponding WBGene index will appear in the returned set ONLY ONCE.

        :return:
        A set of WBGene names.
        """
        if self.df.index.has_duplicates:
            warnings.warn("Warning: this filter object contains multiple rows with the same WBGene index. When "
                          "returning a set or string of features from this DESeqFilter object, each WBGene index will "
                          "appear ONLY ONCE!")
        return set(self.df.index)

    def features_string(self) -> str:
        """
        Returns a string of all WBGene indices in the current DataFrame separated by newline (the WBGene indices which \
        were not filtered out by previously-used filter methods). \
        Warning: if any duplicate features exist in the filter object (same WBGene appears more than once), \
        the corresponding WBGene index will appear in the returned string ONLY ONCE.

        :return:
        A string of WBGene indices separated by newlines (\n). \
        For example, "WBGene00000001\nWBGene00000003\nWBGene12345678".
        """
        return "\n".join(self.features_set())


class DESeqFilter(Filter):
    """
    A class that receives a DESeq output file and can filter it according to various characteristics.

    :param fname: name of the .csv DESeq file to be loaded.
    :type fname: str or pathlib.Path
    """

    def filter_significant(self, alpha: float = 0.1, opposite: bool = False, inplace: bool = True, ):
        """
        Removes all features which did not change significantly, according to the provided alpha.

        :param alpha: the significance threshold to determine which genes will be filtered. between 0 and 1.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current DESeqFilter object. If False, \
        the function will return a new DESeqFilter instance and the current instance will not be affected.
        :return:
        If 'inplace' is False, returns a new instance of DESeqFilter.
        """
        assert isinstance(alpha, float), "alpha must be a float!"
        new_df = self.df[self.df['padj'] <= alpha]
        suffix = f"_sig{alpha}"
        return self._inplace(new_df, opposite, inplace, suffix)

    def filter_top_n(self, n: int = 100, opposite: bool = False, inplace: bool = True, ):
        """
        Removes all features except the 'n' most significantly-changed features.

        :param n: an integer. How many genes to keep in the DESeqFilter instance
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current DESeqFilter object. If False, \
        the function will return a new DESeqFilter instance and the current instance will not be affected.
        :return:
        If 'inplace' is False, returns a new instance of DESeqFilter.
        """
        assert isinstance(n, int), "n must be an integer!"
        self.df.sort_values(by=['padj'])
        new_df = self.df.iloc[0:min(n, self.df.shape[0])]
        suffix = f"_top{n}"
        return self._inplace(new_df, opposite, inplace, suffix)

    def filter_abs_fold_change(self, abslog2fc: float = 1, opposite: bool = False, inplace: bool = True):
        """
        Filters out all features whose absolute log2 fold change is below the indicated threshold. \
        For example: if log2fc is 2.0, all features whose log2 fold change is between 1 and -1 (went up less than \
        two-fold or went down less than two-fold) will be filtered out.

        :param abslog2fc: The threshold absolute log2 fold change for filtering out a feature. Float or int. \
        All features whose absolute log2 fold change is lower than log2fc will be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current DESeqFilter object. If False, \
        the function will return a new DESeqFilter instance and the current instance will not be affected.
        :return:
        If 'inplace' is False, returns a new instance of DESeqFilter.
        """
        assert isinstance(abslog2fc, (float, int)), "abslog2fc must be a number!"
        assert abslog2fc >= 0, "abslog2fc must be non-negative!"
        suffix = f"_changed{abslog2fc}fold"
        new_df = self.df[np.abs(self.df['log2FoldChange']) >= abslog2fc]
        return self._inplace(new_df, opposite, inplace, suffix)

    def filter_fold_change_direction(self, direction: str = 'pos', opposite: bool = False, inplace: bool = True):
        """
        Filters out features according to the direction in which they changed between the two conditions.

        :param direction: 'pos' or 'neg'. If 'pos', will keep only features that have positive log2foldchange. \
        If 'neg', will keep only features that have negative log2foldchange.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current DESeqFilter object. If False, \
        the function will return a new DESeqFilter instance and the current instance will not be affected.
        :return:
        If 'inplace' is False, returns a new instance of DESeqFilter.
        """
        assert isinstance(direction, str), \
            "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. "
        if direction == 'pos':
            new_df = self.df[self.df['log2FoldChange'] > 0]
            suffix = '_PositiveLog2FC'
        elif direction == 'neg':
            new_df = self.df[self.df['log2FoldChange'] < 0]
            suffix = '_NegativeLog2FC'
        else:
            raise ValueError(
                "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. ")
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_fold_change_direction(self):
        """
        Splits the features in the current DESeqFilter object into two complementary, non-overlapping DESeqFilter \
        objects, based on the direction of their log2foldchange. The first object will contain only features with a \
        positive log2foldchange, the second object will contain only features with a negative log2foldchange.

        :return:
        a tuple containing two DESeqFilter objects: the first has only features with positive log2 fold change, \
        and the other has only features with negative log2 fold change.
        """
        return self.filter_fold_change_direction(direction='pos', inplace=False), self.filter_fold_change_direction(
            direction='neg', inplace=False)

    @staticmethod
    def __return_type(index_set: set, return_type: str):
        assert isinstance(return_type, str), "'return_type' must be a string!!"
        if return_type == 'set':
            return index_set
        elif return_type == 'str':
            return "\n".join(index_set)
        else:
            raise ValueError(f"'return type' must be either 'set' or 'str', is instead '{return_type}'!")

    def _set_ops(self, others, return_type, op):
        others = list(others)
        for i, other in enumerate(others):
            if isinstance(other, DESeqFilter):
                others[i] = other.features_set()
            elif isinstance(other, set):
                pass
            else:
                raise TypeError("'other' must be a DESeqFilter object or a set!")
        try:
            op_indices = op(set(self.df.index), *others)
        except TypeError as e:
            if op == set.symmetric_difference:
                raise TypeError(
                    f"Symmetric difference can only be calculated for two objects, {len(others) + 1} were given!")
            else:
                raise e
        return DESeqFilter.__return_type(op_indices, return_type)

    def intersection(self, *others, return_type: str = 'set'):
        """
        Returns a set/string of the WBGene indices that exist in ALL of the given DESeqFilter objects.

        :type others: DESeqFilter or set objects.
        :param others: Objects to calculate intersection with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the intersecting WBGene indices. If 'str', returns a string of \
        the intersecting WBGene indices, delimited by a comma.
        :rtype: set or str
        :return:
        a set/string of the WBGene indices that intersect between two DESeqFilter objects.
        """
        return self._set_ops(others, return_type, set.intersection)

    def union(self, *others, return_type: str = 'set'):
        """
        Returns a set/string of the union of WBGene indices between multiple DESeqFilter objects \
        (the indices that exist in at least one of the DESeqFilter objects).

        :type others: DESeqFilter or set objects.
        :param others: Objects to calculate union with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the union WBGene indices. If 'str', returns a string of \
        the union WBGene indices, delimited by a comma.
        :rtype: set or str
        :return:
         a set/string of the WBGene indices that exist in at least one of the DESeqFilter objects.
        """
        return self._set_ops(others, return_type, set.union)

    def difference(self, *others, return_type: str = 'set'):
        """
        Returns a set/string of the WBGene indices that exist in the first DESeqFilter object but NOT in the others.

        :type others: DESeqFilter or set objects.
        :param others: Objects to calculate difference with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the WBGene indices that exist only in the first DESeqFilter. \
        If 'str', returns a string of the WBGene indices that exist only in the first DESeqFilter, delimited by a comma.
        :rtype: set or str
        :return:
        a set/string of the WBGene indices that that exist only in the first DESeqFilter object (set difference).
        """
        return self._set_ops(others, return_type, set.difference)

    def symmetric_difference(self, other, return_type: str = 'set'):
        """
        Returns a set/string of the WBGene indices that exist either in the first DESeqFilter object OR the second, \
        but NOT in both (set symmetric difference).

        :type other: DESeqFilter or set.
        :param other: a second DESeqFilter object or set to calculate symmetric difference with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the WBGene indices that exist in exactly one DESeqFilter. \
        If 'str', returns a string of the WBGene indices that exist in exactly one DESeqFilter., delimited by a comma.
        :rtype: set or str
        :return:
        a set/string of the WBGene indices that that exist t in exactly one DESeqFilter. (set symmetric difference).
        """
        return self._set_ops([other], return_type, set.symmetric_difference)


class HTCountFilter(Filter):
    """
    A class that receives HTSeq count output files and can filter them according to various characteristics.

    :param fname: name of the .csv HTCount file to be loaded.
    :type fname: str or pathlib.Path
    """

    def pairplot(self, sample_list: list = 'all', log2: bool = False):
        """
        Plot pairwise relationships in the dataset. \
        Can plot both single samples and average multiple replicates. \
        For more information see the documentation of seaborn.pairplot.

        :type sample_list: 'all', list, or nested list.
        :param sample_list: A list of the sample names and/or grouped sample names to be included in the pairplot. \
        All specified samples must be present in the HTCountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :type log2: bool
        :param log2: if True, the pairplot will be calculated with log2 of the dataframe, and not with the raw data. \
        If False (default), the pairplot will be calculated with the raw data.
        :return:
        An instance of seaborn.PairGrid.

        .. figure::  pairplot.png
           :align:   center
           :scale: 40 %

           Example plot of pairplot()
        """
        if sample_list == 'all':
            sample_df = self.df
        else:
            sample_df = self._avg_subsamples(sample_list)
        if log2:
            pairplt = sns.pairplot(np.log2(sample_df))
        else:
            pairplt = sns.pairplot(sample_df)
        plt.show()
        return pairplt

    def _rpm_assertions(self, threshold: float = 1):
        """
        Various assertions for functions that normalize to RPM, or are meant to be used on pre-normalized values.

        :param threshold: optional. A threshold value for filter_low_rpm to be asserted.   
        """
        assert isinstance(threshold, (float, int)), "Threshold must be a number!"
        assert threshold >= 0, "Threshold must be zero or larger!"
        if 'rpm' not in str(self.fname) and 'sizefactor' not in str(self.fname):
            warnings.warn("Warning: using a function meant for normalized values on potentially unnormalized values!")

    def _avg_subsamples(self, sample_list: list):
        """
        Avarages subsamples/replicates according to the specified sample list. \
        Every member in the sample list should be either a name of a single sample (str), \
        or a list of multiple sample names to be averaged (list).

        :param sample_list: A list of the sample names and/or grouped sample names passed by the user. \
        All specified samples must be present in the HTCountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6'] \
        and the resulting output will be a DataFrame containing the following columns: \
        ['SAMPLE1', 'SAMPLE2', 'SAMPLE3', 'SAMPLE6']
        :return:
        a pandas DataFrame containing samples/averaged subsamples according to the specified sample_list.
        """
        samples_df = pd.DataFrame()
        for sample in sample_list:
            if isinstance(sample, str):
                samples_df[sample] = self.df[sample].values
            elif isinstance(sample, (list, str, tuple)):
                samples_df[",".join(sample)] = self.df[sample].mean(axis=1).values

        return samples_df

    def norm_reads_to_rpm(self, all_feature_fname: str, inplace: bool = True):
        """
        Normalizes the reads in the HTCountFilter to reads per million (RPM). \
        Uses a table of feature counts (ambiguous, no feature, not aligned, etc) from HTSeq's output. \
        Divides each column in the HTCountFilter object by (total reads + ambiguous + no feature)*10^-6 .

        :param all_feature_fname: the .csv file which contains feature information about the RNA library \
        (ambiguous, no feature, not aligned, etc).
        :param inplace: If True (default), filtering will be applied to the current HTCountFilter object. If False, \
        the function will return a new HTCountFilter instance and the current instance will not be affected.
        :return:
        If inplace is False, returns a new instance of the Filter object.
        """
        suffix = '_rpm'
        new_df = self.df.copy()
        if isinstance(all_feature_fname, (str, Path)):
            features = general.load_csv(all_feature_fname, 0)
        elif isinstance(all_feature_fname, pd.DataFrame):
            features = all_feature_fname
        else:
            raise TypeError("Invalid type for 'all_feature_fname'!")
        for column in new_df.columns:
            norm_factor = (new_df[column].sum() + features.loc[r'__ambiguous', column] + features.loc[
                r'__no_feature', column]) / (10 ** 6)
            new_df[column] /= norm_factor
        return self._inplace(new_df, opposite=False, inplace=inplace, suffix=suffix)

    def norm_reads_with_size_factor(self, size_factor_fname: str, inplace: bool = True):
        """
        Normalizes the reads in the HTCountFilter using pre-calculated size factors. \
        Such size factors can be calculated using DESeq2's median of ratios method. \
        Receives a table of sample names and their corresponding size factors, \
        and divides each column in the HTCountFilter object dataframe by the corresponding size factor.

        :type size_factor_fname: str or pathlib.Path
        :param size_factor_fname: the .csv file which contains size factors for the different libraries.
        :param inplace: If True (default), filtering will be applied to the current HTCountFilter object. If False, \
        the function will return a new HTCountFilter instance and the current instance will not be affected.
        :return:
        If inplace is False, returns a new instance of the Filter object.
        """
        suffix = '_sizefactor'
        new_df = self.df.copy()
        if isinstance(size_factor_fname, (str, Path)):
            size_factors = general.load_csv(size_factor_fname)
        elif isinstance(size_factor_fname, pd.DataFrame):
            size_factors = size_factor_fname
        else:
            raise TypeError("Invalid type for 'size_factor_fname'!")
        for column in new_df.columns:
            norm_factor = size_factors[column].values
            new_df[column] /= norm_factor
        return self._inplace(new_df, opposite=False, inplace=inplace, suffix=suffix)

    def filter_low_reads(self, threshold: float = 5, opposite: bool = False, inplace: bool = True):
        """
        remove all features which have less then 'threshold' reads per million in all conditions.

        :type threshold: float
        :param threshold: The minimal number of reads (counts, rpm, rpkm, tpm, etc) a feature should have \
        in at least one sample in order not to be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current HTCountFilter object. If False, \
        the function will return a new HTCountFilter instance and the current instance will not be affected.
        :return:
        If 'inplace' is False, returns a new instance of HTCountFilter.
        """
        self._rpm_assertions(threshold=threshold)
        new_df = self.df.loc[[True if max(vals) > threshold else False for gene, vals in self.df.iterrows()]]
        suffix = f"_filt{threshold}reads"
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_by_reads(self, threshold: float = 5):
        """
        Splits the features in the current HTCountFilter object into two complementary, non-overlapping HTCountFilter \
        objects, based on the their maximum expression level. The first object will contain only highly-expressed \
         features (which have reads over the specified threshold in at least one sample). The second object will \
         contain only lowly-expressed features (which have reads below the specified threshold in all samples).

        :param threshold: A float. The minimal number of reads (counts, RPM, RPKM, TPM etc) a feature needs to have \
        in at least one sample in order to be \
        included in the "highly expressed" object and no the "lowly expressed" object.
        :return:
        A tuple containing two HTCountFilter objects: the first has only highly-expressed features, \
        and the second has only lowly-expressed features.
        """
        self._rpm_assertions(threshold=threshold)
        high_expr = self.df.loc[[True if max(vals) > threshold else False for gene, vals in self.df.iterrows()]]
        low_expr = self.df.loc[[False if max(vals) > threshold else True for gene, vals in self.df.iterrows()]]
        return self._inplace(high_expr, opposite=False, inplace=False, suffix=f'_below{threshold}reads'), \
               self._inplace(low_expr, opposite=False, inplace=False, suffix=f'_above{threshold}reads')

    def clustergram(self, sample_names: list = 'all', metric: str = 'euclidean', linkage: str = 'average'):
        """
        Runs and plots a clustergram on the base-2 log of a given set of samples.

        :type sample_names: 'all' or list.
        :param sample_names: the names of the relevant samples in a list. \
        Example input: ["1A_N2_25", "1B_N2_25", "1C_N2_25", "2A_rde4_25", "2B_rde4_25", "2C_rde4_25"]
        :param metric: the distance metric to use in the clustergram. \
        Example inputs: 'euclidean', 'hamming', 'correlation'. \
        For all possible inputs see scipy.spatial.distance.pdist documentation online.
        :param linkage: the linkage method to use in the clustergram. \
        Example inputs: 'single', 'average', 'complete', 'ward'. \
        For all possible inputs see scipy.cluster.hierarchy.linkage documentation online.
        :return:
        A seaborn clustermap object.


        .. figure::  clustergram.png
           :align:   center
           :scale: 40 %

           Example plot of clustergram()
        """
        assert isinstance(metric, str) and isinstance(linkage, str), "Linkage and Metric must be strings!"
        metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean',
                   'hamming', 'jaccard', 'jensenshannon', 'kulsinski', 'mahalanobis', 'matching', 'minkowski',
                   'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
        linkages = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']
        assert metric in metrics and linkage in linkages

        if sample_names == 'all':
            sample_names = list(self.df.columns)
        print('Calculating clustergram...')
        plt.style.use('seaborn-whitegrid')
        clustering = sns.clustermap(np.log2(self.df[sample_names] + 1), method=linkage, metric=metric,
                                    cmap=sns.color_palette("RdBu_r", 10), yticklabels=False)
        plt.show()
        return clustering

    def pca(self, sample_names: list = 'all', n_components=3, sample_grouping: list = None):
        """
        runs and plots a PCA for a given set of samples.

        :type sample_names: 'all' or list.
        :param sample_names: the names of the relevant samples in a list. \
        Example input: ["1_REP_A", "1_REP_B", "1_REP_C", "2_REP_A", "2_REP_B", "2_REP_C", "2_REP_D", "3_REP_A"]
        :param n_components: number of PCA components to return \
        :param sample_grouping: Optional. Indicates which samples are grouped together as replicates, \
        so they will be colored similarly in the PCA plot. A list of indices from 0 and up, that indicates the sample \
         grouping. \
         For example, if sample_names is: \
        ["1_REP_A", "1_REP_B", "1_REP_C", "2_REP_A", "2_REP_B", "2_REP_C", "2_REP_D", "3_REP_A"], \
        then the sample_grouping will be: \
        [0, 0, 0, 1, 1, 1, 1, 2]
        :return:
        A tuple whose first element is an sklearn.decomposition.pca object, \
        and second element is a list of matplotlib.axis objects.

        .. figure::  pca.png
           :align:   center
           :scale: 40 %

           Example plot of pca()
        """
        if sample_names == 'all':
            sample_names = list(self.df.columns)
            srna_data = self.df.transpose()
        else:
            srna_data = self.df[sample_names].transpose()
        srna_data_norm = StandardScaler().fit_transform(srna_data)

        pca_obj = PCA(n_components=n_components)
        pcomps = pca_obj.fit_transform(srna_data_norm)
        columns = [f'Principal component {i + 1}' for i in range(n_components)]
        principal_df = pd.DataFrame(data=pcomps, columns=columns)
        final_df = principal_df
        final_df['lib'] = pd.Series(sample_names)

        pc_var = pca_obj.explained_variance_ratio_
        graphs = 2 if n_components > 2 else 1
        axes = []
        for graph in range(graphs):
            axes.append(HTCountFilter._plot_pca(
                final_df=final_df[['Principal component 1', f'Principal component {2 + graph}', 'lib']],
                pc1_var=pc_var[0], pc2_var=pc_var[1 + graph], sample_grouping=sample_grouping))

        return pca_obj, axes

    @staticmethod
    def _plot_pca(final_df: pd.DataFrame, pc1_var: float, pc2_var: float, sample_grouping: list):
        """
        Internal method, used to plot the results from HTCountFilter.pca. Static class method.

        :param final_df: The DataFrame output from pca
        :param pc1_var: Variance explained by the first PC.
        :param pc2_var: Variance explained by the second PC.
        :param sample_grouping: a list of indices from 0 and up, that indicates what samples are grouped together as \
        biological replicates. For example, if sample_names is: \
        ["1A_N2_25", "1B_N2_25", "1C_N2_25", "2A_rde4_25", "2B_rde4_25", "2C_rde4_25"], \
        then the sample_grouping will be: \
        [0,0,0,1,1,1]
        :return:
        an axis object containing the PCA plot.
        """
        plt.style.use('seaborn-whitegrid')

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.grid(True)
        ax.set_xlabel(f'{final_df.columns[0]} (explained {pc1_var * 100 :.2f}%)', fontsize=15)
        ax.set_ylabel(f'{final_df.columns[1]} (explained {pc2_var * 100 :.2f}%)', fontsize=15)
        ax.set_title('PCA', fontsize=20)

        color_generator = HTCountFilter._color_gen()
        if sample_grouping is None:
            colors = [next(color_generator) for _ in range(len(final_df.columns))]
        else:
            color_opts = [next(color_generator) for _ in range(max(sample_grouping))]
            colors = [color_opts[i - 1] for i in sample_grouping]

        ax.scatter(final_df.iloc[:, 0], final_df.iloc[:, 1], c=colors, s=50)
        for _, row in final_df.iterrows():
            row[0] += 1
            row[1] += 1
            ax.text(*row)
        ax.grid(True)
        return ax

    def scatter_sample_vs_sample(self, sample1: str, sample2: str, xlabel: str = None, ylabel: str = None,
                                 deseq_highlight=None):
        """
        Generate a scatter plot where every dot is a feature, the x value is log10 of reads \
        (counts, RPM, RPKM, TPM, etc) in sample1, the y value is log10 of reads in sample2.

        :param sample1: str/list. Name of the first sample from the HTCountFilter object. \
        If sample1 is a list, they will be avarged as replicates.
        :param sample2: str/list. Name of the second sample from the HTCountFilter object. \
        If sample2 is a list, they will be averaged as replicates.
        :param xlabel: optional. If not specified, sample1 will be used as xlabel.
        :param ylabel: optional. If not specified, sample2 will be used as ylabel.
        :param deseq_highlight: DESeqFilter object or iterable of WBGene indices(optional). \
        If specified, the points in the scatter corresponding to the WBGene indices in deseq_highlight will be \
        highlighted in red.
        :return:
        a matplotlib axis object.

        .. figure::  rpm_vs_rpm.png
           :align:   center
           :scale: 60 %

           Example plot of scatter_sample_vs_sample()

        """
        self._rpm_assertions()
        assert isinstance(sample1, (str, list, tuple, set)) and isinstance(sample2, (str, list, tuple, set))
        xvals = np.log10(self.df[sample1].values + 1) if isinstance(sample1, str) else np.log10(
            self.df[sample1].mean(axis=1).values + 1)
        yvals = np.log10(self.df[sample2].values + 1) if isinstance(sample2, str) else np.log10(
            self.df[sample2].mean(axis=1).values + 1)

        plt.style.use('seaborn-whitegrid')
        if xlabel is None:
            xlabel = f'log10(reads per million) from library {sample1}'
        if ylabel is None:
            ylabel = f'log10(reads per million) from library {sample2}'
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(xlabel, fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        ax.set_title(f'{sample1} vs {sample2}', fontsize=20)
        ax.scatter(xvals, yvals, s=3, c='#6d7178')

        if deseq_highlight is not None:
            highlight_features = deseq_highlight.features_set() if isinstance(deseq_highlight,
                                                                              DESeqFilter) else deseq_highlight
            xvals_highlight = np.log10(self.df[sample1].loc[highlight_features].values + 1) if \
                isinstance(sample1, str) else np.log10(self.df[sample1].loc[highlight_features].mean(axis=1).values + 1)
            yvals_highlight = np.log10(self.df[sample2].loc[highlight_features].values + 1) if \
                isinstance(sample2, str) else np.log10(self.df[sample2].loc[highlight_features].mean(axis=1).values + 1)

            ax.scatter(xvals_highlight, yvals_highlight, s=3, c=np.array([[0.75, 0.1, 0.1]]))
        plt.show()
        return ax

    def violin_plot(self, samples='all'):
        """
        Generates a violin plot of the specified samples in the HTCountFilter object. \
        Can plot both single samples and average multiple replicates. \
        It is recommended to use this function on normalized values and not on absolute read values. \
        Box inside the violin plot indicates 25% and 75% percentiles, and the white dot indicates the median.

        :type samples: 'all' or list.
        :param samples: A list of the sample names and/or grouped sample names to be plotted in the violin plot. \
        All specified samples must be present in the HTCountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :return:
        a seaborn violin object.

        .. figure::  violin.png
           :align:   center
           :scale: 60 %

           Example plot of violin_plot()

        """
        self._rpm_assertions()
        if samples == 'all':
            samples_df = self.df
        else:
            samples_df = self._avg_subsamples(samples)

        samples_df = np.log10(samples_df + 1)
        fig = plt.figure(figsize=(8, 8))

        violin = sns.violinplot(data=samples_df)
        plt.style.use('seaborn-whitegrid')
        plt.xlabel("Samples")
        plt.ylabel("Log10 RPM")
        plt.show()
        return violin

    # TODO: add ranksum test
    # TODO: fix futureWarning scipy

    @staticmethod
    def from_folder(folder_path: str, save_csv: bool = True, norm_to_rpm: bool = False, save_reads_fname: str = None,
                    save_not_counted_fname: str = None):
        """
        Iterates over HTSeq count .txt files in a given folder and combines them into a single HTCountFilter table. \
        Can also save the count data table and the uncounted data table to .csv files, and normalize the HTCountFilter \
        table to reads per million (RPM). Note that the saved data will always be count data, and not normalized data, \
        regardless if the HTCountFilter table was normalized or not.

        :param folder_path: str or pathlib.Path. Full path of the folder that contains individual htcount .txt files.
        :param save_csv: bool. If True (default), the joint DataFrame of count data and uncounted data will be saved \
        to two separate .csv files. The files will be saved in 'folder_path', and named according to the parameters \
        'save_reads_fname' for the count data, and 'save_not_counted_fname' for the uncounted data (unaligned, \
        alignment not unique, etc).
        :param norm_to_rpm: bool. If True, the HTCountFilter table will be automatically normalized to reads per \
        million (RPM). If False (defualt), the HTCountFilter object will not be normalized, and will instead contain \
        absolute count data (as in the original htcount .txt files). \
        Note that if save_csv is True, the saved .csv fill will contain ABSOLUTE COUNT DATA, as in the original \
        htcount .txt files, and NOT normalized data.
        :param save_reads_fname: str. Name under which to save the combined count data table. Does not need to include \
        the '.csv' suffix.
        :param save_not_counted_fname: save_reads_fname: str. Name under which to save the combined uncounted data. \
        Does not need to include the '.csv' suffix.
        :return:
        an HTCountFilter object containing the combined count data from all individual htcount .txt files in the \
        specified folder.
        """
        file_suffix = '.csv'
        if save_csv:
            assert isinstance(save_reads_fname, str)
            assert isinstance(save_not_counted_fname, str)

            if not save_reads_fname.endswith(file_suffix):
                save_reads_fname += file_suffix
            if not save_not_counted_fname.endswith(file_suffix):
                save_not_counted_fname += file_suffix

            save_reads_fname = f"{folder_path}\\{save_reads_fname}"
            save_not_counted_fname = f"{folder_path}\\{save_not_counted_fname}"

        folder = Path(folder_path)
        df = pd.DataFrame()
        for item in folder.iterdir():
            if item.is_file() and item.suffix == '.txt':
                df = pd.concat([df, pd.read_csv(item, sep='\t', index_col=0, names=[item.name])], axis=1)

        uncounted = df.loc['__no_feature':'__alignment_not_unique']
        counts = pd.concat([df, uncounted]).drop_duplicates(keep=False)

        if save_csv:
            general.save_to_csv(df=counts, filename=save_reads_fname, suffix='')
            general.save_to_csv(df=uncounted, filename=save_not_counted_fname, suffix='')

        fname = save_reads_fname if save_csv else folder.name + file_suffix
        h = HTCountFilter((Path(fname), counts))
        if norm_to_rpm:
            h.norm_reads_to_rpm(uncounted)
        return h

# TODO: a function that receives a dataframe, and can plot correlation with the BigTable instead of just enrichment
# TODO: add option for mask in clustergram
# TODO: fix no sample grouping in pca returning error
# TODO: function that prints all biotypes in the sample
