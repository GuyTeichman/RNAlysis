"""
This module can filter, normalize, intersect and visualize tabular data such as read counts and differential expression data.
Any tabular data saved in a csv format can be imported. \
Use this module to perform various filtering operations on your data,  normalize your data, \
perform set operations (union, intersection, etc), run basic exploratory analyses and plots \
(such as PCA, clustergram, violin plots, scatter, etc), \
save the filtered data to your computer, and more.
When you save filtered/modified data, its new file name will include by default \
 all of the operations performed on it, in the order they were performed, to allow easy traceback of your analyses.

"""

import os
import types
import warnings
from itertools import tee
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple, Type, Union, Callable

import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pairwisedist as pwdist
import seaborn as sns
from grid_strategy import strategies
from numba import jit
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances, silhouette_score

from rnalysis.utils import io, validation, preprocessing, ref_tables, clustering


class Filter:
    """
    An all-purpose Filter object.


    **Attributes**

    df: pandas DataFrame
        A DataFrame that contains the DESeq output file contents. \
        The DataFrame is modified upon usage of filter operations.
    shape: tuple (rows, columns)
        The dimensions of df.
    columns: list
        The columns of df.
    fname: pathlib.Path
        The path and filename for the purpose of saving df as a csv file. \
        Updates automatically when filter operations are applied.
    index_set: set
        All of the indices in the current DataFrame (which were not removed by previously used filter methods) \
        as a set.
    index_string: string
        A string of all feature indices in the current DataFrame separated by newline.
    """
    __slots__ = {'fname': 'filename with full path', 'df': 'pandas.DataFrame with the data', 'shape': '(rows, columns)'}

    def __init__(self, fname: Union[str, Path, tuple], drop_columns: Union[str, List[str]] = False):

        """
        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded DataFrame.
        :type drop_columns: str, list of str, or False (default False)

        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/counted.csv")

        """
        # init from a tuple (fname, DataFrame/Series). Used in self.__copy__(), self._inplace
        if isinstance(fname, tuple):
            assert isinstance(fname[1], (pd.DataFrame, pd.Series)) and isinstance(fname[0], (str, Path))
            self.fname = fname[0]
            self.df = fname[1]
        # init from a file (load the csv into a DataFrame/Series)
        else:
            assert isinstance(fname, (str, Path))
            self.fname = Path(fname)
            self.df = io.load_csv(fname, 0, squeeze=True, drop_columns=drop_columns)
        # check for duplicate indices
        if self.df.index.has_duplicates:
            warnings.warn("This Filter object contains multiple rows with the same name/index.")
        self.shape = self.df.shape

    def __repr__(self):
        return f"{type(self).__name__}('{self.fname}')"

    def __str__(self):
        return f"{type(self).__name__} of file {self.fname.stem}{self.fname.suffix}"

    def __len__(self):
        """
        Returns the number of rows in the Filter object

        """
        return self.shape[0]

    def __eq__(self, other):
        if self.df.equals(other.df) and self.shape == other.shape:
            return True
        return False

    def __contains__(self, item):
        return True if item in self.df.index else False

    def __copy__(self):

        """
        :rtype: Filter

        """
        return type(self)((self.fname, self.df.copy(deep=True)))

    @property
    def columns(self) -> list:
        """
        The columns of df.

        :return: a list of the columns in the Filter object.
        :rtype: list
        """
        return list(self.df.columns)

    def _inplace(self, new_df: pd.DataFrame, opposite: bool, inplace: bool, suffix: str,
                 printout_operation: str = 'filter'):

        """
        Executes the user's choice whether to filter in-place or create a new instance of the Filter object.

        :param new_df: the post-filtering DataFrame
        :type new_df: pd.DataFrame
        :param opposite: Determines whether to return the filtration ,or its opposite.
        :type opposite: bool
        :param inplace: Determines whether to filter in-place or not.
        :type inplace: bool
        :param suffix: The suffix to be added to the filename
        :type suffix: str
        :return: If inplace is False, returns a new instance of the Filter object.

        """
        assert isinstance(inplace, bool), "'inplace' must be True or False!"
        assert isinstance(opposite, bool), "'opposite' must be True or False!"
        assert printout_operation.lower() in ['filter', 'normalize', 'sort'], \
            f"Invalid input for variable 'printout_operation': {printout_operation}"
        # when user requests the opposite of a filter, return the Set Difference between the filtering result and self
        if opposite:
            new_df = self.df.loc[self.df.index.difference(new_df.index)]
            suffix += 'opposite'

        # update filename with the suffix of the operation that was just performed
        new_fname = Path(os.path.join(str(self.fname.parent), f"{self.fname.stem}{suffix}{self.fname.suffix}"))

        # generate printout for user ("Filtered X features, leaving Y... filtered inplace/not inplace")
        printout = ''
        if printout_operation.lower() == 'filter':
            printout += f"Filtered {self.df.shape[0] - new_df.shape[0]} features, leaving {new_df.shape[0]} " \
                        f"of the original {self.df.shape[0]} features. "
        elif printout_operation.lower() == 'normalize':
            printout += f"Normalized the values of {new_df.shape[0]} features. "
        elif printout_operation.lower() == 'sort':
            printout += f"Sorted {new_df.shape[0]} features. "
        # if inplace, modify the df, fname and shape properties of self
        if inplace:
            if printout_operation.lower() == 'filter':
                printout += 'Filtered inplace.'
            elif printout_operation.lower() == 'normalize':
                printout += 'Normalized inplace.'
            elif printout_operation.lower() == 'sort':
                printout += 'Sorted inplace.'
            print(printout)
            self.df, self.fname = new_df, new_fname
            self.shape = self.df.shape
        # if not inplace, copy self, modify the df/fname/shape properties of the copy, and return it
        else:
            if printout_operation.lower() == 'filter':
                printout += 'Filtering result saved to new object.'
            elif printout_operation.lower() == 'normalize':
                printout += 'Normalization result saved to a new object.'
            elif printout_operation.lower() == 'sort':
                printout += 'Sorting result saved to a new object.' \
                            '' \
                            ''
            print(printout)
            tmp_df, tmp_fname = self.df, self.fname
            self.df, self.fname = new_df, new_fname
            new_obj = self.__copy__()
            self.df, self.fname = tmp_df, tmp_fname
            return new_obj

    def save_csv(self, alt_filename: Union[None, str, Path] = None):

        """
        Saves the current filtered data to a .csv file.

        :param alt_filename: If None, file name will be generated automatically \
        according to the filtering methods used. \
        If it's a string, it will be used as the name of the saved file. Example input: 'myfilename'
        :type alt_filename: str, pathlib.Path, or None (default)

        """
        suffix = '.csv'
        # save with the default filename if no alternative filename was given
        if alt_filename is None:
            alt_filename = self.fname
        else:
            assert isinstance(alt_filename, (str, Path)), \
                f"'alt_filename' must be a string or Path object. Instead got {type(alt_filename)}."
            # make sure we don't add another suffix on top of an existing suffix
            if (isinstance(alt_filename, str) and alt_filename.endswith(suffix)) or \
                (isinstance(alt_filename, Path) and alt_filename.suffix == suffix):
                suffix = ''
            alt_filename = os.path.join(str(self.fname.parent), f"{alt_filename}{suffix}")
        io.save_csv(self.df, alt_filename)

    @staticmethod
    def _color_gen():

        """
        A generator that yields distinct colors up to a certain limit, and then yields randomized RGB values.

        :return: a color name string (like 'black', \
        or a numpy.ndarray of size (3,) containing three random values each between 0 and 1.

        """
        preset_colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink',
                         'tab:gray', 'tab:olive', 'tab:cyan', 'gold', 'maroon', 'mediumslateblue', 'fuchsia',
                         'mediumblue', 'black', 'lawngreen']
        for color in preset_colors:
            yield color
        while True:
            yield np.random.random(3)

    @staticmethod
    def _from_string(msg: str = '', delimiter: str = '\n'):

        """
        Takes a manual string input from the user, and then splits it using a delimiter into a list of values. \

        :param msg: a promprt to be printed to the user
        :param delimiter: the delimiter used to separate the values. Default is '\n'
        :return: A list of the comma-seperated values the user inserted.

        """
        string = input(msg)
        split = string.split(sep=delimiter)
        if split[-1] == '':
            split = split[:-1]
        return split

    def head(self, n: int = 5):

        """
        Return the first n rows of the Filter object. See pandas.DataFrame.head documentation.

        :type n: positive int, default 5
        :param n: Number of rows to show.
        :return: returns the first n rows of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> d.head()
                               baseMean  log2FoldChange  ...         pvalue           padj
            WBGene00000002  6820.755327        7.567762  ...   0.000000e+00   0.000000e+00
            WBGene00000003  3049.625670        9.138071  ...  4.660000e-302  4.280000e-298
            WBGene00000004  1432.911791        8.111737  ...  6.400000e-237  3.920000e-233
            WBGene00000005  4028.154186        6.534112  ...  1.700000e-228  7.800000e-225
            WBGene00000006  1230.585240        7.157428  ...  2.070000e-216  7.590000e-213
            <BLANKLINE>
            [5 rows x 6 columns]

            >>> d.head(3) # return only the first 3 rows
                               baseMean  log2FoldChange  ...         pvalue           padj
            WBGene00000002  6820.755327        7.567762  ...   0.000000e+00   0.000000e+00
            WBGene00000003  3049.625670        9.138071  ...  4.660000e-302  4.280000e-298
            WBGene00000004  1432.911791        8.111737  ...  6.400000e-237  3.920000e-233
            <BLANKLINE>
            [3 rows x 6 columns]

        """
        return self.df.head(n)

    def tail(self, n: int = 5) -> pd.DataFrame:

        """
        Return the last n rows of the Filter object. See pandas.DataFrame.tail documentation.

        :type n: positive int, default 5
        :param n: Number of rows to show.
        :rtype: pandas.DataFrame
        :return: returns the last n rows of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> d.tail()
                               baseMean  log2FoldChange  ...        pvalue          padj
            WBGene00000025  2236.185837        2.477374  ...  1.910000e-81  1.460000e-78
            WBGene00000026   343.648987       -4.037191  ...  2.320000e-75  1.700000e-72
            WBGene00000027   175.142856        6.352044  ...  1.580000e-74  1.120000e-71
            WBGene00000028   219.163200        3.913657  ...  3.420000e-72  2.320000e-69
            WBGene00000029  1066.242402       -2.811281  ...  1.420000e-70  9.290000e-68
            <BLANKLINE>
            [5 rows x 6 columns]


            >>> d.tail(8) # returns the last 8 rows
                               baseMean  log2FoldChange  ...        pvalue          padj
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

        """
        return self.df.tail(n)

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
        :return: If inplace is False, returns a new and filtered instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> # keep only the rows whose value in the column 'log2FoldChange' is below the 75th percentile
            >>> d.filter_percentile(0.75,'log2FoldChange')
            Filtered 7 features, leaving 21 of the original 28 features. Filtered inplace.

            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> # keep only the rows vulse value in the column 'log2FoldChange' is above the 25th percentile
            >>> d.filter_percentile(0.25,'log2FoldChange',opposite=True)
            Filtered 7 features, leaving 21 of the original 28 features. Filtered inplace.


        """
        assert isinstance(percentile, (float, int)) and 0 <= percentile <= 1, \
            "percentile must be a float between 0 and 1!"
        assert isinstance(column, str) and column in self.df, "Invalid column name!"
        suffix = f'_below{percentile}percentile'
        new_df = self.df[self.df[column] <= self.df[column].quantile(percentile)]
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_by_percentile(self, percentile: float, column: str) -> tuple:

        """
        Splits the Filter object into two Filter objects: \
        below and above the specified percentile in the spcfieid column.

        :type percentile: float between 0 and 1
        :param percentile: The percentile that all features above it will be filtered out.
        :type column: str
        :param column: Name of the DataFrame column according to which the filtering will be performed.
        :rtype: Tuple[filtering.Filter, filtering.Filter]
        :return: a tuple of two Filter objects: the first contains all of the features below the specified percentile, \
        and the second contains all of the features above and equal to the specified percentile.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> below, above = d.split_by_percentile(0.75,'log2FoldChange')
            Filtered 7 features, leaving 21 of the original 28 features. Filtering result saved to new object.
            Filtered 21 features, leaving 7 of the original 28 features. Filtering result saved to new object.

        """
        return self.filter_percentile(percentile=percentile, column=column, opposite=False,
                                      inplace=False), self.filter_percentile(percentile=percentile, column=column,
                                                                             opposite=True, inplace=False)

    def filter_biotype(self, biotype: Union[str, List[str]] = 'protein_coding',
                       ref: str = 'predefined', opposite: bool = False, inplace: bool = True):

        """
        Filters out all features that do not match the indicated biotype/biotypes. \
        Legal inputs: 'protein_coding','pseudogene','piRNA','miRNA','ncRNA','lincRNA','rRNA','snRNA','snoRNA'.

        :type biotype: string or list of strings
        :param biotype: the biotypes which will not be filtered out.
        :param ref: Name of the biotype reference file used to determine biotypes. \
        Default is the path defined by the user in the settings.yaml file.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows whose biotype is 'protein_coding'
            >>> counts.filter_biotype('protein_coding',ref='tests/biotype_ref_table_for_tests.csv')
            Filtered 9 features, leaving 13 of the original 22 features. Filtered inplace.

            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows whose biotype is 'protein_coding' or 'pseudogene'
            >>> counts.filter_biotype(['protein_coding','pseudogene'],ref='tests/biotype_ref_table_for_tests.csv')
            Filtered 0 features, leaving 22 of the original 22 features. Filtered inplace.

        """
        suffix = f"_{'_'.join(biotype)}"
        # make sure 'biotype' is a list of strings
        assert isinstance(biotype, (str, list)), "biotype must be a string or a list!"
        if isinstance(biotype, str):
            biotype = [biotype]
        # load the biotype reference table
        ref = ref_tables.get_biotype_ref_path(ref)
        ref_df = io.load_csv(ref)
        validation.validate_biotype_table(ref_df)
        ref_df.set_index('gene', inplace=True)
        # generate set of legal inputs
        legal_inputs = set(ref_df['biotype'].unique())
        # generate an empty boolean mask for filtering down the line
        mask = pd.Series(np.zeros_like(ref_df['biotype'], dtype=bool), index=ref_df['biotype'].index, name='biotype')
        # perform bitwise 'OR' between each biotype in the 'biotype' list and the mask
        for bio in biotype:
            assert bio in legal_inputs, f"biotype {bio} is not a legal string!"
            mask = mask | (ref_df['biotype'] == bio)
        # gene names which remain after filtering are the 'True' in the mask AND were previously in the DataFrame
        gene_names = ref_df[mask].index.intersection(self.df.index)
        new_df = self.df.loc[gene_names]
        return self._inplace(new_df, opposite, inplace, suffix)

    def filter_by_attribute(self, attributes: Union[str, List[str]] = None, mode='union', ref: str = 'predefined',
                            opposite: bool = False, inplace: bool = True):

        """
        Filters features according to user-defined attributes from an Attribute Reference Table. \
        When multiple attributes are given, filtering can be done in 'union' mode \
        (where features that belong to at least one attribute are not filtered out), or in 'intersection' mode \
        (where only features that belong to ALL attributes are not filtered out). \
        To learn more about user-defined attributes and Attribute Reference Tables, read the user guide.

        :type attributes: string or list of strings, \
        which are column titles in the user-defined Attribute Reference Table.
        :param attributes: attributes to filter by.
        :type mode: 'union' or 'intersection'.
        :param mode: If 'union', filters out every genomic feature that does not belong to one or more \
        of the indicated attributes. If 'intersection', \
        filters out every genomic feature that does not belong to ALL of the indicated attributes.
        :type ref: str or pathlib.Path (default 'predefined')
        :param ref: filename/path of the attribute reference table to be used as reference.
        :type opposite: bool (default False)
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new and filtered instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows that belong to the attribute 'attribute1'
            >>> counts.filter_by_attribute('attribute1',ref='tests/attr_ref_table_for_examples.csv')
            Filtered 15 features, leaving 7 of the original 22 features. Filtered inplace.

            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows that belong to the attributes 'attribute1' OR 'attribute3' (union)
            >>> counts.filter_by_attribute(['attribute1','attribute3'],ref='tests/attr_ref_table_for_examples.csv')
            Filtered 14 features, leaving 8 of the original 22 features. Filtered inplace.

            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows that belong to both attributes 'attribute1' AND 'attribute3' (intersection)
            >>> counts.filter_by_attribute(['attribute1','attribute3'],mode='intersection',
            ... ref='tests/attr_ref_table_for_examples.csv')
            Filtered 19 features, leaving 3 of the original 22 features. Filtered inplace.

            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows that DON'T belong to either 'attribute1','attribute3' or both
            >>> counts.filter_by_attribute(['attribute1','attribute3'],ref='tests/attr_ref_table_for_examples.csv',
            ... opposite=True)
            Filtered 8 features, leaving 14 of the original 22 features. Filtered inplace.

            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows that DON'T belong to both 'attribute1' AND 'attribute3'
            >>> counts.filter_by_attribute(['attribute1','attribute3'],mode='intersection',
            ... ref='tests/attr_ref_table_for_examples.csv',opposite=True)
            Filtered 3 features, leaving 19 of the original 22 features. Filtered inplace.

        """
        suffix = '_filtbyattr'
        # get attributes as input if none were supplied
        if attributes is None:
            attributes = self._from_string(
                "Please insert attributes separated by newline "
                "(for example: \n'epigenetic_related_genes\nnrde-3 targets\nALG-3/4 class small RNAs')")
        # make sure 'attributes' is a set/list/tuple of strings
        elif isinstance(attributes, str):
            attributes = [attributes]
        else:
            assert isinstance(attributes, (list, tuple, set))
        # make sure 'mode' is legal
        assert isinstance(mode, str), "'mode' must be a string!"
        mode = mode.lower()
        assert mode in {'union', 'intersection'}, \
            f"Illegal input {mode}: mode must be either 'union' or 'intersection'"
        # load the Attribute Reference Table
        ref = ref_tables.get_attr_ref_path(ref)
        attr_ref_table = io.load_csv(ref)
        validation.validate_attr_table(attr_ref_table)
        attr_ref_table.set_index('gene', inplace=True)
        # generate list of the indices that are positive for each attribute
        attr_indices_list = [attr_ref_table[attr_ref_table[attr].notnull()].index for attr in attributes]
        indices = pd.Index([])
        # if in union mode, calculate union between the indices of all attributes
        if mode == 'union':
            suffix += 'Union'
            for idx in attr_indices_list:
                indices = indices.union(idx)
            indices = indices.intersection(self.df.index)
        # if in intersection mode, calculate intersection between the indices of all attributes
        elif mode == 'intersection':
            suffix += 'Intersection'
            indices = self.df.index
            for idx in attr_indices_list:
                indices = indices.intersection(idx)

        new_df = self.df.loc[set(indices)]
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_by_attribute(self, attributes: Union[List[str], Tuple[str]], ref: str = 'predefined') -> tuple:

        """
        Splits the Filter object into multiple Filter objects, \
        each corresponding to one of the specified Attribute Reference Table attributes. \
        Each object contains only features that match its Attribute Reference Table attribute.

        :param attributes: list of attribute names from the Attribute Reference Table to filter by.
        :type attributes: list of strings
        :param ref: filename/path of the reference table to be used as reference.
        :rtype: Tuple[filtering.Filter]
        :return: A tuple of Filter objects, each containing only features that match one Attribute Reference Table attribute; \
        the Filter objects are returned in the same order the attributes were given in.


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> attribute1,attribute2 = counts.split_by_attribute(['attribute1','attribute2'],
            ... ref='tests/attr_ref_table_for_examples.csv')
            Filtered 15 features, leaving 7 of the original 22 features. Filtering result saved to new object.
            Filtered 20 features, leaving 2 of the original 22 features. Filtering result saved to new object.

        """
        assert isinstance(attributes, (list, tuple)), \
            f"'attributes' must be a list or a tuple. Got {type(attributes)} instead. "
        for attr in attributes:
            assert isinstance(attr, str), f"All attributes in 'split_by_attribute()' must be of type str. " \
                                          f"Attribute '{attr}' is of type {type(attr)}"
        return tuple([self.filter_by_attribute(att, mode='union', ref=ref, inplace=False) for att in attributes])

    def describe(self, percentiles: Union[list, Tuple, np.ndarray] = (0.01, 0.25, 0.5, 0.75, 0.99)):

        """
        Generate descriptive statistics that summarize the central tendency, dispersion and shape \
        of the dataset’s distribution, excluding NaN values. \
        For more information see the documentation of pandas.DataFrame.describe.

        :type percentiles: list-like of numbers, optional
        :param percentiles: The percentiles to include in the output. \
        All should fall between 0 and 1. \
        The default is [.25, .5, .75], which returns the 25th, 50th, and 75th percentiles.
        :return: Summary statistics of the dataset.
        :rtype: Series or DataFrame


        :Examples:
            >>> from rnalysis import filtering
            >>> import numpy as np
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> counts.describe()
                          cond1         cond2         cond3         cond4
            count     22.000000     22.000000     22.000000     22.000000
            mean    2515.590909   2209.227273   4230.227273   3099.818182
            std     4820.512674   4134.948493   7635.832664   5520.394522
            min        0.000000      0.000000      0.000000      0.000000
            1%         0.000000      0.000000      0.000000      0.000000
            25%        6.000000      6.250000      1.250000      0.250000
            50%       57.500000     52.500000     23.500000     21.000000
            75%     2637.000000   2479.000000   6030.500000   4669.750000
            99%    15054.950000  12714.290000  21955.390000  15603.510000
            max    15056.000000  12746.000000  22027.000000  15639.000000

            >>> # show the deciles (10%, 20%, 30%... 90%) of the columns
            >>> counts.describe(percentiles=np.arange(0.1, 1, 0.1))
                          cond1         cond2         cond3         cond4
            count     22.000000     22.000000     22.000000     22.000000
            mean    2515.590909   2209.227273   4230.227273   3099.818182
            std     4820.512674   4134.948493   7635.832664   5520.394522
            min        0.000000      0.000000      0.000000      0.000000
            10%        0.000000      0.200000      0.000000      0.000000
            20%        1.400000      3.200000      1.000000      0.000000
            30%       15.000000     15.700000      2.600000      1.000000
            40%       28.400000     26.800000     14.000000      9.000000
            50%       57.500000     52.500000     23.500000     21.000000
            60%       82.000000    106.800000     44.000000     33.000000
            70%      484.200000    395.500000    305.000000    302.500000
            80%     3398.600000   3172.600000   7981.400000   6213.000000
            90%     8722.100000   7941.800000  16449.500000  12129.900000
            max    15056.000000  12746.000000  22027.000000  15639.000000

        """
        return self.df.describe(percentiles=percentiles)

    @property
    def index_set(self) -> set:

        """
        Returns all of the features in the current DataFrame (which were not removed by previously used filter methods) \
        as a set. \
        if any duplicate features exist in the filter object (same WBGene appears more than once), \
        the corresponding WBGene index will appear in the returned set ONLY ONCE.

        :return: A set of WBGene names.


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> myset = counts.index_set
            >>> print(myset)
            {'WBGene00044022', 'WBGene00077504', 'WBGene00007079', 'WBGene00007069', 'WBGene00007063',
            'WBGene00007067', 'WBGene00077503', 'WBGene00007078', 'WBGene00007064', 'WBGene00077502', 'WBGene00044951',
            'WBGene00007077', 'WBGene00007066', 'WBGene00007076', 'WBGene00014997', 'WBGene00043990', 'WBGene00007074',
            'WBGene00043987', 'WBGene00007071', 'WBGene00043989', 'WBGene00043988', 'WBGene00007075'}

        """
        if self.df.index.has_duplicates:
            warnings.warn(" this filter object contains multiple rows with the same WBGene index. When "
                          "returning a set or string of features from this DESeqFilter object, each WBGene index will "
                          "appear ONLY ONCE!")
        return set(self.df.index)

    @property
    def index_string(self) -> str:
        r"""
        Returns a string of all feature indices in the current DataFrame, \
        sorted by their current order in the FIlter object, and separated by newline. \

        This includes all of the feature indices which were not filtered out by previously-used filter methods. \
         if any duplicate features exist in the filter object (same index appears more than once), \
        the corresponding index will appear in the returned string ONLY ONCE.

        :return: A string of WBGene indices separated by newlines (\\n). \
        For example, "WBGene00000001\\nWBGene00000003\\nWBGene12345678".


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> counts.sort(by='cond1',ascending=False)
            >>> mystring = counts.index_string
            >>> print(mystring)
            WBGene00007075
            WBGene00043988
            WBGene00043990
            WBGene00007079
            WBGene00007076
            WBGene00043989
            WBGene00007063
            WBGene00007077
            WBGene00007078
            WBGene00007071
            WBGene00007064
            WBGene00007066
            WBGene00007074
            WBGene00043987
            WBGene00007067
            WBGene00014997
            WBGene00044022
            WBGene00077503
            WBGene00077504
            WBGene00077502
            WBGene00007069
            WBGene00044951

        """
        return "\n".join(list(self.df.index))

    def print_features(self):

        """
        Print the feature indices in the Filter object, sorted by their current order in the FIlter object, \
        and separated by newline.


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> counts.sort(by='cond1',ascending=False)
            >>> counts.print_features()
            WBGene00007075
            WBGene00043988
            WBGene00043990
            WBGene00007079
            WBGene00007076
            WBGene00043989
            WBGene00007063
            WBGene00007077
            WBGene00007078
            WBGene00007071
            WBGene00007064
            WBGene00007066
            WBGene00007074
            WBGene00043987
            WBGene00007067
            WBGene00014997
            WBGene00044022
            WBGene00077503
            WBGene00077504
            WBGene00077502
            WBGene00007069
            WBGene00044951

        """
        print(self.index_string)

    def biotypes(self, return_format: str = 'short', ref: str = 'predefined'):

        """
        Returns a DataFrame of the biotypes in the Filter object and their count.
        :type return_format: 'short' or 'long' (default 'short')
        :param return_format: 'short' returns a short-form DataFrame, which states the biotypes \
        in the Filter object and their count. 'long' returns a long-form DataFrame,
        which also provides descriptive statistics of each column per biotype.
        :param ref: Name of the biotype reference table used to determine biotype. Default is ce11 (included in the package).
        :rtype: pandas.DataFrame
        :returns: a pandas DataFrame showing the number of values belonging to each biotype, \
        as well as additional descriptive statistics of format=='long'.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> # short-form view
            >>> d.biotypes(ref='tests/biotype_ref_table_for_tests.csv')
                            gene
            biotype
            protein_coding    26
            pseudogene         1
            unknown            1

            >>> # long-form view
            >>> d.biotypes(return_format='long', ref='tests/biotype_ref_table_for_tests.csv')
                           baseMean               ...           padj
                              count         mean  ...            75%            max
            biotype                               ...
            protein_coding     26.0  1823.089609  ...   1.005060e-90   9.290000e-68
            pseudogene          1.0  2688.043701  ...   1.800000e-94   1.800000e-94
            unknown             1.0  2085.995094  ...  3.070000e-152  3.070000e-152
            <BLANKLINE>
            [3 rows x 48 columns]

        """
        # validate 'return_format'
        assert isinstance(return_format, str), f"'return_format' must be a string!"
        return_format = return_format.lower()
        assert return_format in {'short', 'long'}, f"Invalid format '{return_format}'. Must be 'short' or 'long'."
        # load the Biotype Reference Table
        ref = ref_tables.get_biotype_ref_path(ref)
        ref_df = io.load_csv(ref)
        validation.validate_biotype_table(ref_df)
        # find which genes from tne Filter object don't appear in the Biotype Reference Table
        not_in_ref = self.df.index.difference(ref_df['gene'])
        if len(not_in_ref) > 0:
            warnings.warn(
                f'{len(not_in_ref)} of the features in the Filter object do not appear in the Biotype Reference Table. ')
            ref_df = ref_df.append(pd.DataFrame({'gene': not_in_ref, 'biotype': 'not_in_biotype_reference'}))
        # return just the number of genes/indices belonging to each biotype
        if return_format == 'short':
            return ref_df.set_index('gene', drop=False).loc[self.df.index].groupby('biotype').count()
        # additionally return descriptive statistics for each biotype
        elif return_format == 'long':
            self_df = self.df.__deepcopy__()
            self_df['biotype'] = ref_df.set_index('gene').loc[self.df.index]
            return self_df.groupby('biotype').describe()

    def number_filters(self, column: str, operator: str, value, opposite: bool = False, inplace: bool = True):

        """
        Applay a number filter (greater than, equal, lesser than) on a particular column in the Filter object.

        :type column: str
        :param column: name of the column to filter by
        :type operator: str: 'gt' / 'greater than' / '>', 'eq' / 'equals' / '=', 'lt' / 'lesser than' / '<'
        :param operator: the operator to filter the column by (greater than, equal or lesser than)
        :type value: number (int or float)
        :param value: the value to filter by
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> filt = filtering.Filter('tests/test_files/test_deseq.csv')
            >>> #keep only rows that have a value greater than 5900 in the column 'baseMean'.
            >>> filt.number_filters('baseMean','gt',5900)
            Filtered 26 features, leaving 2 of the original 28 features. Filtered inplace.

            >>> filt = filtering.Filter('tests/test_files/test_deseq.csv')
            >>> #keep only rows that have a value greater than 5900 in the column 'baseMean'.
            >>> filt.number_filters('baseMean','greater than',5900)
            Filtered 26 features, leaving 2 of the original 28 features. Filtered inplace.

            >>> filt = filtering.Filter('tests/test_files/test_deseq.csv')
            >>> #keep only rows that have a value greater than 5900 in the column 'baseMean'.
            >>> filt.number_filters('baseMean','>',5900)
            Filtered 26 features, leaving 2 of the original 28 features. Filtered inplace.

        """
        # determine whether operator is valid
        operator_dict = {'gt': 'gt', 'greater than': 'gt', '>': 'gt', 'eq': 'eq', 'equals': 'eq', '=': 'eq', 'lt': 'lt',
                         'lesser than': 'lt', '<': 'lt', 'equal': 'eq'}
        operator = operator.lower()
        assert operator in operator_dict, f"Invalid operator {operator}"
        op = operator_dict[operator]
        # determine that 'value' is a number
        assert isinstance(value, (int, float)), f"'value' must be a number!"
        # determine that the column is legal
        assert column in self.columns, f"column {column} not in DataFrame!"

        suffix = f"_{column}{op}{value}"
        # perform operation according to operator
        if op == 'eq':
            new_df = self.df[self.df[column] == value]
        elif op == 'gt':
            new_df = self.df[self.df[column] > value]
        elif op == 'lt':
            new_df = self.df[self.df[column] < value]

        # noinspection PyUnboundLocalVariable
        return self._inplace(new_df, opposite, inplace, suffix)

    def text_filters(self, column: str, operator: str, value: str, opposite: bool = False, inplace: bool = True):

        """
        Applay a text filter (equals, contains, starts with, ends with) on a particular column in the Filter object.

        :type column: str
        :param column: name of the column to filter by
        :type operator: str: 'eq' / 'equals' / '=', 'ct' / 'contains' / 'in', 'bw' / 'begins with', 'ew' / 'ends with'
        :param operator: the operator to filter the column by (equals, contains, starts with, ends with)
        :type value: number (int or float)
        :param value: the value to filter by
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> filt = filtering.Filter('tests/test_files/text_filters.csv')
            >>> # keep only rows that have a value that starts with 'AC3' in the column 'name'.
            >>> filt.text_filters('name','sw','AC3')
            Filtered 17 features, leaving 5 of the original 22 features. Filtered inplace.

        """
        # determine whether operator is valid
        operator_dict = {'eq': 'eq', 'equals': 'eq', '=': 'eq', 'ct': 'ct', 'in': 'ct', 'contains': 'ct', 'sw': 'sw',
                         'starts with': 'sw', 'ew': 'ew', 'ends with': 'ew', 'equal': 'eq', 'begins with': 'sw'}
        operator = operator.lower()
        assert operator in operator_dict, f"Invalid operator {operator}"
        op = operator_dict[operator]
        # determine that 'value' is a string
        assert isinstance(value, str), f"'value' must be a string!"
        # determine that the column is legal
        assert column in self.columns, f"column {column} not in DataFrame!"

        suffix = f"_{column}{op}{value}"
        # perform operation according to operator
        if op == 'eq':
            new_df = self.df[self.df[column] == value]
        elif op == 'ct':
            new_df = self.df[self.df[column].str.contains(value)]
        elif op == 'ew':
            new_df = self.df[self.df[column].str.endswith(value)]
        elif op == 'sw':
            new_df = self.df[self.df[column].str.startswith(value)]

        # noinspection PyUnboundLocalVariable
        return self._inplace(new_df, opposite, inplace, suffix)

    def filter_missing_values(self, columns: Union[str, List[str]] = 'all', opposite: bool = False,
                              inplace: bool = True):
        """
        Remove all rows whose values in the specified columns are missing (NaN).

        :type columns: str or list of str (default 'all')
        :param columns:name/names of the columns to check for missing values.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of the Filter object.

        :Examples:
            >>> from rnalysis import filtering
            >>> filt = filtering.Filter('tests/test_files/test_deseq_with_nan.csv')
            >>> filt_no_nan = filt.filter_missing_values(inplace=False)
            Filtered 3 features, leaving 25 of the original 28 features. Filtering result saved to new object.
            >>> filt_no_nan_basemean = filt.filter_missing_values(columns='baseMean', inplace=False)
            Filtered 1 features, leaving 27 of the original 28 features. Filtering result saved to new object.
            >>> filt_no_nan_basemean_pval = filt.filter_missing_values(columns=['baseMean','pval'], inplace=False)
            Filtered 2 features, leaving 26 of the original 28 features. Filtering result saved to new object.

        """
        subset = None
        suffix = '_removemissingvals'
        if columns == 'all':
            # on the rare off-chance that there is a column named 'all', decide not to decide
            if 'all' in self.columns:
                raise IndexError("Filter object contains a column named 'all'. RNAlysis cannot decide whether "
                                 "to filter based on the column 'all' or based on all columns. ")

        elif isinstance(columns, str):
            assert columns in self.columns, f"Column '{columns}' does not exist in the Filter object."
            if len(self.columns) > 1:
                subset = [columns]
                suffix += columns
        elif isinstance(columns, (list, tuple, set, np.ndarray)):
            for col in columns:
                assert isinstance(col, str), f"Column name {col} is of type {type(col)} instead of str. "
                assert col in self.columns, f"Column '{col}' does not exist in the Filter object."
                suffix += col
            subset = list(columns)
        else:
            raise TypeError(f"Invalid type for 'columns': {type(columns)}")
        if subset is not None:
            new_df = self.df.dropna(axis=0, subset=subset, inplace=False)
        else:
            new_df = self.df.dropna(axis=0, inplace=False)

        return self._inplace(new_df, opposite, inplace, suffix)

    def sort(self, by: Union[str, List[str]], ascending: Union[bool, List[bool]] = True, na_position: str = 'last',
             inplace: bool = True):

        """
        Sort the rows by the values of specified column or columns.

        :type by: str or list of str
        :param by: Names of the column or columns to sort by.
        :type ascending: bool or list of bool, default True
        :param ascending: Sort ascending vs. descending. Specify list for multiple sort orders. \
        If this is a list of bools, it must have the same length as 'by'.
        :type na_position: 'first' or 'last', default 'last'
        :param na_position: If 'first', puts NaNs at the beginning; if 'last', puts NaNs at the end.
        :type inplace: bool, default True
        :param inplace: If True, perform operation in-place. \
        Otherwise, returns a sorted copy of the Filter object without modifying the original.
        :return: None if inplace=True, a sorted Filter object otherwise.


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> counts.head()
                            cond1  cond2  cond3  cond4
            WBGene00007063    633    451    365    388
            WBGene00007064     60     57     20     23
            WBGene00044951      0      0      0      1
            WBGene00007066     55    266     46     39
            WBGene00007067     15     13      1      0
            >>> counts.sort(by='cond1',ascending=True)
            >>> counts.head()
                            cond1  cond2  cond3  cond4
            WBGene00044951      0      0      0      1
            WBGene00077504      0      0      0      0
            WBGene00007069      0      2      1      0
            WBGene00077502      0      0      0      0
            WBGene00077503      1      4      2      0

        """
        suffix = f'_sortedby{by}ascending{ascending}na{na_position}'
        new_df = self._sort(by=by, ascending=ascending, inplace=inplace, na_position=na_position)
        if inplace:
            new_df = self.df
        return self._inplace(new_df, False, inplace, suffix, 'sort')

    def _sort(self, by: Union[str, List[str]], ascending: Union[bool, List[bool]] = True, na_position: str = 'last',
              inplace: bool = True):
        """
        Sort the rows by the values of specified column or columns.

        :type by: str or list of str
        :param by: Names of the column or columns to sort by.
        :type ascending: bool or list of bool, default True
        :param ascending: Sort ascending vs. descending. Specify list for multiple sort orders. \
        If this is a list of bools, it must have the same length as 'by'.
        :type na_position: 'first' or 'last', default 'last'
        :param na_position: If 'first', puts NaNs at the beginning; if 'last', puts NaNs at the end.
        :type inplace: bool, default True
        :param inplace: If True, perform operation in-place. \
        Otherwise, returns a sorted copy of the Filter object without modifying the original.
        :return: None if inplace=True, a sorted Filter object otherwise.
        """

        return self.df.sort_values(by=by, axis=0, ascending=ascending, inplace=inplace, na_position=na_position)

    def filter_top_n(self, by: Union[str, List[str]], n: int = 100, ascending: Union[bool, List[bool]] = True,
                     na_position: str = 'last', opposite: bool = False, inplace: bool = True, ):

        """
        Sort the rows by the values of specified column or columns, then keep only the top 'n' rows.

        :type by: string or list of strings
        :param by: Names of the column or columns to sort and then filter by.
        :type n: int
        :param n: How many features to keep in the Filter object.
        :type ascending: bool or list of bools (default True)
        :param ascending: Sort ascending vs. descending. Specify list for multiple sort orders. \
        If this is a list of bools, it must have the same length as 'by'.
        :type na_position: 'first' or 'last', default 'last'
        :param na_position: If 'first', puts NaNs at the beginning; if 'last', puts NaNs at the end.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of Filter.


        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only the 10 rows with the highest values in the columns 'cond1'
            >>> counts.filter_top_n(by='cond1',n=10, ascending=False)
            Filtered 12 features, leaving 10 of the original 22 features. Filtered inplace.

            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only the 10 rows which have the lowest values in the columns 'cond1'
            >>> # and then the highest values in the column 'cond2'
            >>> counts.filter_top_n(by=['cond1','cond2'],n=10, ascending=[True,False])
            Filtered 12 features, leaving 10 of the original 22 features. Filtered inplace.

        """
        order = 'asc' if ascending else 'desc'
        suffix = f"_top{n}{by}{order}"
        assert isinstance(n, int), "n must be an integer!"
        assert n > 0, "n must be a positive integer!"
        if isinstance(by, list):
            for col in by:
                assert col in self.columns, f"{col} is not a column in the Filter object!"
        else:
            assert by in self.columns, f"{by} is not a column in the Filter object!"
        # sort the DataFrame by the specified column/columns, in the specified order
        self._sort(by=by, ascending=ascending, na_position=na_position, inplace=True)
        # keep only the top n values in the DataFrame after the sort (or the top len(self) items, if n>len(self))
        if n > self.df.shape[0]:
            warnings.warn(f'Current number of rows {self.df.shape[0]} is smaller than the specified n={n}. '
                          f'Therefore output Filter object will only have {self.df.shape[0]} rows. ')
        new_df = self.df.iloc[0:min(n, self.df.shape[0])]
        return self._inplace(new_df, opposite, inplace, suffix)

    @staticmethod
    def __return_type(index_set: set, return_type: str):
        assert isinstance(return_type, str), "'return_type' must be a string!!"
        if return_type == 'set':
            return index_set
        elif return_type == 'str':
            return "\n".join(index_set)
        else:
            raise ValueError(f"'return type' must be either 'set' or 'str', is instead '{return_type}'!")

    def _set_ops(self, others, return_type, op: Any):
        """
        Apply the supplied set operation (union/intersection/difference/symmetric difference) to the supplied objects.

        :param others: the other objects to apply the set operation to
        :type others: Filter or set objects.
        :param return_type: the return type of the output
        :type return_type: 'set' or 'str'
        :param op: the set operation
        :type op: function (set.union, set.intersection, set.difference or set.symmetric_difference)
        :return: a set/string of indices resulting from the set operation
        :rtype: set or str
        """
        others = list(others)
        for i, other in enumerate(others):
            if isinstance(other, Filter):
                others[i] = other.index_set
            elif isinstance(other, set):
                pass
            else:
                raise TypeError("'other' must be a Filter object or a set!")
        try:
            op_indices = op(set(self.df.index), *others)
        except TypeError as e:
            if op == set.symmetric_difference:
                raise TypeError(
                    f"Symmetric difference can only be calculated for two objects, {len(others) + 1} were given!")
            else:
                raise e
        return Filter.__return_type(op_indices, return_type)

    def intersection(self, *others, return_type: str = 'set', inplace: bool = False):

        """
        Keep only the features that exist in ALL of the given Filter objects/sets. \
        Can be done inplace on the first Filter object, or return a set/string of features.

        :type others: Filter or set objects.
        :param others: Objects to calculate intersection with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the intersecting WBGene indices. If 'str', returns a string of \
        the intersecting indices, delimited by a comma.
        :type inplace: bool, default False
        :param inplace: If True, filtering will be applied to the current Filter object. If False (default), \
        the function will return a set/str that contains the intersecting indices.
        :rtype: set or str
        :return: If inplace=False, returns a set/string of the WBGene indices that intersect between two Filter objects.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> a_set = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> # calculate intersection and return a set
            >>> d.intersection(a_set)
            {'WBGene00000002', 'WBGene00000003'}

            # calculate intersection and filter in-place
            >>> d.intersection(a_set, inplace=True)
            Filtered 26 features, leaving 2 of the original 28 features. Filtered inplace.

        """
        # if intersection is performed inplace, apply it to the self
        if inplace:
            suffix = f"_intersection"
            new_set = self._set_ops(others, 'set', set.intersection)
            return self._inplace(self.df.loc[new_set], opposite=False, inplace=inplace, suffix=suffix)
        # if intersection is not performed inplace, return a set/string according to user's request
        else:
            new_set = self._set_ops(others, return_type, set.intersection)
            return new_set

    def union(self, *others, return_type: str = 'set'):

        """
        Returns a set/string of the union of WBGene indices between multiple Filter objects \
        (the indices that exist in at least one of the Filter objects/sets).

        :type others: Filter or set objects.
        :param others: Objects to calculate union with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the union WBGene indices. If 'str', returns a string of \
        the union WBGene indices, delimited by a comma.
        :rtype: set or str
        :return:  a set/string of the WBGene indices that exist in at least one of the Filter objects.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # calculate union and return a set
            >>> d.union(counts)
            {'WBGene00000017', 'WBGene00000021', 'WBGene00044022', 'WBGene00077504', 'WBGene00000012',
            'WBGene00000024', 'WBGene00007079', 'WBGene00000010', 'WBGene00000020', 'WBGene00000005', 'WBGene00007069',
            'WBGene00007063', 'WBGene00007067', 'WBGene00077503', 'WBGene00007078', 'WBGene00000026', 'WBGene00000029',
            'WBGene00000002', 'WBGene00000003', 'WBGene00000006', 'WBGene00007064', 'WBGene00077502', 'WBGene00044951',
            'WBGene00000007', 'WBGene00000008', 'WBGene00000019', 'WBGene00007077', 'WBGene00000004', 'WBGene00007066',
            'WBGene00007076', 'WBGene00000013', 'WBGene00014997', 'WBGene00000023', 'WBGene00043990', 'WBGene00007074',
            'WBGene00000025', 'WBGene00000011', 'WBGene00043987', 'WBGene00007071', 'WBGene00000015', 'WBGene00000018',
            'WBGene00043989', 'WBGene00043988', 'WBGene00000014', 'WBGene00000016', 'WBGene00000027', 'WBGene00000028',
            'WBGene00007075', 'WBGene00000022', 'WBGene00000009'}

        """
        # union always returns a set/string (What is the meaning of in-place union between two different tables?)
        return self._set_ops(others, return_type, set.union)

    def difference(self, *others, return_type: str = 'set', inplace: bool = False):

        """
        Keep only the features that exist in the first Filter object/set but NOT in the others. \
        Can be done inplace on the first Filter object, or return a set/string of features.

        :type others: Filter or set objects.
        :param others: Objects to calculate difference with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the WBGene indices that exist only in the first Filter object. \
        If 'str', returns a string of the WBGene indices that exist only in the first Filter object, \
        delimited by a comma.
        :type inplace: bool, default False
        :param inplace: If True, filtering will be applied to the current Filter object. If False (default), \
        the function will return a set/str that contains the intersecting indices.
        :rtype: set or str
        :return: If inplace=False, returns a set/string of the WBGene indices\
         that exist only in the first Filter object/set (set difference).


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")
            >>> counts = filtering.CountFilter('tests/test_files/counted.csv')
            >>> a_set = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> # calculate difference and return a set
            >>> d.difference(counts, a_set)
            {'WBGene00007063', 'WBGene00007064', 'WBGene00007066', 'WBGene00007067', 'WBGene00007069', 'WBGene00007071',
             'WBGene00007074', 'WBGene00007075', 'WBGene00007076', 'WBGene00007077', 'WBGene00007078', 'WBGene00007079',
             'WBGene00014997', 'WBGene00043987', 'WBGene00043988', 'WBGene00043989', 'WBGene00043990', 'WBGene00044022',
             'WBGene00044951', 'WBGene00077502', 'WBGene00077503', 'WBGene00077504'}

            # calculate difference and filter in-place
            >>> d.difference(counts, a_set, inplace=True)
            Filtered 2 features, leaving 26 of the original 28 features. Filtered inplace.

        """
        # if difference is performed inplace, apply it to the self
        if inplace:
            suffix = f"_difference"
            new_set = self._set_ops(others, 'set', set.difference)
            return self._inplace(self.df.loc[new_set], opposite=False, inplace=inplace, suffix=suffix)
        # if difference is not performed inplace, return a set/string according to user's request
        else:
            new_set = self._set_ops(others, return_type, set.difference)
            return new_set

    def symmetric_difference(self, other, return_type: str = 'set'):

        """
        Returns a set/string of the WBGene indices that exist either in the first Filter object/set OR the second, \
        but NOT in both (set symmetric difference).

        :type other: Filter or set.
        :param other: a second Filter object/set to calculate symmetric difference with.
        :type return_type: 'set' or 'str.
        :param return_type: If 'set', returns a set of the WBGene indices that exist in exactly one Filter object. \
        If 'str', returns a string of the WBGene indices that exist in exactly one Filter object, delimited by a comma.
        :rtype: set or str
        :return: a set/string of the WBGene indices that that exist t in exactly one Filter. (set symmetric difference).


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter("tests/test_files/test_deseq.csv")
            >>> counts = filtering.CountFilter('tests/test_files/counted.csv')
            >>> # calculate difference and return a set
            >>> d.symmetric_difference(counts)
            {'WBGene00000017', 'WBGene00077504', 'WBGene00000024', 'WBGene00000010', 'WBGene00000020',
            'WBGene00007069', 'WBGene00007063', 'WBGene00007067', 'WBGene00007078', 'WBGene00000029', 'WBGene00000006',
            'WBGene00007064', 'WBGene00000019', 'WBGene00000004', 'WBGene00007066', 'WBGene00014997', 'WBGene00000023',
            'WBGene00007074', 'WBGene00000025', 'WBGene00043989', 'WBGene00043988', 'WBGene00000014', 'WBGene00000027',
            'WBGene00000021', 'WBGene00044022', 'WBGene00007079', 'WBGene00000012', 'WBGene00000005', 'WBGene00077503',
            'WBGene00000026', 'WBGene00000003', 'WBGene00000002', 'WBGene00077502', 'WBGene00044951', 'WBGene00007077',
            'WBGene00000007', 'WBGene00000008', 'WBGene00007076', 'WBGene00000013', 'WBGene00043990', 'WBGene00043987',
            'WBGene00007071', 'WBGene00000011', 'WBGene00000015', 'WBGene00000018', 'WBGene00000016', 'WBGene00000028',
            'WBGene00007075', 'WBGene00000022', 'WBGene00000009'}

        """
        # symmetric difference always returns a set/string
        # (What is the meaning of in-place symmetric difference between two different tables?)
        return self._set_ops([other], return_type, set.symmetric_difference)


class FoldChangeFilter(Filter):
    """
    A class that contains a single column, representing the gene-specific fold change between two conditions. \

     this class does not support 'inf' and '0' values, and importing a file with such values could lead \
    to incorrect filtering and statistical analyses.


    **Attributes**

    df: pandas Series
        A Series that contains the fold change values. \
        The Series is modified upon usage of filter operations.
    shape: tuple (rows, columns)
        The dimensions of df.
    columns: list
        The columns of df.
    fname: pathlib.Path
        The path and filename for the purpose of saving df as a csv file. \
        Updates automatically when filter operations are applied.
    index_set: set
        All of the indices in the current DataFrame (which were not removed by previously used filter methods) \
        as a set.
    index_string: string
        A string of all feature indices in the current DataFrame separated by newline.
    numerator: str
        Name of the numerator used to calculate the fold change.
    denominator: str
        Name of the denominator used to calculate the fold change.
    """
    __slots__ = {'numerator': 'name of the numerator', 'denominator': 'name of the denominator'}

    def __init__(self, fname: Union[str, Path, tuple], numerator_name: str, denominator_name: str):
        super().__init__(fname)
        self.numerator = numerator_name
        self.denominator = denominator_name
        self.df.name = 'Fold Change'
        # inf/0 can be problematic for functions down the line (like randomization test)
        if np.inf in self.df or 0 in self.df:
            warnings.warn(
                " FoldChangeFilter does not support 'inf' or '0' values! "
                "Unexpected results may occur during filtering or statistical analyses. ")

    def __repr__(self):
        return f"{type(self).__name__}('{self.fname}', '{self.numerator}', '{self.denominator}')"

    def __str__(self):
        return f"{type(self).__name__} (numerator: '{self.numerator}', denominator: '{self.denominator}') " \
               f"of file {self.fname.name}"

    def __copy__(self):
        return type(self)((self.fname, self.df.copy(deep=True)), numerator_name=self.numerator,
                          denominator_name=self.denominator)

    @property
    def columns(self) -> list:
        """
        The columns of df.

        :return: a list of the columns in the Filter object.
        :rtype: list
        """
        return [self.df.name]

    def randomization_test(self, ref, alpha: float = 0.05, reps=10000, save_csv: bool = False, fname=None,
                           random_seed: int = None):

        """
        Perform a randomization test to examine whether the fold change of a group of specific genomic features \
        is significantly different than the fold change of a background set of genomic features.

        :type ref: FoldChangeFilter
        :param ref: A reference FoldChangeFilter object which contains the fold change for every reference gene. \
        Will be used to calculate the expected score and to perform randomizations.
        :type alpha: float between 0 and 1
        :param alpha: Indicates the threshold for significance (alpha).
        :type reps: int larger than 0
        :param reps: How many repetitions to run the randomization for. \
        10,000 is the default. Recommended 10,000 or higher.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type random_seed: non-negative integer (default None)
        :type random_seed: The random seed used to initialize the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs.
        :rtype: pandas DataFrame
        :return: A Dataframe with the number of given genes, the observed fold change for the given group of genes, \
        the expected fold change for a group of genes of that size and the p value for the comparison.


        :Examples:
            >>> from rnalysis import filtering
            >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv' , 'numerator' , 'denominator')
            >>> f_background = f.filter_biotype('protein_coding', ref='tests/biotype_ref_table_for_tests.csv', inplace=False) #keep only protein-coding genes as reference
            Filtered 9 features, leaving 13 of the original 22 features. Filtering result saved to new object.
            >>> f_test = f_background.filter_by_attribute('attribute1', ref='tests/attr_ref_table_for_examples.csv', inplace=False)
            Filtered 6 features, leaving 7 of the original 13 features. Filtering result saved to new object.
            >>> rand_test_res = f_test.randomization_test(f_background)
            Calculating...
               group size  observed fold change  ...      pval  significant
            0           7              2.806873  ...  0.360264        False

            [1 rows x 5 columns]

        """
        # calculate observed and expected mean fold-change, and the set size (n)
        obs_fc = self.df.mean(axis=0)
        exp_fc = ref.df.mean()
        n = self.df.shape[0]
        # set random seed if requested
        if random_seed is not None:
            assert isinstance(random_seed, int) and random_seed >= 0, f"random_seed must be a non-negative integer. " \
                                                                      f"Value {random_seed} invalid."
            np.random.seed(random_seed)
        # run randomization test
        print('Calculating...')
        res = self._foldchange_randomization(ref.df.values, reps, obs_fc, exp_fc, n)
        # format the output DataFrame
        res_df = pd.DataFrame(res, columns=['group size', 'observed fold change', 'expected fold change', 'pval'],
                              index=[0])
        res_df['significant'] = res_df['pval'] <= alpha
        # save the output DataFrame if requested
        if save_csv:
            io.save_csv(res_df, fname)
        print(res_df)

        return res_df

    @staticmethod
    @jit(nopython=True)
    def _foldchange_randomization(vals: np.ndarray, reps: int, obs_fc: float, exp_fc: float, n: int):
        success = 0
        # determine the randomization test's direction (is observed greater/lesser than expected)
        if obs_fc > exp_fc:
            for _ in range(reps):
                # count the number of succeccess (mean(random subset) >= mean(observed)) over 'reps' repeats
                success += np.mean(np.random.choice(vals, n, replace=False)) >= obs_fc
        else:
            for _ in range(reps):
                # count the number of succeccess (mean(random subset) <= mean(observed)) over 'reps' repeats
                success += np.mean(np.random.choice(vals, n, replace=False)) <= obs_fc
        # calculate the positively-biased (less sig.) estimator of p-value (the unbiased estimator is (success / reps))
        # the positive bias corrects for the finite (and perhaps small) number of random permutations tested
        pval = (success + 1) / (reps + 1)
        return [[n, obs_fc, exp_fc, pval]]

    def filter_abs_log2_fold_change(self, abslog2fc: float = 1, opposite: bool = False, inplace: bool = True):

        """
        Filters out all features whose absolute log2 fold change is below the indicated threshold. \
        For example: if log2fc is 1.0, all features whose log2 fold change is between 1 and -1 (went up less than \
        two-fold or went down less than two-fold) will be filtered out.

        :param abslog2fc: The threshold absolute log2 fold change for filtering out a feature. Float or int. \
        All features whose absolute log2 fold change is lower than log2fc will be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current FoldChangeFilter object. If False, \
        the function will return a new FoldChangeFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of FoldChangeFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv','numerator name','denominator name')
            >>> f.filter_abs_log2_fold_change(2) # keep only rows whose log2(fold change) is >=2 or <=-2
            Filtered 18 features, leaving 4 of the original 22 features. Filtered inplace.

        """
        assert isinstance(abslog2fc, (float, int)), "abslog2fc must be a number!"
        assert abslog2fc >= 0, "abslog2fc must be non-negative!"
        suffix = f"_{abslog2fc}abslog2foldchange"
        new_df = self.df[np.abs(np.log2(self.df)) >= abslog2fc].dropna()
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
        :param inplace: If True (default), filtering will be applied to the current FoldChangeFilter object. If False, \
        the function will return a new FoldChangeFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of FoldChangeFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv','numerator name','denominator name')
            >>> # keep only rows with a positive log2(fold change) value
            >>> f.filter_fold_change_direction('pos')
            Filtered 10 features, leaving 12 of the original 22 features. Filtered inplace.

            >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv','numerator name','denominator name')
            >>>  # keep only rows with a negative log2(fold change) value
            >>> f.filter_fold_change_direction('neg')
            Filtered 14 features, leaving 8 of the original 22 features. Filtered inplace.

            >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv','numerator name','denominator name')
            >>> # keep only rows with a non-positive log2(fold change) value
            >>> f.filter_fold_change_direction('pos', opposite=True)
            Filtered 12 features, leaving 10 of the original 22 features. Filtered inplace.

        """
        assert isinstance(direction, str), \
            "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. "
        if direction == 'pos':
            new_df = self.df[self.df > 1]
            suffix = '_PositiveLog2FC'
        elif direction == 'neg':
            new_df = self.df[self.df < 1]
            suffix = '_NegativeLog2FC'
        else:
            raise ValueError(
                "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. ")
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_fold_change_direction(self) -> tuple:

        """
        Splits the features in the current FoldChangeFilter object into two complementary, non-overlapping \
        FoldChangeFilter objects, based on the direction of their log2(fold change). \
        The first object will contain only features with a positive log2(fold change), \
        the second object will contain only features with a negative log2(fold change). \
        Features with log2(fold change) = 0 will be ignored.

        :rtype: Tuple[filtering.FoldChangeFilter, filtering.FoldChangeFilter]
        :return: a tuple containing two FoldChangeFilter objects: the first has only features with positive log2 fold change, \
        and the other has only features with negative log2 fold change.


        :Examples:
            >>> from rnalysis import filtering
            >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv','numerator name','denominator name')
            >>> pos_log2fc, neg_log2fc = f.split_fold_change_direction()
            Filtered 10 features, leaving 12 of the original 22 features. Filtering result saved to new object.
            Filtered 14 features, leaving 8 of the original 22 features. Filtering result saved to new object.

        """
        return self.filter_fold_change_direction(direction='pos', inplace=False), \
               self.filter_fold_change_direction(direction='neg', inplace=False)


class DESeqFilter(Filter):
    """
    A class that receives a DESeq output file and can filter it according to various characteristics.

    **Attributes**

    df: pandas DataFrame
        A DataFrame that contains the DESeq output file contents. \
        The DataFrame is modified upon usage of filter operations.
    shape: tuple (rows, columns)
        The dimensions of df.
    columns: list
        The columns of df.
    fname: pathlib.Path
        The path and filename for the purpose of saving df as a csv file. \
        Updates automatically when filter operations are applied.
    index_set: set
        All of the indices in the current DataFrame (which were not removed by previously used filter methods) \
        as a set.
    index_string: string
        A string of all feature indices in the current DataFrame separated by newline.

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
        :return: If 'inplace' is False, returns a new instance of DESeqFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_significant(0.1) # keep only rows whose adjusted p-value is <=0.1
            Filtered 4 features, leaving 25 of the original 29 features. Filtered inplace.

             >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_significant(0.1, opposite=True) # keep only rows whose adjusted p-value is >0.1
            Filtered 25 features, leaving 4 of the original 29 features. Filtered inplace.

        """
        assert isinstance(alpha, float), "alpha must be a float!"
        new_df = self.df[self.df['padj'] <= alpha]
        suffix = f"_sig{alpha}"
        return self._inplace(new_df, opposite, inplace, suffix)

    def filter_abs_log2_fold_change(self, abslog2fc: float = 1, opposite: bool = False, inplace: bool = True):

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
        :return: If 'inplace' is False, returns a new instance of DESeqFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_abs_log2_fold_change(2) # keep only rows whose log2(fold change) is >=2 or <=-2
            Filtered 1 features, leaving 28 of the original 29 features. Filtered inplace.

        """
        assert isinstance(abslog2fc, (float, int)), "abslog2fc must be a number!"
        assert abslog2fc >= 0, "abslog2fc must be non-negative!"
        suffix = f"_{abslog2fc}abslog2foldchange"
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
        :return: If 'inplace' is False, returns a new instance of DESeqFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_fold_change_direction('pos') # keep only rows with a positive log2(fold change) value
            Filtered 3 features, leaving 26 of the original 29 features. Filtered inplace.

            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_fold_change_direction('neg') # keep only rows with a negative log2(fold change) value
            Filtered 27 features, leaving 2 of the original 29 features. Filtered inplace.

            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_fold_change_direction('pos', opposite=True) # keep only rows with a non-positive log2(fold change) value
            Filtered 26 features, leaving 3 of the original 29 features. Filtered inplace.

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

    def split_fold_change_direction(self) -> tuple:

        """
        Splits the features in the current DESeqFilter object into two complementary, non-overlapping DESeqFilter \
        objects, based on the direction of their log2foldchange. The first object will contain only features with a \
        positive log2foldchange, the second object will contain only features with a negative log2foldchange.

        :rtype: Tuple[filtering.DESeqFilter, filteirng.DESeqFilter]
        :return: a tuple containing two DESeqFilter objects: the first has only features with positive log2 fold change, \
        and the other has only features with negative log2 fold change.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/test_deseq.csv')
            >>> pos, neg = d.split_fold_change_direction()
            Filtered 2 features, leaving 26 of the original 28 features. Filtering result saved to new object.
            Filtered 26 features, leaving 2 of the original 28 features. Filtering result saved to new object.

        """
        return self.filter_fold_change_direction(direction='pos', inplace=False), self.filter_fold_change_direction(
            direction='neg', inplace=False)

    def volcano_plot(self, alpha: float = 0.1):

        """
        Plots a volcano plot (log2(fold change) vs -log10(adj. p-value)) of the DESeqFilter object. \
        Significantly upregulated features are colored in red, \
        and significantly downregulated features are colored in blue.

        :type alpha: float between 0 and 1
        :param alpha: the significance threshold to color data points as significantly up/down-regulated.

        .. figure::  volcano.png
           :align:   center
           :scale: 70 %

           Example plot of volcano_plot()

        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.style.use('seaborn-white')
        colors = pd.Series(index=self.df.index)
        colors.loc[(self.df['padj'] <= alpha) & (self.df['log2FoldChange'] > 0)] = 'tab:red'
        colors.loc[(self.df['padj'] <= alpha) & (self.df['log2FoldChange'] < 0)] = 'tab:blue'
        colors.fillna('grey', inplace=True)
        ax.scatter(self.df['log2FoldChange'], -np.log10(self.df['padj']), c=colors, s=1)
        ax.set_title(f"Volcano plot of {self.fname.stem}", fontsize=18)
        ax.set_xlabel('Log2(fold change)', fontsize=15)
        ax.set_ylabel('-Log10(adj. p-value)', fontsize=15)
        plt.show()
        return fig


class CountFilter(Filter):
    """
    A class that receives a count matrix and can filter it according to various characteristics.

    **Attributes**

    df: pandas DataFrame
        A DataFrame that contains the count matrix contents. \
        The DataFrame is modified upon usage of filter operations.
    shape: tuple (rows, columns)
        The dimensions of df.
    columns: list
        The columns of df.
    fname: pathlib.Path
        The path and filename for the purpose of saving df as a csv file. \
        Updates automatically when filter operations are applied.
    index_set: set
        All of the indices in the current DataFrame (which were not removed by previously used filter methods) \
        as a set.
    index_string: string
        A string of all feature indices in the current DataFrame separated by newline.
    triplicates: list
        Returns a nested list of the column names in the CountFilter, grouped by alphabetical order into triplicates. \
        For example, if counts.columns is ['A_rep1','A_rep2','A_rep3','B_rep1','B_rep2',_B_rep3'], then \
        counts.triplicates will be  [['A_rep1','A_rep2','A_rep3'],['B_rep1','B_rep2',_B_rep3']]

    """
    _precomputed_metrics = {'spearman': pwdist.spearman_distance, 'pearson': pwdist.pearson_distance,
                            'ys1': pwdist.ys1_distance, 'yr1': pwdist.yr1_distance,
                            'jackknife': pwdist.jackknife_distance}
    _transforms = {True: preprocessing.standard_box_cox, False: preprocessing.standardize}

    @property
    def triplicates(self):

        """
        Returns a nested list of the column names in the CountFilter, grouped by alphabetical order into triplicates. \
        For example, if counts.columns is ['A_rep1','A_rep2','A_rep3','B_rep1','B_rep2',_B_rep3'], then \
        counts.triplicates will be  [['A_rep1','A_rep2','A_rep3'],['B_rep1','B_rep2',_B_rep3']]

        """

        mltplr = 3
        triplicate = [self.columns[i * mltplr:(1 + i) * mltplr] for i in range(self.shape[1] // mltplr)]
        if len(self.columns[(self.shape[1] // mltplr) * mltplr::]) > 0:
            triplicate.append([self.columns[(self.shape[1] // mltplr) * mltplr::]])
            warnings.warn(
                f'Number of samples {self.shape[1]} is not divisible by 3. Appending the remaining {self.shape[1] % mltplr} as an inncomplete triplicate')
        return triplicate

    def fold_change(self, numerator, denominator, numer_name: str = 'default', denom_name: str = 'default'):

        """
        Calculate the fold change between the numerator condition and the denominator condition, \
        and return it as a FoldChangeFilter object.

        :type numerator: str, or list of strs
        :param numerator: the CountFilter columns to be used as the numerator. If multiple arguments are given \
        in a list, they will be averaged.
        :type denominator: str, or list of strs
        :param denominator: the CountFilter columns to be used as the denominator. If multiple arguments are given \
        in a list, they will be averaged.
        :type numer_name: str or 'default'
        :param numer_name: name to give the numerator condition. If 'default', the name will be generarated \
        automatically from the names of numerator columns.
        :type denom_name: str or 'default'
        :param denom_name: name to give the denominator condition. If 'default', the name will be generarated \
        automatically from the names of denominator columns.
        :rtype: FoldChangeFilter
        :return: A new instance of FoldChangeFilter


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted_fold_change.csv')
            >>> # calculate the fold change of mean(cond1_rep1,cond1_rep2)/mean(cond2_rep1,cond_2rep2)
            >>> f = c.fold_change(['cond1_rep1','cond1_rep2'],['cond2_rep1','cond2_rep2'])
            >>> f.numerator
            "Mean of ['cond1_rep1', 'cond1_rep2']"
            >>> f.denominator
            "Mean of ['cond2_rep1', 'cond2_rep2']"
            >>> type(f)
            rnalysis.filtering.FoldChangeFilter

        """
        assert isinstance(numer_name, str), "numerator name must be a string or 'default'!"
        assert isinstance(denom_name, str), "denominator name must be a string or 'default'!"
        numer_name = f"Mean of {numerator}" if numer_name == 'default' else numer_name
        denom_name = f"Mean of {denominator}" if denom_name == 'default' else denom_name
        # numerator and denominator should end up as a list of column names (strings)
        assert isinstance(numerator, (str, list, tuple)), "numerator must be a string or a list!"
        assert isinstance(denominator, (str, list, tuple)), "denominator must be a string or a list!"
        if isinstance(numerator, str):
            numerator = [numerator]
        elif isinstance(numerator, tuple):
            numerator = list(numerator)
        if isinstance(denominator, str):
            denominator = [denominator]
        elif isinstance(denominator, tuple):
            denominator = list(denominator)
        for num in numerator:
            assert num in self.df, f"all numerator arguments must be columns in the CountFilter object! ({num})"
        for den in denominator:
            assert den in self.df, f"all denominator arguments must be columns in the CountFilter object! ({den})"

        srs = (self.df[numerator].mean(axis=1) + 1) / (self.df[denominator].mean(axis=1) + 1)
        new_fname = Path(f"{str(self.fname.parent)}/{self.fname.stem}'_fold_change_'"
                         f"{numer_name}_over_{denom_name}_{self.fname.suffix}")
        # init the FoldChangeFilter object from an existing Series
        fcfilt = FoldChangeFilter((new_fname, srs), numerator_name=numer_name, denominator_name=denom_name)

        return fcfilt

        pass

    def pairplot(self, sample_list: list = 'all', log2: bool = False):

        """
        Plot pairwise relationships in the dataset. \
        Can plot both single samples and average multiple replicates. \
        For more information see the documentation of seaborn.pairplot.

        :type sample_list: 'all', list, or nested list.
        :param sample_list: A list of the sample names and/or grouped sample names to be included in the pairplot. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :type log2: bool
        :param log2: if True, the pairplot will be calculated with log2 of the dataframe, and not with the raw data. \
        If False (default), the pairplot will be calculated with the raw data.
        :return: An instance of seaborn.PairGrid.

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

    def _threshold_assertions(self, threshold: float = 1):

        """
        Assertions for functions that filter normalized values by some threshold, \
        or are meant to be used on normalized values.

        :param threshold: optional. A threshold value for filter_low_rpm to be asserted.

        """
        assert isinstance(threshold, (float, int)), "Threshold must be a number!"
        assert threshold >= 0, "Threshold must be zero or larger!"
        if '_norm' not in str(self.fname):
            warnings.warn("This function is meant for normalized values, and your values may not be normalized. ")

    def _avg_subsamples(self, sample_list: list):

        """
        Avarages subsamples/replicates according to the specified sample list. \
        Every member in the sample list should be either a name of a single sample (str), \
        or a list of multiple sample names to be averaged (list).

        :param sample_list: A list of the sample names and/or grouped sample names passed by the user. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6'] \
        and the resulting output will be a DataFrame containing the following columns: \
        ['SAMPLE1', 'SAMPLE2', 'SAMPLE3', 'SAMPLE6']
        :return: a pandas DataFrame containing samples/averaged subsamples according to the specified sample_list.

        """
        samples_df = pd.DataFrame()
        for sample in sample_list:
            if isinstance(sample, str):
                samples_df[sample] = self.df[sample].values
            elif isinstance(sample, (list, str, tuple)):
                samples_df[",".join(sample)] = self.df[sample].mean(axis=1).values

        return samples_df

    def normalize_to_rpm(self, special_counter_fname: str, inplace: bool = True):

        """
        Normalizes the reads in the CountFilter to Reads Per Million (RPM). \
        Uses a table of feature counts (ambiguous, no feature, not aligned, etc) from HTSeq's output. \
        Divides each column in the CountFilter object by (total reads + ambiguous + no feature)*10^-6 .

        :param special_counter_fname: the .csv file which contains feature information about the RNA library \
        (ambiguous, no feature, not aligned, etc).
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_to_rpm("tests/test_files/uncounted.csv")
            Normalized the values of 22 features. Normalized inplace.

        """
        suffix = '_normtorpm'
        new_df = self.df.copy()
        if isinstance(special_counter_fname, (str, Path)):
            features = io.load_csv(special_counter_fname, 0)
        elif isinstance(special_counter_fname, pd.DataFrame):
            features = special_counter_fname
        else:
            raise TypeError("Invalid type for 'special_counter_fname'!")
        for column in new_df.columns:
            norm_factor = (new_df[column].sum() + features.loc[r'__ambiguous', column] + features.loc[
                r'__no_feature', column] + features.loc[r'__alignment_not_unique', column]) / (10 ** 6)
            new_df[column] /= norm_factor
        return self._inplace(new_df, opposite=False, inplace=inplace, suffix=suffix, printout_operation='normalize')

    def normalize_with_scaling_factors(self, scaling_factor_fname: Union[str, Path], inplace: bool = True):

        """
        Normalizes the reads in the CountFilter using pre-calculated scaling factors. \
        Receives a table of sample names and their corresponding size factors, \
        and divides each column in the CountFilter by the corresponding scaling factor.

        :type scaling_factor_fname: str or pathlib.Path
        :param scaling_factor_fname: the .csv file which contains size factors for the different libraries.
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_with_scaling_factors("tests/test_files/scaling_factors.csv")
            Normalized the values of 22 features. Normalized inplace.

        """
        suffix = '_normwithscalingfactors'
        new_df = self.df.copy()
        if isinstance(scaling_factor_fname, (str, Path)):
            size_factors = io.load_csv(scaling_factor_fname)
        elif isinstance(scaling_factor_fname, pd.DataFrame):
            size_factors = scaling_factor_fname
        else:
            raise TypeError("Invalid type for 'scaling_factor_fname'!")
        for column in new_df.columns:
            norm_factor = size_factors[column].values
            new_df[column] /= norm_factor
        return self._inplace(new_df, opposite=False, inplace=inplace, suffix=suffix, printout_operation='normalize')

    def filter_low_reads(self, threshold: float = 5, opposite: bool = False, inplace: bool = True):

        """
        remove features which have less then 'threshold' reads all columns.

        :type threshold: float
        :param threshold: The minimal number of reads (counts, rpm, rpkm, tpm, etc) a feature should have \
        in at least one sample in order not to be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of CountFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> c.filter_low_reads(5) # remove all rows whose values in all columns are all <5
            Filtered 6 features, leaving 16 of the original 22 features. Filtered inplace.

        """
        self._threshold_assertions(threshold=threshold)
        new_df = self.df.loc[[True if max(vals) > threshold else False for gene, vals in self.df.iterrows()]]
        suffix = f"_filt{threshold}reads"
        return self._inplace(new_df, opposite, inplace, suffix)

    def split_by_reads(self, threshold: float = 5) -> tuple:

        """
        Splits the features in the current CountFilter object into two complementary, non-overlapping CountFilter \
        objects, based on the their maximum expression level. The first object will contain only highly-expressed \
         features (which have reads over the specified threshold in at least one sample). The second object will \
         contain only lowly-expressed features (which have reads below the specified threshold in all samples).

        :param threshold: The minimal number of reads (counts, RPM, RPKM, TPM etc) a feature needs to have \
        in at least one sample in order to be \
        included in the "highly expressed" object and no the "lowly expressed" object.
        :type threshold: float (default 5)
        :rtype: Tuple[filtering.CountFilter, filtering.CountFilter]
        :return: A tuple containing two CountFilter objects: the first has only highly-expressed features, \
        and the second has only lowly-expressed features.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> low_expr,high_expr = c.split_by_reads(5)
            Filtered 6 features, leaving 16 of the original 22 features. Filtering result saved to new object.
            Filtered 16 features, leaving 6 of the original 22 features. Filtering result saved to new object.

        """
        self._threshold_assertions(threshold=threshold)
        high_expr = self.df.loc[[True if max(vals) > threshold else False for gene, vals in self.df.iterrows()]]
        low_expr = self.df.loc[[False if max(vals) > threshold else True for gene, vals in self.df.iterrows()]]
        return self._inplace(high_expr, opposite=False, inplace=False, suffix=f'_below{threshold}reads'), self._inplace(
            low_expr, opposite=False, inplace=False, suffix=f'_above{threshold}reads')

    def filter_by_row_sum(self, threshold: float = 5, opposite: bool = False, inplace: bool = True):

        """
        Removes features/rows whose sum is belove 'threshold'.

        :type threshold: float
        :param threshold: The minimal sum a row should have in order not to be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of CountFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> c.filter_by_row_sum(5) # remove all rows whose sum is <5
            Filtered 4 features, leaving 18 of the original 22 features. Filtered inplace.

        """
        self._threshold_assertions(threshold=threshold)
        new_df = self.df.loc[self.df.sum(axis=1) >= threshold]
        suffix = f"_filt{threshold}sum"
        return self._inplace(new_df, opposite, inplace, suffix)

    def _plot_clustering(self, n_clusters: int, data, labels: np.ndarray, centers: np.ndarray, title: str,
                         plot_style: str, split_plots: bool):
        assert plot_style.lower() in {'all', 'std_bar', 'std_area'}, \
            f"Invalid value for 'plot_style': '{plot_style}'. 'plot_style' must be 'all', 'std_bar' or 'std_area'. "
        assert isinstance(split_plots, bool), f"'split_plots' must be of type bool, instead got {type(split_plots)}"
        if split_plots:
            subplots = [plt.subplots(constrained_layout=True, figsize=(6.4, 5.6)) for _ in range(n_clusters)]
            ylabel_fontsize = 14
            ticklabels_fontsize = 9
        else:
            grid = strategies.SquareStrategy()
            subplots = grid.get_grid(n_clusters)
            plt.close()
            fig = plt.figure(figsize=(14, 9))
            fig.suptitle(title, fontsize=18)
            ylabel_fontsize = 10
            ticklabels_fontsize = 5
        axes = []
        min_y, max_y = 0, 0
        color_generator = CountFilter._color_gen()
        for i, subplot in enumerate(subplots):
            if split_plots:
                axes.append(subplot[1])
                fig = subplot[0]
                fig.suptitle(title, fontsize=18)
            else:
                axes.append(fig.add_subplot(subplot))
            mean = centers[i, :]
            stdev = data[labels == i, :].T.std(axis=1, ddof=1)
            x = np.arange(self.shape[1]) + 0.5
            color = next(color_generator)
            max_y = max(max_y, np.max(stdev + mean))
            min_y = min(min_y, np.min(mean - stdev))
            axes[-1].axhline(color='grey', linewidth=2, linestyle='--')
            if plot_style.lower() == 'std_area':
                axes[-1].plot(x, mean, marker='o', linewidth=2, color=color, markeredgecolor='black')
                axes[-1].fill_between(x, mean + stdev, mean - stdev, alpha=0.2, color=color)
            else:
                axes[-1].errorbar(x, mean, marker='o', linewidth=2, color=color, yerr=stdev, markeredgecolor='black')
                if plot_style.lower() == 'all':
                    vals = data[labels == i, :].T
                    axes[-1].plot(x, vals, color=color, alpha=0.05, linewidth=0.35)

            axes[-1].set_title(f"Cluster number {i + 1} ({np.count_nonzero(labels == i)} genes)")
            axes[-1].set_ylabel('Standardized expression', fontsize=ylabel_fontsize)
            axes[-1].set_xticks(x)
            axes[-1].set_xticklabels(self.columns, fontsize=ticklabels_fontsize)
        for ax in axes:
            ax.set_ylim((min_y, max_y))
        if not split_plots:
            fig.tight_layout(rect=[0, 0.03, 1, 0.92])
        plt.show()
        return fig, axes

    def _gap_statistic(self, clusterer_class: type, transform: Callable, random_state: int, clusterer_kwargs: dict,
                       n_refs: int = 10, max_clusters: int = 20):
        # determine the range of k values to be tested
        k_range = np.arange(1, max_clusters + 1)
        print(f"Calculating optimal K using the Gap Statistic method in range {2}:{max_clusters}...")
        # set the random state, if one was supplied
        if random_state is not None:
            np.random.seed(random_state)
        # calculate the SVD of the observed data (X), and transform via X' = x_tag = dot(X, V)
        # note: the data is centered by sklearn.decomposition.PCA, no need to pre-center it.
        pca = PCA(random_state=random_state).fit(self.df.values)
        x_tag = pca.transform(self.df.values)
        # determine the ranges of the columns of X' (x_tag)
        a, b = x_tag.min(axis=0, keepdims=True), x_tag.max(axis=0, keepdims=True)
        # transform the observed data using Box-Cox, and then standardize it
        data = transform(self.df.values)
        # generate 'n_refs' random reference arrays:
        # draw uniform features Z' over the ranges of the columns of X', and back-transform via Z = dot(Z', V.T), then
        # transform the random reference data using Box-Cox, and then standardize it
        refs = [transform(pca.inverse_transform(np.random.random_sample(size=self.df.values.shape) * (b - a) + a))
                for _ in range(n_refs)]
        # allocate empty arrays for observed/expected log(inertia), gap scores Gap(K) and gap error S(K)
        log_disp_obs = np.zeros((len(k_range)))
        log_disp_exp = np.zeros((len(k_range)))
        gap_scores = np.zeros((len(k_range)))
        gap_error = np.zeros((len(k_range)))

        # iterate over K values in the range
        for ind, this_k in enumerate(k_range):
            # init the clusterer with given arguments
            clusterer = clusterer_class(n_clusters=this_k, **clusterer_kwargs)
            # cluster each of the n_refs reference arrays into Ki clusters
            ref_disps = np.zeros(n_refs)
            for ref_ind, ref in enumerate(refs):
                # calculate dispersion/inertia for reference data clustered into Ki clusters
                ref_disps[ref_ind] = self._compute_dispersion(clusterer.fit(ref), ref, this_k)
            # calculate dispersion/inertia for the observed data clustered into Ki clusters
            disp = self._compute_dispersion(clusterer.fit(data), data, this_k)
            # calculate the mean of the log(reference dispersion) values
            log_disp_exp[ind] = np.mean(np.log(ref_disps))
            # calculate the log(observed dispersion) value
            log_disp_obs[ind] = np.log(disp)
            # calculate gap score for K = Ki
            gap_scores[ind] = log_disp_exp[ind] - log_disp_obs[ind]
            # calculate the gap standard deviation Sd(K) for K = Ki
            stdev = np.sqrt(np.mean((np.log(ref_disps) - log_disp_exp[ind]) ** 2.0))
            # calculate the gap standard error (Sk) for K = Ki
            gap_error[ind] = stdev * np.sqrt(1 + (1 / n_refs))

        # Find the smallest K such that Gap(K) >= Gap(K+1) - S(K+1), AKA Diff(K) = Gap(k) - Gap(k+1) + S(K+1) >= 0
        diff = gap_scores[:-1] - gap_scores[1::] + gap_error[1::]
        # Return all K values that fulfilled the condition as 'potentially good K values'
        k_candidates = k_range[:-1][diff >= 0]
        # if no K fulfilled the condition (the gap function is monotonically increasing),
        # return the last value of K as the best K
        best_k = k_candidates[0] if len(k_candidates) > 0 else k_range[-1]
        best_k_ind = np.argmax(k_range == best_k)
        print(f"Using the Gap Statistic method, {best_k} was chosen as the best number of clusters (K).")
        if len(k_candidates) >= 2:
            print(f"Other potentially good values of K that were found: {list(k_candidates[1:])}. ")

        # plot results
        fig = self._plot_gap(k_range, log_disp_obs, log_disp_exp, gap_scores, gap_error, best_k, best_k_ind)

        return best_k, fig

    @staticmethod
    def _compute_dispersion(clusterer, data, k):
        # known as Wk in the Gap Statistic paper
        if hasattr(clusterer, 'inertia_'):
            return clusterer.inertia_
        labels = clusterer.labels_
        dispersion = 0
        for r in range(k):
            this_label = labels == r
            n = np.sum(this_label)
            # don't calculate within-cluster dispersion for empty clusters or clusters with a single member
            if n > 1:
                dispersion += np.sum(pairwise_distances(data[this_label, :]) ** 2) / (2 * n)
        return dispersion

    @staticmethod
    def _plot_gap(k_range, log_disp_obs, log_disp_exp, gap_scores, gap_error, best_k, best_k_ind):
        fig, (ax_inertia, ax) = plt.subplots(1, 2, figsize=(14, 9))
        ax_inertia.plot(k_range, log_disp_obs, '-o')
        ax_inertia.plot(k_range, log_disp_exp, '-o')
        ax_inertia.legend(['Observed', 'Expected'])
        ax_inertia.set_ylabel(r"$\ln$(dispersion)", fontsize=15)
        ax_inertia.set_xlabel("Number of clusters (K)", fontsize=15)
        ax_inertia.set_xticks(k_range)
        ax.errorbar(k_range, gap_scores, yerr=gap_error, marker='o', color='r')
        ax.set_title("Gap Statistic method for optimal K selection")
        ax.set_ylabel('Gap Value', fontsize=15)
        ax.set_xlabel("Number of clusters (K)", fontsize=15)
        ax.annotate(f'Best K={best_k}', xy=(best_k, (gap_scores[best_k_ind] - gap_error[best_k_ind]) / 1.05),
                    xytext=(best_k, (gap_scores[best_k_ind] - gap_error[best_k_ind]) / 1.2),
                    arrowprops=dict(facecolor='black'))
        ax.set_xticks(k_range)
        sns.despine()
        plt.show()
        return fig

    def _silhouette(self, clusterer_class: type, transform: Callable, clusterer_kwargs: dict, max_clusters: int = 20):
        """

        :param clusterer_class:
        :type clusterer_class:
        :param clusterer_kwargs:
        :type clusterer_kwargs:
        :param max_clusters:
        :type max_clusters:
        :return:
        :rtype:
        """
        print(f"Calculating optimal k using the Silhouette method in range {2}:{max_clusters}...")
        data = transform(self.df.values)
        sil_scores = []
        k_range = np.arange(2, max_clusters + 1)
        for n_clusters in k_range:
            clusterer = clusterer_class(n_clusters=n_clusters, **clusterer_kwargs)
            sil_scores.append(silhouette_score(data, clusterer.fit_predict(data)))
        fig = plt.figure(figsize=(7, 9))
        ax = fig.add_subplot(111)
        ax.plot(k_range, sil_scores, '--o', color='r')
        ax.set_ylim(-1.0, 1.0)
        ax.set_title("Silhouette method for optimal k selection")
        ax.set_ylabel('Average silhouette score')
        ax.set_xlabel("Number of clusters (k)")
        ax.axhline(0, color='grey', linewidth=2, linestyle='--')
        best_k = k_range[int(np.argmax(sil_scores))]
        ax.annotate(f'Best k={best_k}', xy=(best_k, max(sil_scores)), xytext=(best_k + 1, 0.75),
                    arrowprops=dict(facecolor='black', shrink=0.15))
        ax.set_xticks(k_range)
        print(f"Using the Silhouette method, {best_k} was chosen as the best number of clusters (k).")
        return best_k, fig

    def _parse_k(self, k: Union[int, List[int], str], clusterer_class: type, transform: Callable,
                 random_state: int, clusterer_kwargs: dict, max_clusters: int):

        max_clusters = min(20, self.shape[0] // 4) if max_clusters == 'default' else max_clusters
        if isinstance(k, str) and k.lower() == 'silhouette':
            best_k, _ = self._silhouette(clusterer_class=clusterer_class, transform=transform,
                                         clusterer_kwargs=clusterer_kwargs, max_clusters=max_clusters)
            k = [best_k]
        elif isinstance(k, str) and k.lower() == 'gap':
            best_k, _ = self._gap_statistic(clusterer_class=clusterer_class, transform=transform,
                                            random_state=random_state, clusterer_kwargs=clusterer_kwargs,
                                            max_clusters=max_clusters)
            k = [best_k]
        else:
            if not isinstance(k, Iterable):
                k = [k]
            k, k_copy = tee(k)

            assert np.all([isinstance(item, int) and item >= 2 for item in k_copy]), \
                f"Invalid value for k: '{k}'. k must be an integer>=2, Iterable of integers, 'gap', or 'silhouette'. "
        return k

    def _clustering_get_transform(self, power_transform: bool, metric: str) -> Tuple[Callable, str]:
        if metric in self._precomputed_metrics:
            def transform(x):
                return self._precomputed_metrics[metric](self._transforms[power_transform](x))

            return_metric = 'precomputed'
        else:
            transform = self._transforms[power_transform]
            return_metric = metric
        return transform, return_metric

    @staticmethod
    def _clustering_assert_metric(metric: str, linkage: str = None):
        legal_metrics = {'euclidean', 'spearman', 'pearson', 'manhattan', 'cosine', 'ys1', 'yr1', 'jackknife'}
        legal_linkages = {'single', 'average', 'complete', 'ward'}
        assert isinstance(metric, str), f"'metric' must be a string. Instead got '{type(metric)}'."
        metric = metric.lower()
        assert metric in legal_metrics, f"'metric' must be one of {legal_metrics}. Instead got '{metric}'."

        if linkage is not None:
            assert isinstance(linkage, str), f"'linkage' must be a string. Instead got '{type(linkage)}'."
            linkage = linkage.lower()
            assert linkage in legal_linkages, f"'linkage' must be in {legal_linkages}. Instead got '{linkage}'."
            return metric, linkage
        return metric

    def split_kmeans(self, k: Union[int, List[int], str], random_state: int = None, n_init: int = 3,
                     max_iter: int = 300, power_transform: bool = False,
                     plot_style: str = 'all', split_plots: bool = False, max_clusters: int = 'default'):
        # get the transform function
        transform, _ = self._clustering_get_transform(power_transform, 'euclidean')
        # set keyworded arguments for clusterer type = KMeans
        clusterer_kwargs = dict(init='k-means++', random_state=random_state, n_init=n_init, max_iter=max_iter)
        # parse k value
        k = self._parse_k(k=k, clusterer_class=KMeans, transform=transform, random_state=random_state,
                          clusterer_kwargs=clusterer_kwargs, max_clusters=max_clusters)
        # generate standardized data for plots
        data_for_plot = transform(self.df)
        # iterate over all K values, generating clustering results and plotting them
        filt_obj_tuples = []
        for this_k in k:
            this_k = int(this_k)
            clusterer = KMeans(n_clusters=this_k, **clusterer_kwargs).fit(transform(self.df.values))
            # plot results
            self._plot_clustering(n_clusters=this_k, data=data_for_plot, labels=clusterer.labels_,
                                  centers=clusterer.cluster_centers_,
                                  title=f"Results of K-Means Clustering for K={this_k} and "
                                        f"power_transform={power_transform}", plot_style=plot_style,
                                  split_plots=split_plots)
            # split the CountFilter object
            filt_obj_tuples.append(
                tuple([self._inplace(self.df.loc[clusterer.labels_ == i], opposite=False, inplace=False,
                                     suffix=f'_kmeanscluster{i + 1}') for i in range(this_k)]))
        # if only a single K was calculated, don't return it as a list of length
        return filt_obj_tuples[0] if len(filt_obj_tuples) == 1 else filt_obj_tuples

    def split_hierarchical(self, n_clusters: Union[int, List[int], str, None], metric: str = 'euclidean',
                           linkage: str = 'average', power_transform: bool = False, distance_threshold: float = None,
                           plot_style: str = 'all', split_plots: bool = False,
                           max_clusters: int = 'default', gap_random_state: int = None):

        # assert metrics/linkages are legal
        metric, linkage = self._clustering_assert_metric(metric, linkage)
        # set keyworded arguments for clusterer type = AgglomerativeClustering
        clusterer_kwargs = dict(affinity=metric, linkage=linkage)
        # get the transform/pairwise-distance function
        transform, clusterer_kwargs['affinity'] = self._clustering_get_transform(power_transform, metric)
        # generate standardized data for plots
        data_for_plot = preprocessing.standardize(self.df)
        # if k was supplied, parse k value
        if n_clusters is not None:
            k = self._parse_k(k=n_clusters, clusterer_class=AgglomerativeClustering, transform=transform,
                              random_state=gap_random_state, clusterer_kwargs=clusterer_kwargs,
                              max_clusters=max_clusters)
        # if k was not supplied, determined k by AgglomerativeClustering with distance_threshold
        else:
            clusterer = AgglomerativeClustering(n_clusters=None, distance_threshold=distance_threshold,
                                                **clusterer_kwargs).fit(transform(self.df.values))
            k = clusterer.n_clusters_
            # calculate cluster centers
            centers = np.array([data_for_plot[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(k)])
            # plot results
            self._plot_clustering(n_clusters=k, data=data_for_plot, labels=clusterer.labels_, centers=centers,
                                  title=f"Results of Hierarchical Clustering for n_clusters={k}, "
                                        f"metric='{metric}', \nlinkage='{linkage}' and "
                                        f"power_transform={power_transform}", plot_style=plot_style,
                                  split_plots=split_plots)
            # split the CountFilter object
            return tuple([self._inplace(self.df.loc[clusterer.labels_ == i], opposite=False, inplace=False,
                                        suffix=f'_hierarchicalcluster{i + 1}') for i in range(k)])
        # iterate over all K values, generating clustering results and plotting them
        filt_obj_tuples = []
        for this_k in k:
            this_k = int(this_k)
            clusterer = AgglomerativeClustering(n_clusters=this_k, **clusterer_kwargs).fit(transform(self.df.values))
            # calculate cluster centers
            centers = np.array([data_for_plot[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(this_k)])
            # plot results
            self._plot_clustering(n_clusters=this_k, data=data_for_plot, labels=clusterer.labels_, centers=centers,
                                  title=f"Results of Hierarchical Clustering for n_clusters={this_k}, "
                                        f"metric='{metric}', \nlinkage='{linkage}' and "
                                        f"power_transform={power_transform}", plot_style=plot_style,
                                  split_plots=split_plots)
            # split the CountFilter object
            filt_obj_tuples.append(
                tuple([self._inplace(self.df.loc[clusterer.labels_ == i], opposite=False, inplace=False,
                                     suffix=f'_hierarchicalcluster{i + 1}') for i in range(this_k)]))
        # if only a single K was calculated, don't return it as a list of length
        return filt_obj_tuples[0] if len(filt_obj_tuples) == 1 else filt_obj_tuples

    def split_kmedoids(self, k: Union[int, List[int], str], random_state: int = None, n_init: int = 3,
                       max_iter: int = 300, metric: str = 'euclidean', power_transform: bool = False,
                       plot_style: str = 'all', split_plots: bool = False, max_clusters: int = 'default'):
        # assert metric is legal
        metric = self._clustering_assert_metric(metric)
        # set keyworded arguments for clusterer type = KMedoidsIter
        clusterer_kwargs = dict(init='k-medoids++', metric=metric, random_state=random_state, n_init=n_init,
                                max_iter=max_iter)
        # get the transform/pairwise-distance function
        transform, clusterer_kwargs['metric'] = self._clustering_get_transform(power_transform, metric)
        # generate standardized data for plots
        data_for_plot = preprocessing.standardize(self.df)
        # parse K value
        k = self._parse_k(k=k, clusterer_class=clustering.KMedoidsIter, transform=transform, random_state=random_state,
                          clusterer_kwargs=clusterer_kwargs, max_clusters=max_clusters)
        # iterate over all K values, generating clustering results and plotting them
        filt_obj_tuples = []
        for this_k in k:
            this_k = int(this_k)
            clusterer = clustering.KMedoidsIter(n_clusters=this_k, **clusterer_kwargs).fit(transform(self.df.values))
            # get cluster centers
            centers = data_for_plot[clusterer.medoid_indices_, :]
            # plot results
            self._plot_clustering(n_clusters=this_k, data=data_for_plot, labels=clusterer.labels_, centers=centers,
                                  title=f"Results of K-Medoids Clustering for K={this_k}, metric='{metric}', "
                                        f"power_transform={power_transform}",
                                  plot_style=plot_style, split_plots=split_plots)
            # split the CountFilter object
            filt_obj_tuples.append(
                tuple([self._inplace(self.df.loc[clusterer.labels_ == i], opposite=False, inplace=False,
                                     suffix=f'_kmedoidscluster{i + 1}') for i in range(this_k)]))
        # if only a single K was calculated, don't return it as a list of length
        return filt_obj_tuples[0] if len(filt_obj_tuples) == 1 else filt_obj_tuples

    def split_hdbscan(self, min_cluster_size: int = 5, min_samples: int = 1, metric='euclidean',
                      cluster_selection_epsilon: float = 0, cluster_selection_method: str = 'eom',
                      power_transform: bool = False,
                      plot_style: str = 'all', split_plots: bool = False, return_prob: bool = False):

        assert isinstance(min_cluster_size, int) and min_cluster_size >= 2, \
            f"'min_cluster_size' must be an integer >=2. Instead got {min_cluster_size}, type={type(min_cluster_size)}."
        assert min_cluster_size <= self.shape[0], \
            "'min_cluster_size' cannot be larger than the number of features in the Filter object. "
        assert isinstance(cluster_selection_method, str), \
            f"'cluster_selection_method' must be a string. Instead got {type(cluster_selection_method)}."
        assert isinstance(metric, str), f"'metric' must be a string. Instead got {type(metric)}."
        cluster_selection_method = cluster_selection_method.lower()
        metric = metric.lower()
        transform, metric = self._clustering_get_transform(power_transform, metric)
        data = transform(self.df.values)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                                    cluster_selection_epsilon=cluster_selection_epsilon, metric=metric,
                                    cluster_selection_method=cluster_selection_method)
        clusterer.fit(data)
        n_clusters = clusterer.labels_.max() + 1
        probabilities = clusterer.probabilities_
        unclustered = np.count_nonzero(clusterer.labels_ == -1)
        if n_clusters == 0:
            print("Found 0 clusters with the given parameters. Please try again with different parameters. ")
        else:
            print(f"Found {n_clusters} clusters of average size "
                  f"{(len(clusterer.labels_) - unclustered) / n_clusters  :.2f}. "
                  f"Number of unclustered genes is {unclustered}, "
                  f"which are {100 * (unclustered / len(clusterer.labels_)) :.2f}% of the genes.")
            means = np.array([data[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(n_clusters)])
            self._plot_clustering(n_clusters=n_clusters, data=data, labels=clusterer.labels_, centers=means,
                                  title=f"Results of HDBSCAN Clustering for min_cluster_size={min_cluster_size}, "
                                        f"min_samples = {min_samples}, metric='{metric}', \n"
                                        f"epsilon={cluster_selection_epsilon}, "
                                        f"method='{cluster_selection_method}'and "
                                        f"power_transform={power_transform}", plot_style=plot_style,
                                  split_plots=split_plots)

        filt_objs = tuple([self._inplace(self.df.loc[clusterer.labels_ == i], opposite=False, inplace=False,
                                         suffix=f'_hdbscancluster{i + 1}') for i in range(n_clusters)])
        return [filt_objs, probabilities] if return_prob else filt_objs

    def clustergram(self, sample_names: list = 'all', metric: str = 'euclidean', linkage: str = 'average'):

        """
        Performs hierarchical clustering and plots a clustergram on the base-2 log of a given set of samples.

        :type sample_names: 'all' or list.
        :param sample_names: the names of the relevant samples in a list. \
        Example input: ["condition1_rep1", "condition1_rep2", "condition1_rep3", \
        "condition2_rep1", "condition3_rep1", "condition3_rep2"]
        :type metric: 'euclidean', 'hamming', 'correlation', or any other \
        distance metric available in scipy.spatial.distance.pdist
        :param metric: the distance metric to use in the clustergram. \
        For all possible inputs and their meaning see scipy.spatial.distance.pdist documentation online.
        :type linkage: 'single', 'average', 'complete', 'weighted', 'centroid', 'median' or 'ward'.
        :param linkage: the linkage method to use in the clustergram. \
        For all possible inputs and their meaning see scipy.cluster.hierarchy.linkage documentation online.

        :return: A seaborn clustermap object.


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
                                    cmap=sns.color_palette("RdBu_r", 12), yticklabels=False)
        plt.show()
        return clustering

    def plot_expression(self, features: Union[List[str], str],
                        sample_grouping: Union[Dict[str, List[int]], Dict[str, List[str]]],
                        count_unit: str = 'Reads per million'):

        """
        Plot the average expression and standard error of the specified features under the specified conditions.

        :type features: str or list of strings
        :param features: the feature/features to plot expression for.
        :type sample_grouping: dict, with condition names as keys \
        and list of the sample numbers or names for each condition as a list
        :param sample_grouping: a dictionary of the conditions to plot expression for. \
        Each key should be a name of a conditions, and the value for each key is \
        a list of the numbers of columns to be used as samples of that condition. \
        For example, if the first 3 columns are replicates of the condition 'condition 1' and \
        the last 3 column are replicates of the condition 'condition 2', then sample_grouping should be: \
        {'condition 1':[0, 1, 2], 'condition 2':[3, 4, 5]}
        :type count_unit: str, default 'Reads per million'
        :param count_unit: The unit of the count data. Will be displayed in the y axis.

        .. figure::  plot_expression.png
           :align:   center
           :scale: 40 %

           Example plot of plot_expression()

        """
        plt.style.use('seaborn-white')
        if isinstance(features, str):
            features = [features]
        assert isinstance(features, list), "'features' must be a string or list of strings!"

        g = strategies.SquareStrategy()
        subplots = g.get_grid(len(features))
        plt.close()
        fig = plt.figure()
        axes = []
        ylims = []
        for subplot, feature in zip(subplots, features):
            axes.append(fig.add_subplot(subplot))
            mean = [self.df.loc[feature].iloc[ind].mean() if np.all([isinstance(i, int) for i in ind]) else
                    self.df.loc[feature][ind].mean()
                    for ind in sample_grouping.values()]
            sem = [self.df.loc[feature].iloc[ind].sem() if np.all([isinstance(i, int) for i in ind]) else
                   self.df.loc[feature][ind].sem() for
                   ind in sample_grouping.values()]
            axes[-1].bar(np.arange(len(sample_grouping)), mean, yerr=sem)
            axes[-1].set_xticks(np.arange(len(sample_grouping)))
            axes[-1].set_xticklabels(list(sample_grouping.keys()))
            axes[-1].set_title(feature)
            plt.ylabel(count_unit)
            sns.despine()
            ylims.append(axes[-1].get_ylim()[1])
        for ax in axes:
            ax.set_ylim((0.0, max(ylims)))
        fig.tight_layout()
        plt.show()
        return fig

    def pca(self, sample_names: list = 'all', n_components=3, sample_grouping: list = None, labels: bool = True):

        """
        Performs Principal Component Analysis (PCA), visualizing the principal components that explain the most\
         variance between the different samples. The function will automatically plot Principal Component #1 \
         with every other Principal Components calculated.

        :type sample_names: 'all' or list.
        :param sample_names: the names of the relevant samples in a list. \
        Example input: ["1_REP_A", "1_REP_B", "1_REP_C", "2_REP_A", "2_REP_B", "2_REP_C", "2_REP_D", "3_REP_A"]
        :type n_components: positive int (default 3)
        :param n_components: number of PCA components to return.
        :type sample_grouping: list of positive integers, 'triplicates' or None (default)
        :param sample_grouping: Optional. Indicates which samples are grouped together as replicates, \
        so they will be colored similarly in the PCA plot. A list of indices from 1 and up, that indicates the sample \
         grouping. \
         For example, if sample_names is: \
        ["1_REP_A", "1_REP_B", "1_REP_C", "2_REP_A", "2_REP_B", "2_REP_C", "2_REP_D", "3_REP_A"], \
        then the sample_grouping will be: \
        [1, 1, 1, 2, 2, 2, 2, 3]. \
        If 'triplicate', then sample_groupins will automatically group samples into triplicates. For example: \
        [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4].
        :type labels: bool (default True)
        :param labels: if True, labels the points on the PCA plot.
        :return: A tuple whose first element is an sklearn.decomposition.pca object, \
        and second element is a list of matplotlib.axis objects.

        .. figure::  pca.png
           :align:   center
           :scale: 40 %

           Example plot of pca()

        """
        assert isinstance(n_components, int) and n_components >= 2, \
            f"'n_components' must be an integer >=2. Instead got {n_components}."
        if sample_names == 'all':
            sample_names = list(self.df.columns)
            data = self.df.transpose()
        else:
            data = self.df[sample_names].transpose()
        data_standardized = preprocessing.standardize(data)

        pca_obj = PCA(n_components=n_components)
        pcomps = pca_obj.fit_transform(data_standardized)
        columns = [f'Principal component {i + 1}' for i in range(n_components)]
        principal_df = pd.DataFrame(data=pcomps, columns=columns)
        final_df = principal_df
        final_df['lib'] = pd.Series(sample_names)

        pc_var = pca_obj.explained_variance_ratio_
        graphs = 2 if n_components > 2 else 1
        axes = []
        for graph in range(graphs):
            axes.append(CountFilter._plot_pca(
                final_df=final_df[['Principal component 1', f'Principal component {2 + graph}', 'lib']],
                pc1_var=pc_var[0], pc2_var=pc_var[1 + graph], sample_grouping=sample_grouping, labels=labels))

        return pca_obj, axes

    @staticmethod
    def _plot_pca(final_df: pd.DataFrame, pc1_var: float, pc2_var: float, sample_grouping: list, labels: bool):

        """
        Internal method, used to plot the results from CountFilter.pca().

        :param final_df: The DataFrame output from pca
        :param pc1_var: Variance explained by the first PC.
        :param pc2_var: Variance explained by the second PC.
        :param sample_grouping: a list of indices from 0 and up, that indicates what samples are grouped together as \
        biological replicates. For example, if sample_names is: \
        ["1A_N2_25", "1B_N2_25", "1C_N2_25", "2A_rde4_25", "2B_rde4_25", "2C_rde4_25"], \
        then the sample_grouping will be: \
        [0,0,0,1,1,1]
        :return: an axis object containing the PCA plot.

        """
        plt.style.use('seaborn-whitegrid')

        if sample_grouping is None:
            sample_grouping = [i + 1 for i in range(final_df.shape[0])]
        elif sample_grouping == 'triplicate' or sample_grouping == 'triplicates':
            sample_grouping = [1 + i // 3 for i in range(final_df.shape[0])]

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.grid(True)
        ax.set_xlabel(f'{final_df.columns[0]} (explained {pc1_var * 100 :.2f}%)', fontsize=15)
        ax.set_ylabel(f'{final_df.columns[1]} (explained {pc2_var * 100 :.2f}%)', fontsize=15)
        ax.set_title('PCA', fontsize=20)

        color_generator = CountFilter._color_gen()
        color_opts = [next(color_generator) for _ in range(max(sample_grouping))]
        colors = [color_opts[i - 1] for i in sample_grouping]

        ax.scatter(final_df.iloc[:, 0], final_df.iloc[:, 1], c=colors, s=50)
        if labels:
            for _, row in final_df.iterrows():
                row[0] += 1
                row[1] += 1
                ax.text(*row)
        ax.grid(True)
        return ax

    def scatter_sample_vs_sample(self, sample1: Union[str, List[str]], sample2: Union[str, List[str]],
                                 xlabel: str = None, ylabel: str = None, title: str = None, highlight=None):

        """
        Generate a scatter plot where every dot is a feature, the x value is log10 of reads \
        (counts, RPM, RPKM, TPM, etc) in sample1, the y value is log10 of reads in sample2.

        :param sample1: Name of the first sample from the CountFilter object. \
        If sample1 is a list, they will be avarged as replicates.
        :type sample1: string or list of strings
        :param sample2: Name of the second sample from the CountFilter object. \
        If sample2 is a list, they will be averaged as replicates.
        :type sample2: string or list of strings
        :param xlabel: optional. If not specified, sample1 will be used as xlabel.
        :type xlabel: str
        :param ylabel: optional. If not specified, sample2 will be used as ylabel.
        :type ylabel: str
        :param title: optional. If not specified, a title will be generated automatically.
        :type title: str
        :param highlight: If specified, the points in the scatter corresponding to the names/features in 'highlight' \
        will be highlighted in red.
        :type highlight: Filter object or iterable of strings
        :return: a matplotlib axis object.

        .. figure::  rpm_vs_rpm.png
           :align:   center
           :scale: 60 %

           Example plot of scatter_sample_vs_sample()


        """
        self._threshold_assertions()
        assert isinstance(sample1, (str, list, tuple, set)) and isinstance(sample2, (str, list, tuple, set))
        xvals = np.log10(self.df[sample1].values + 1) if isinstance(sample1, str) else np.log10(
            self.df[sample1].mean(axis=1).values + 1)
        yvals = np.log10(self.df[sample2].values + 1) if isinstance(sample2, str) else np.log10(
            self.df[sample2].mean(axis=1).values + 1)

        plt.style.use('seaborn-whitegrid')
        if xlabel is None:
            xlabel = f'log10(reads per million + 1) from library {sample1}'
        if ylabel is None:
            ylabel = f'log10(reads per million + 1) from sample {sample2}'
        if title is None:
            title = f'{sample1} vs {sample2}'
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(xlabel, fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        ax.set_title(title, fontsize=17)
        ax.scatter(xvals, yvals, s=3, c='#6d7178')

        if highlight is not None:
            highlight_features = highlight.index_set if issubclass(highlight.__class__,
                                                                   Filter) else highlight
            highlight_intersection = highlight_features.intersection(self.index_set)
            if len(highlight_intersection) < len(highlight_features):
                warnings.warn(
                    f'Out of {len(highlight_features)} features to be highlighted, '
                    f'{len(highlight_features) - len(highlight_intersection)} features are missing from the CountFilter '
                    f'object and will not be highlighted.')
                highlight_features = highlight_intersection

            xvals_highlight = np.log10(self.df[sample1].loc[highlight_features].values + 1) if \
                isinstance(sample1, str) else np.log10(self.df[sample1].loc[highlight_features].mean(axis=1).values + 1)
            yvals_highlight = np.log10(self.df[sample2].loc[highlight_features].values + 1) if \
                isinstance(sample2, str) else np.log10(self.df[sample2].loc[highlight_features].mean(axis=1).values + 1)

            ax.scatter(xvals_highlight, yvals_highlight, s=3, c=np.array([[0.75, 0.1, 0.1]]))
        plt.show()
        return ax

    def box_plot(self, samples='all', notch: bool = True, scatter: bool = False, ylabel: str = 'log10(RPM + 1)'):

        """
        Generates a box plot of the specified samples in the CountFilter object in log10 scale. \
        Can plot both single samples and average multiple replicates. \
        It is recommended to use this function on normalized values and not on absolute read values. \
        The box indicates 25% and 75% percentiles, and the white dot indicates the median.

        :type samples: 'all' or list.
        :param samples: A list of the sample names and/or grouped sample names to be plotted. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :type notch: bool (default True)
        :param notch: if True, adds a confidence-interval notch to the box-plot.
        :type scatter: bool (default False)
        :param scatter: if True, adds a scatter-plot on top of the box-plot.
        :type ylabel: str (default 'Log10(RPM + 1)')
        :param ylabel: the label of the Y axis.
        :return: a seaborn boxplot object.

        .. figure::  ???.png
           :align:   center
           :scale: 60 %

           Example plot of box_plot()


        """
        self._threshold_assertions()
        if samples == 'all':
            samples_df = self.df
        else:
            samples_df = self._avg_subsamples(samples)

        samples_df = np.log10(samples_df + 1)
        _ = plt.figure(figsize=(8, 8))

        box = sns.boxplot(data=np.log10(samples_df + 1), notch=notch)
        if scatter:
            _ = sns.stripplot(data=np.log10(samples_df + 1), color='gray', size=2)
        plt.style.use('seaborn-whitegrid')
        plt.xlabel("Samples")
        plt.ylabel(ylabel)
        plt.show()
        return box

    def enhanced_box_plot(self, samples='all', scatter: bool = False, ylabel: str = 'log10(RPM + 1)'):

        """
        Generates an enhanced box-plot of the specified samples in the CountFilter object in log10 scale. \
        Can plot both single samples and average multiple replicates. \
        It is recommended to use this function on normalized values and not on absolute read values. \
        The box indicates 25% and 75% percentiles, and the white dot indicates the median.

        :type samples: 'all' or list.
        :param samples: A list of the sample names and/or grouped sample names to be plotted. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :type scatter: bool (default False)
        :param scatter: if True, adds a scatter-plot on top of the box-plot.
        :type ylabel: str (default 'Log10(RPM + 1)')
        :param ylabel: the label of the Y axis.
        :return: a seaborn enhanced box-plot object.

        .. figure::  ???.png
           :align:   center
           :scale: 60 %

           Example plot of enhanced_box_plot()

        """
        self._threshold_assertions()
        if samples == 'all':
            samples_df = self.df
        else:
            samples_df = self._avg_subsamples(samples)

        samples_df = np.log10(samples_df + 1)
        _ = plt.figure(figsize=(8, 8))

        boxen = sns.boxenplot(data=samples_df)
        if scatter:
            _ = sns.stripplot(data=samples_df, color='gray', size=2)
        plt.style.use('seaborn-whitegrid')
        plt.xlabel("Samples")
        plt.ylabel(ylabel)
        plt.show()
        return boxen

    def violin_plot(self, samples='all', ylabel: str = 'log10(RPM + 1)'):

        """
        Generates a violin plot of the specified samples in the CountFilter object in log10 scale. \
        Can plot both single samples and average multiple replicates. \
        It is recommended to use this function on normalized values and not on absolute read values. \
        Box inside the violin plot indicates 25% and 75% percentiles, and the white dot indicates the median.

        :type samples: 'all' or list.
        :param samples: A list of the sample names and/or grouped sample names to be plotted. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :type ylabel: str (default 'Log10(RPM + 1)')
        :param ylabel: the label of the Y axis.
        :return: a seaborn violin object.

        .. figure::  violin.png
           :align:   center
           :scale: 60 %

           Example plot of violin_plot()


        """
        self._threshold_assertions()
        if samples == 'all':
            samples_df = self.df
        else:
            samples_df = self._avg_subsamples(samples)

        samples_df = np.log10(samples_df + 1)
        _ = plt.figure(figsize=(8, 8))

        violin = sns.violinplot(data=samples_df)
        plt.style.use('seaborn-whitegrid')
        plt.xlabel("Samples")
        plt.ylabel(ylabel)
        plt.show()
        return violin

    # TODO: add ranksum test

    @staticmethod
    def from_folder(folder_path: str, norm_to_rpm: bool = False, save_csv: bool = False, counted_fname: str = None,
                    uncounted_fname: str = None, input_format: str = '.txt'):

        """
        Iterates over HTSeq count .txt files in a given folder and combines them into a single CountFilter table. \
        Can also save the count data table and the uncounted data table to .csv files, and normalize the CountFilter \
        table to reads per million (RPM). Note that the saved data will always be count data, and not normalized data, \
        regardless if the CountFilter table was normalized or not.

        :param folder_path: str or pathlib.Path. Full path of the folder that contains individual htcount .txt files.
        :param norm_to_rpm: bool. If True, the CountFilter table will be automatically normalized to reads per \
        million (RPM). If False (defualt), the CountFilter object will not be normalized, and will instead contain \
        absolute count data (as in the original htcount .txt files). \
        Note that if save_csv is True, the saved .csv fill will contain ABSOLUTE COUNT DATA, as in the original \
        htcount .txt files, and NOT normalized data.
        :param save_csv: bool. If True, the joint DataFrame of count data and uncounted data will be saved \
        to two separate .csv files. The files will be saved in 'folder_path', and named according to the parameters \
        'counted_fname' for the count data, and 'uncounted_fname' for the uncounted data (unaligned, \
        alignment not unique, etc).
        :param counted_fname: str. Name under which to save the combined count data table. Does not need to include \
        the '.csv' suffix.
        :param uncounted_fname: counted_fname: str. Name under which to save the combined uncounted data. \
        Does not need to include the '.csv' suffix.
        :param input_format: the file format of the input files. Default is '.txt'.
        :return: an CountFilter object containing the combined count data from all individual htcount .txt files in the \
        specified folder.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder')

            >>> c = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder', norm_to_rpm=True) # This will also normalize the CountFilter to reads-per-million (RPM).
            Normalized the values of 10 features. Normalized inplace.

            >>> c = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder', save_csv=True, counted_fname='name_for_reads_csv_file', uncounted_fname='name_for_uncounted_reads_csv_file') # This will also save the counted reads and uncounted reads as separate .csv files

        """
        file_suffix = '.csv'
        if save_csv:
            assert isinstance(counted_fname, str)
            assert isinstance(uncounted_fname, str)

            if not counted_fname.endswith(file_suffix):
                counted_fname += file_suffix
            if not uncounted_fname.endswith(file_suffix):
                uncounted_fname += file_suffix

            counted_fname = os.path.join(folder_path, counted_fname)
            uncounted_fname = os.path.join(folder_path, uncounted_fname)

        folder = Path(folder_path)
        df = pd.DataFrame()
        for item in folder.iterdir():
            if item.is_file() and item.suffix == input_format:
                df = pd.concat([df, pd.read_csv(item, sep='\t', index_col=0, names=[item.stem])], axis=1)
        assert not df.empty, f"Error: no valid files with the suffix '{input_format}' were found in '{folder_path}'."

        uncounted = df.loc[
            ['__no_feature', '__ambiguous', '__alignment_not_unique', '__too_low_aQual', '__not_aligned']]
        counts = df.drop(uncounted.index, inplace=False)

        if save_csv:
            io.save_csv(df=counts, filename=counted_fname)
            io.save_csv(df=uncounted, filename=uncounted_fname)

        fname = counted_fname if save_csv else os.path.join(folder.absolute(), folder.name + file_suffix)
        h = CountFilter((Path(fname), counts))
        if norm_to_rpm:
            h.normalize_to_rpm(uncounted)
        return h


class Pipeline:
    """
    An all-purpose Filter object.


    **Attributes**

    functions: list
        A list of the functions in the Pipeline.
    params: list
        A list of the parameters of the functions in the Pipeline.
    filter_type: Filter object
        The type of Filter objects to which the Pipeline can be applied
    """
    __slots__ = {'functions': 'list of functions to perform', 'params': 'list of function parameters',
                 'filter_type': 'type of filter objects to which the Pipeline can be applied'}

    def __init__(self, filter_type: Union[str,
                                          Type[Filter], Type[DESeqFilter], Type[FoldChangeFilter],
                                          Type[CountFilter]] = Filter):
        """
        :param filter_type: the type of Filter object the Pipeline can be applied to.
        :type filter_type: str or Filter object (default=filtering.Filter)

        :Examples:
            >>> from rnalysis import filtering
            >>> pipe = filtering.Pipeline()
            >>> deseq_pipe = filtering.Pipeline('deseqfilter')

        """
        self.functions = []
        self.params = []

        filter_types = {'filter': Filter, 'deseqfilter': DESeqFilter, 'foldchangefilter': FoldChangeFilter,
                        'countfilter': CountFilter}
        assert isinstance(filter_type,
                          (type, str)), f"'filter_type' must be type of a Filter object, is instead {type(filter_type)}"
        if isinstance(filter_type, str):
            assert filter_type.lower() in filter_types, f"Invalid filter_type {filter_type}. "
            filter_type = filter_types[filter_type.lower()]
        else:
            assert filter_type in filter_types.values(), f"Invalid filter_type {filter_type}"
        self.filter_type = filter_type

    def __str__(self):
        string = f"Pipeline for {self.filter_type.__name__} objects"
        if len(self) > 0:
            string += f":\n\t" + '\n\t'.join(
                self._func_signature(func, params[0], params[1]) for func, params in zip(self.functions, self.params))
        return string

    def __repr__(self):
        string = f"Pipeline('{self.filter_type.__name__}')"
        if len(self) > 0:
            string += ": " + "-->".join(
                self._func_signature(func, params[0], params[1]) for func, params in zip(self.functions, self.params))
        return string

    def __len__(self):
        return len(self.functions)

    @staticmethod
    def _param_string(args: tuple, kwargs: dict):

        """
        Returns a formatted string of the given arguments and keyworded arguments.

        :param args: arguments to format as string
        :type args: tuple
        :param kwargs: keyworded arguments to format as string
        :type kwargs: dict
        :return: a formatted string of arguments and keyworded argumentss
        :rtype: str

        """
        args_str = ', '.join([f"'{arg}'" if isinstance(arg, str) else f"{arg}" for arg in args])
        kwargs_str = ', '.join(
            [f"{key}='{arg}'" if isinstance(arg, str) else f"{key}={arg}" for key, arg in kwargs.items()])
        if len(args_str) == 0:
            return kwargs_str
        elif len(kwargs_str) == 0:
            return args_str
        else:
            return f"{args_str}, {kwargs_str}"

    def _func_signature(self, func: types.FunctionType, args: tuple, kwargs: dict):
        """
        Returns a string functions signature for the given function and arguments.

        :param func: the function or method to generate signature for
        :type func: function
        :param args: arguments given for the function
        :type args: tuple
        :param kwargs: keyworded arguments given for the function
        :type kwargs: dict
        :return: function signature string
        :rtype: str
        """
        return f"{self.filter_type.__name__}.{func.__name__}({self._param_string(args, kwargs)})"

    def add_function(self, func: Union[types.FunctionType, str], *args, **kwargs):

        """
        Add a function to the pipeline. Arguments can be stated with or without the correspoding keywords. \
        For example: Pipeline.add_function('sort', 'columnName', ascending=False, na_position='first'). \
        Pipelines support virtually all functions in the filtering module, including \
        filtering functions, normalization functions, splitting functions and plotting functions. \
        Do not include the 'inplace' argument when adding functions to a pipeline; instead, \
        include it when applying the pipeline to a Filter object using Pipeline.apply_to(). \

        :param func: function to be added
        :type func: function or name of function from the filtering module
        :param args: unkeyworded arguments for the added function in their natural order. For example: 0.1, True
        :param kwargs: keyworded arguments for the added function. For example: opposite=True

        """
        assert isinstance(func, (str, types.FunctionType)), f"'func' must be a function/str, is {type(func)} instead."
        if isinstance(func, str):
            func = func.lower()  # function names are always expected to be lowercase. This prevents capitalized typos.
            assert hasattr(self.filter_type, func), \
                f"Function {func} does not exist for filter_type {self.filter_type}."
            func = getattr(self.filter_type, func)
        else:
            assert hasattr(self.filter_type, func.__name__) and getattr(self.filter_type, func.__name__) == func, \
                f"Function {func.__name__} does not exist for filter_type {self.filter_type}. "
        if 'inplace' in kwargs:
            warnings.warn(
                'The "inplace" argument supplied to this function will be ignored. '
                'To apply the pipeline inplace, state "inplace=True" when calling Pipeline.apply_to().')
        self.functions.append(func)
        self.params.append((args, kwargs))
        print(
            f"Added function '{self._func_signature(func, args, kwargs)}' to the pipeline.")

    def _apply_filter_norm_sort(self, func: types.FunctionType,
                                filter_object: Union[Type[Filter], Type[DESeqFilter], Type[FoldChangeFilter],
                                                     Type[CountFilter]], args: tuple, kwargs: dict, inplace: bool):
        """
        Apply a filtering/normalizing/sorting function.

        :param func: function to apply
        :type func: function
        :param filter_object: Filter object to apply function to
        :type filter_object: Filter, CountFilter, DESeqFilter, or FoldChangeFilter.
        :param args: arguments for the function
        :type args: tuple
        :param kwargs: keyworded arguments for the function
        :type kwargs: dict
        :param inplace: if True, function will be applied inplace.
        :type inplace: bool
        :return: Filter object to which the function was applied.
        """
        kwargs = kwargs.copy()
        if not inplace:
            kwargs['inplace'] = False
            try:
                if isinstance(filter_object, tuple):
                    temp_object = []
                    for obj in filter_object:
                        temp_object.append(func(obj, *args, **kwargs))
                    filter_object = tuple(temp_object)
                else:
                    filter_object = func(filter_object, *args, **kwargs)
            except (ValueError, AssertionError, TypeError) as e:
                raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")
        else:
            kwargs['inplace'] = True
            try:
                if isinstance(filter_object, tuple):
                    for obj in filter_object:
                        func(obj, *args, **kwargs)
                else:
                    func(filter_object, *args, **kwargs)
            except (ValueError, AssertionError, TypeError) as e:
                raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")
        return filter_object

    def _apply_split(self, func: types.FunctionType,
                     filter_object: Union[Type[Filter], Type[DESeqFilter], Type[FoldChangeFilter], Type[CountFilter]],
                     args: tuple, kwargs: dict, other_outputs, other_cnt):
        """
         Apply a splitting function.

        :param func: function to apply
        :type func: function
        :param filter_object: Filter object to apply function to
        :type filter_object: Filter, CountFilter, DESeqFilter, or FoldChangeFilter.
        :param args: arguments for the function
        :type args: tuple
        :param kwargs: keyworded arguments for the function
        :type kwargs: dict
        :param other_outputs: dictionary with additional function outputs
        :type other_outputs: dict
        :param other_cnt: counter for how many times each function was already called
        :type other_cnt: dict
        :return: Filter object to which the function was applied.
        """
        try:
            if isinstance(filter_object, tuple):
                temp_object = []
                for obj in filter_object:
                    temp_res = func(obj, *args, **kwargs)
                    if isinstance(temp_res, tuple):
                        temp_object.extend(temp_res)
                    elif isinstance(temp_res, list):
                        for item in temp_res:
                            if isinstance(item, tuple):
                                temp_object.extend(item)
                            elif item is not None:
                                other_outputs[f"{func.__name__}_{other_cnt.setdefault(func.__name__, 1)}"] = item
                                other_cnt[func.__name__] += 1

                    else:
                        raise ValueError(f"Unrecognized output type {type(temp_res)} from function {func}:\n{temp_res}")
                filter_object = tuple(temp_object)
            else:
                temp_res = func(filter_object, *args, **kwargs)
                if isinstance(temp_res, tuple):
                    filter_object = temp_res
                elif isinstance(temp_res, list):
                    temp_object = []
                    for item in temp_res:
                        if isinstance(item, tuple):
                            temp_object.extend(item)
                        elif item is not None:
                            other_outputs[f"{func.__name__}_{other_cnt.setdefault(func.__name__, 1)}"] = item
                            other_cnt[func.__name__] += 1
                    filter_object = tuple(temp_object)
        except (ValueError, AssertionError, TypeError) as e:
            raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")
        return filter_object

    def _apply_other(self, func, filter_object, args, kwargs, other_outputs, other_cnt):
        """
        Apply a non filtering/splitting/normalizing/sorting function.

        :param func: function to apply
        :type func: function
        :param filter_object: Filter object to apply function to
        :type filter_object: Filter, CountFilter, DESeqFilter, or FoldChangeFilter.
        :param args: arguments for the function
        :type args: tuple
        :param kwargs: keyworded arguments for the function
        :type kwargs: dict
        :param other_outputs: dictionary with additional function outputs
        :type other_outputs: dict
        :param other_cnt: counter for how many times each function was already called
        :type other_cnt: dict
        :return: Filter object to which the function was applied.
        """
        try:
            if isinstance(filter_object, tuple):
                tmp_outputs = []
                for obj in filter_object:
                    this_output = func(obj, *args, **kwargs)
                    if this_output is not None:
                        tmp_outputs.append(this_output)
                if len(tmp_outputs) > 0:
                    other_outputs[f"{func.__name__}_{other_cnt.setdefault(func.__name__, 1)}"] = tuple(tmp_outputs)
                    other_cnt[func.__name__] += 1
            else:
                this_output = func(filter_object, *args, **kwargs)
                if this_output is not None:
                    other_outputs[f"{func.__name__}_{other_cnt.setdefault(func.__name__, 1)}"] = this_output
                    other_cnt[func.__name__] += 1
        except (ValueError, AssertionError, TypeError) as e:
            raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")

    def apply_to(self, filter_object, inplace: bool = True):

        """
        Sequentially apply all functions in the Pipeline to a given Filter object.

        :param filter_object: filter object to apply the Pipeline to. \
        Type of filter_object must be identical to `Pipeline.filter_type`.
        :type filter_object: Filter, CountFilter, DESeqFilter, or FoldChangeFilter
        :param inplace: Determines whether to apply operations in-place or not.
        :type inplace: bool (default=True)
        :return: If inplace=False, a Filter object/tuple of Filter objects will be returned. \
        If the functions in the Pipeline return any additional outputs, \
        they will also be returned in a dictionary. Otherwise, nothing will be returned.
        :rtype: Filter object, Tuple[Filter, dict], dict, or None

        """
        # noinspection PyTypeHints
        assert issubclass(filter_object.__class__,
                          self.filter_type), f"Supplied filter object of type {type(filter_object)} " \
                                             f"mismatches the specified filter_type {self.filter_type}. "
        assert len(self.functions) > 0 and len(self.params) > 0, "Cannot apply an empty pipeline!"

        other_outputs = dict()
        other_cnt = dict()
        # iterate over all functions and arguments
        for func, (args, kwargs) in zip(self.functions, self.params):
            if 'filter' in func.__name__ or func.__name__.startswith('normalize') or func.__name__ == 'sort':
                filter_object = self._apply_filter_norm_sort(func, filter_object, args, kwargs, inplace)
            elif func.__name__.startswith('split'):
                assert not inplace, f"Cannot apply the split function {self._func_signature(func, args, kwargs)} " \
                                    f"when inplace={inplace}!"
                filter_object = self._apply_split(func, filter_object, args, kwargs, other_outputs, other_cnt)
            else:
                self._apply_other(func, filter_object, args, kwargs, other_outputs, other_cnt)

        if not inplace or isinstance(filter_object, tuple):
            if len(other_outputs) > 0:
                return filter_object, other_outputs
            return filter_object
        if len(other_outputs) > 0:
            return other_outputs

    def remove_last_function(self):

        """
        Removes from the Pipeline the last function that was added to it. Removal is in-place.

        """
        assert len(self.functions) > 0 and len(self.params) > 0, "Pipeline is empty, no functions to remove!"
        func = self.functions.pop(-1)
        args, kwargs = self.params.pop(-1)
        print(
            f"Removed function {func.__name__} with parameters [{self._param_string(args, kwargs)}] from the pipeline.")

# TODO: a function that receives a dataframe, and can plot correlation with the ref. table instead of just enrichment
# TODO: add option for mask in clustergram
# TODO: heat map plot of multiple DESEQ files
