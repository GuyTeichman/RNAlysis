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
import copy
import os
import re
import types
from pathlib import Path
from typing import Any, Iterable, List, Tuple, Union, Callable, Sequence, Literal

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import polars.selectors as cs
import seaborn as sns
from grid_strategy import strategies
from scipy.stats import gstd, spearmanr, sem
from scipy.stats.mstats import gmean
from sklearn.decomposition import PCA
from sklearn.preprocessing import PowerTransformer, StandardScaler
from tqdm.auto import tqdm

from rnalysis.utils import clustering, generic, ontology, settings, validation, differential_expression, \
    param_typing, genome_annotation, parsing, io
from rnalysis.utils.generic import readable_name
from rnalysis.utils.param_typing import *
import warnings

@readable_name('Generic table')
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
    __slots__ = {'fname': 'filename with full path', 'df': 'pandas.DataFrame with the data'}

    def __init__(self, fname: Union[str, Path], drop_columns: Union[str, List[str]] = None,
                 suppress_warnings: bool = False):

        """
        Load a table.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded table.
        :type drop_columns: str, list of str, or None (default=None)

        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/counted.csv")

        """
        # init from a file (load the csv into a DataFrame/Series)
        assert isinstance(fname, (str, Path))
        self.fname = Path(fname)
        self.df = io.load_table(fname, squeeze=True, drop_columns=drop_columns)
        if isinstance(self.df, pl.Series):
            self.df = self.df.to_frame()
        # check for duplicate indices
        if not suppress_warnings:
            self._init_warnings()

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame, name: Union[str, Path]) -> 'Filter':
        obj = cls.__new__(cls)
        obj.df = df.clone()
        obj.fname = Path(name)
        obj._init_warnings()
        return obj

    def _init_warnings(self):
        if isinstance(self.df, pl.DataFrame):
            cond = self.df.select(pl.first()).is_duplicated().any()
        elif isinstance(self.df, pl.Series):
            cond = self.df.is_duplicated().any()
        else:
            raise TypeError("The 'df' attribute must be a polars DataFrame or Series object.")
        if cond:
            warnings.warn("This Filter object contains multiple rows with the same name/index.")

    def __repr__(self):
        return f"{type(self).__name__}('{self.fname.as_posix()}')"

    def __str__(self):
        return f"{self.readable_name} named '{self.fname.stem}{self.fname.suffix}'"

    def __len__(self):
        """
        Returns the number of rows in the Filter object

        """
        return self.shape[0]

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        if self.df.equals(other.df) and self.shape == other.shape:
            return True
        return False

    def __contains__(self, item):
        return self.df.filter(pl.first() == item).shape[0] > 0

    def __iter__(self):
        yield from self.df.select(pl.first())

    def __copy__(self):
        return type(self).from_dataframe(self.df.clone(), self.fname)

    @property
    def columns(self) -> list:
        """
        The columns of df.

        :return: a list of the columns in the Filter object.
        :rtype: list
        """
        return list(self.df.columns[1:])

    @property
    def shape(self) -> Tuple[int, int]:
        return self.df.shape[0], self.df.shape[1] - 1

    def _update(self, **kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self, key, val)
            except AttributeError:
                raise AttributeError(f"Cannot update attribute {key} for {type(self)} object: attribute does not exist")

    def _inplace(self, new_df: pl.DataFrame, opposite: bool, inplace: bool, suffix: str,
                 printout_operation: str = 'filter', **filter_update_kwargs):

        """
        Executes the user's choice whether to filter in-place or create a new instance of the Filter object.

        :param new_df: the post-filtering DataFrame
        :type new_df: pl.DataFrame
        :param opposite: Determines whether to return the filtration ,or its opposite.
        :type opposite: bool
        :param inplace: Determines whether to filter in-place or not.
        :type inplace: bool
        :param suffix: The suffix to be added to the filename
        :type suffix: str
        :return: If inplace is False, returns a new instance of the Filter object.

        """
        legal_operations = {'filter': 'Filtering', 'normalize': 'Normalization', 'sort': 'Sorting',
                            'transform': 'Transformation', 'translate': 'Translation'}
        assert isinstance(inplace, bool), "'inplace' must be True or False!"
        assert isinstance(opposite, bool), "'opposite' must be True or False!"
        assert printout_operation.lower() in legal_operations, \
            f"Invalid input for variable 'printout_operation': {printout_operation}"
        # when user requests the opposite of a filter, return the Set Difference between the filtering result and self
        if opposite:
            new_df = self.df.filter(~pl.first().is_in(new_df.select(pl.first())))
            suffix += 'opposite'

        # update filename with the suffix of the operation that was just performed
        new_fname = Path(os.path.join(str(self.fname.parent), f"{self.fname.stem}{suffix}{self.fname.suffix}"))

        # generate printout for user ("Filtered X features, leaving Y... filtered inplace/not inplace")
        printout = ''
        if printout_operation.lower() == 'filter':
            printout += f"Filtered {self.shape[0] - new_df.shape[0]} features, leaving {new_df.shape[0]} " \
                        f"of the original {self.shape[0]} features. "
        elif printout_operation.lower() == 'translate':
            printout += f"Translated the gene IDs of {new_df.shape[0]} features. "
            if self.shape[0] != new_df.shape[0]:
                printout += f"Filtered {self.shape[0] - new_df.shape[0]} unmapped features, " \
                            f"leaving {new_df.shape[0]} of the original {self.shape[0]} features. "

        else:
            printout += f'{printout_operation.capitalize().rstrip("e")}ed {new_df.shape[0]} features. '
        # if inplace, modify the df, fname and other attributes of self
        if inplace:
            printout += f'{printout_operation.capitalize().rstrip("e")}ed inplace.'
            print(printout)
            self._update(df=new_df, fname=new_fname, **filter_update_kwargs)
        # if not inplace, copy self, modify the df/fname properties of the copy, and return it
        else:
            printout += f'{legal_operations[printout_operation]} result saved to new object.'
            print(printout)
            new_obj = self.__copy__()
            new_obj._update(df=new_df, fname=new_fname, **filter_update_kwargs)
            return new_obj

    def save_table(self, suffix: Literal['.csv', '.tsv', '.parquet'] = '.csv',
                   alt_filename: Union[None, str, Path] = None):

        """
        Save the current filtered data table.

        :param suffix: the file suffix
        :type suffix: '.csv', '.tsv', or '.parquet' (default='.csv')
        :param alt_filename: If None, file name will be generated automatically \
        according to the filtering methods used. \
        If it's a string, it will be used as the name of the saved file. Example input: 'myfilename'
        :type alt_filename: str, pathlib.Path, or None (default)

        """
        # save with the default filename if no alternative filename was given
        if alt_filename is None:
            alt_filename = self.fname.with_suffix(suffix)
        else:
            assert isinstance(alt_filename, (str, Path)), \
                f"'alt_filename' must be a string or Path object. Instead got {type(alt_filename)}."
            alt_filename = self.fname.parent.joinpath(alt_filename).with_suffix(suffix)
        io.save_table(self.df, alt_filename)

    def save_csv(self, alt_filename: Union[None, str, Path] = None):

        """
        Saves the current filtered data to a .csv file.

        :param alt_filename: If None, file name will be generated automatically \
        according to the filtering methods used. \
        If it's a string, it will be used as the name of the saved file. Example input: 'myfilename'
        :type alt_filename: str, pathlib.Path, or None (default)

        """
        self.save_table('.csv', alt_filename)

    def save_parquet(self, alt_filename: Union[None, str, Path] = None):

        """
        Saves the current filtered data to a .parquet file.

        :param alt_filename: If None, file name will be generated automatically \
        according to the filtering methods used. \
        If it's a string, it will be used as the name of the saved file. Example input: 'myfilename'
        :type alt_filename: str, pathlib.Path, or None (default)

        """
        self.save_table('.parquet', alt_filename)

    @staticmethod
    def _from_string(msg: str = '', delimiter: str = '\n'):

        """
        Takes a manual string input from the user, and then splits it using a delimiter into a list of values.

        :param msg: a promprt to be printed to the user
        :param delimiter: the delimiter used to separate the values. Default is '\n'
        :return: A list of the comma-seperated values the user inserted.
        """
        string = input(msg)
        split = string.split(sep=delimiter)
        if split[-1] == '':
            split = split[:-1]
        return split

    @readable_name('Table head')
    def head(self, n: PositiveInt = 5) -> pl.DataFrame:

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

    @readable_name('Table tail')
    def tail(self, n: PositiveInt = 5) -> pl.DataFrame:

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

    @readable_name('Concatenate two tables (based on column names)')
    def concatenate(self, other: Sequence[str]) -> 'Filter':
        assert isinstance(other, Filter), f"Expected a Filter object, got {type(other)} instead."
        assert sorted(self.columns) == sorted(other.columns), "The two tables do not have the same columns!"
        assert len({item[0] for item in self.df.select(pl.first()).iter_rows()}.intersection(
            {item[0] for item in
             other.df.select(pl.first()).iter_rows()})) == 0, "The two tables have overlapping indices!"
        new_df = pl.concat([self.df, other.df], how='vertical')
        return self.from_dataframe(new_df, Path(self.fname.stem + '_' + other.fname.stem + '.csv'))

    @readable_name('Filter rows with duplicate names/IDs')
    def filter_duplicate_ids(self, keep: Literal['first', 'last', 'neither'] = 'first', opposite: bool = False,
                             inplace: bool = True):
        """
        Filter out rows with duplicate names/IDs (index).

        :param keep: determines which of the duplicates to keep for each group of duplicates. 'first' will keep the \
        first duplicate found for each group; 'last' will keep the last; \
        and 'neither' will remove *all* of the values in the group.
        :type keep: 'first', 'last', or 'neither' (default='first')
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and filtered instance of the Filter object.
        """
        suffix = f'_dropduplicateskeep{keep}'
        if keep == 'neither':
            keep = 'none'
        new_df = self.df.unique(cs.first(), keep=keep)
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter specific rows by name')
    def filter_by_row_name(self, row_names: Union[str, List[str]], opposite: bool = False, inplace: bool = True):
        """
        Filter out specific rows from the table by their name (index).

        :param row_names: list of row names to be removed from the table.
        :type row_names: str or list of str
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and filtered instance of the Filter object.
        """
        suffix = '_filterbyrowname'
        row_names = parsing.data_to_list(row_names)
        new_df = self.df.filter(~pl.first().is_in(row_names))
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Drop columns from the table')
    def drop_columns(self, columns: param_typing.ColumnNames, inplace: bool = True):
        """
        Drop specific columns from the table.

        :param columns: The names of the column/columns to be dropped fro mthe table.
        :type columns: str or list of str
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and filtered instance of the Filter object.
        """
        suffix = '_dropcolumns'
        columns = parsing.data_to_list(columns)
        for col in columns:
            assert col in self.columns, f"column '{col}' does not exist!"
        new_df = self.df.drop(columns)
        return self._inplace(new_df, False, inplace, suffix, 'transform')

    @readable_name('Histogram of a column')
    def histogram(self, column: param_typing.ColumnName, bins: int = 100,
                  x_label: Union[str, Literal['auto']] = 'auto', y_logscale: bool = False, x_logscale: bool = False):
        assert column in self.columns, f"column '{column}' does not exist!"
        fig, ax = plt.subplots()
        ax.hist(self.df[column], bins=bins)
        ax.set_xlabel(column if x_label == 'auto' else x_label)
        ax.set_ylabel('Frequency')
        ax.set_title(f'Histogram of column "{column}"')
        if y_logscale:
            ax.set_yscale('log')
        if x_logscale:
            ax.set_xscale('log')
        plt.show()
        return fig

    @staticmethod
    def _create_one2many_table(one2many_dict: io.OrthologDict, mode: Literal['ortholog', 'paralog'] = 'ortholog'):
        pairs = []
        for key, vals in one2many_dict.mapping_dict.items():
            for val in vals:
                pairs.append((key, val))
        one2many_table = pl.DataFrame(pairs, schema=['gene', mode])
        return one2many_table

    @readable_name('Find paralogs within species (using PantherDB)')
    def find_paralogs_panther(self, organism: Union[Literal['auto'], str, int, Literal[get_panther_taxons()]] = 'auto',
                              gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto'):
        """
        Find paralogs within the same species using the PantherDB database.

        :param organism: organism name or NCBI taxon ID of the input genes' source species.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :return: DataFrame describing all discovered paralog mappings.
        """
        (taxon_id, organism), gene_id_type = io.get_taxon_and_id_type(organism, gene_id_type, self.index_set)

        mapper = io.PantherOrthologMapper(taxon_id, taxon_id, gene_id_type)
        ids = parsing.data_to_tuple(self.index_set)
        ortho_dict_one2many = mapper.get_paralogs(ids)

        one2many_table = self._create_one2many_table(ortho_dict_one2many, 'paralog')
        return one2many_table

    @readable_name('Find paralogs within species (using Ensembl)')
    def find_paralogs_ensembl(self, organism: Union[Literal['auto'], str, int, Literal[get_ensembl_taxons()]] = 'auto',
                              gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                              filter_percent_identity: bool = True):
        """
        Find paralogs within the same species using the Ensembl database.

        :param organism: organism name or NCBI taxon ID of the input genes' source species.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :param filter_percent_identity: if True (default), when encountering non-unique ortholog mappings, \
        *RNAlysis* will only keep the mappings with the highest percent_identity score.
        :type filter_percent_identity: bool (default=True)
        :return: DataFrame describing all discovered paralog mappings.
        """
        (taxon_id, organism), gene_id_type = io.get_taxon_and_id_type(organism, gene_id_type, self.index_set)

        mapper = io.EnsemblOrthologMapper(taxon_id, taxon_id, gene_id_type)
        ids = parsing.data_to_tuple(self.index_set)
        ortho_dict_one2many = mapper.get_paralogs(ids, filter_percent_identity)
        one2many_table = self._create_one2many_table(ortho_dict_one2many, 'paralog')
        return one2many_table

    @readable_name('Map genes to nearest orthologs (using PantherDB)')
    def map_orthologs_panther(self, map_to_organism: Union[str, int, Literal[get_panther_taxons()]],
                              map_from_organism: Union[
                                  Literal['auto'], str, int, Literal[get_panther_taxons()]] = 'auto',
                              gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                              filter_least_diverged: bool = True,
                              non_unique_mode: Literal[ORTHOLOG_NON_UNIQUE_MODES] = 'first',
                              remove_unmapped_genes: bool = False, inplace: bool = True):
        """
        Map genes to their nearest orthologs in a different species using the PantherDB database. \
        This function generates a table describing all matching discovered ortholog pairs (both unique and non-unique) \
        and returns it, and can also translate the genes in this data table into their nearest ortholog, \
        as well as remove unmapped genes.

        :param map_to_organism: organism name or NCBI taxon ID of the target species for ortholog mapping.
        :type map_to_organism: str or int
        :param map_from_organism: organism name or NCBI taxon ID of the input genes' source species.
        :type map_from_organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :param filter_least_diverged: if True (default), *RNAlysis* will only fetch ortholog mappings that were \
        flagged as a 'least diverged ortholog' on the PantherDB database. \
        You can read more about this flag on the PantherDB website: https://www.pantherdb.org/genes/
        :type filter_least_diverged: bool (default=True)
        :param non_unique_mode: How to handle non-unique mappings. 'first' will keep the first mapping found for each \
        gene; 'last' will keep the last; 'random' will keep a random mapping; \
        and 'none' will discard all non-unique mappings.
        :type non_unique_mode: 'first', 'last', 'random', or 'none' (default='first')
        :param remove_unmapped_genes: if True, rows with gene names/IDs that could not be mapped to an ortholog \
        will be dropped from the table. \
        Otherwise, they will remain in the table with their original gene name/ID.
        :type remove_unmapped_genes: bool (default=False)
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: DataFrame describing all discovered mappings (unique and otherwise). If inplace=True, \
        returns a filtered instance of the Filter object as well.
        """
        taxon_id_to = io.map_taxon_id(map_to_organism)[0]
        (taxon_id_from, organism), gene_id_type = io.get_taxon_and_id_type(map_from_organism, gene_id_type,
                                                                           self.index_set)

        mapper = io.PantherOrthologMapper(taxon_id_to, taxon_id_from, gene_id_type)
        ids = parsing.data_to_tuple(self.index_set)
        ortho_dict_one2one, ortho_dict_one2many = mapper.get_orthologs(ids, non_unique_mode, filter_least_diverged)

        one2many_table = self._create_one2many_table(ortho_dict_one2many)

        new_df = self._apply_dict_map(ortho_dict_one2one, remove_unmapped_genes)
        suffix = f'_orthologsPanther{taxon_id_to}'
        return self._inplace(new_df, False, inplace, suffix, 'translate'), one2many_table

    def _apply_dict_map(self, one2one_map: Union[io.OrthologDict, io.GeneIDDict], remove_unmapped: bool):
        map_df = pl.DataFrame([list(one2one_map.mapping_dict.keys()), list(one2one_map.mapping_dict.values())],
                              schema=[self.df.columns[0], 'mappedVals'])
        new_df = self.df.join(map_df, left_on=self.df.columns[0], right_on=map_df.columns[0], how='left')

        if remove_unmapped:
            new_df = new_df.filter(pl.col('mappedVals').is_not_null())
        else:
            new_df = new_df.with_columns(pl.col('mappedVals').fill_null(pl.col(self.df.columns[0])).alias('mappedVals'))

        new_df = new_df.with_columns(pl.col('mappedVals').alias(self.df.columns[0])).drop('mappedVals')
        return new_df

    @readable_name('Map genes to nearest orthologs (using Ensembl)')
    def map_orthologs_ensembl(self, map_to_organism: Union[str, int, Literal[get_ensembl_taxons()]],
                              map_from_organism: Union[
                                  Literal['auto'], str, int, Literal[get_ensembl_taxons()]] = 'auto',
                              gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                              filter_percent_identity: bool = True,
                              non_unique_mode: Literal[ORTHOLOG_NON_UNIQUE_MODES] = 'first',
                              remove_unmapped_genes: bool = False, inplace: bool = True):
        """
        Map genes to their nearest orthologs in a different species using the Ensembl database. \
        This function generates a table describing all matching discovered ortholog pairs (both unique and non-unique) \
        and returns it, and can also translate the genes in this data table into their nearest ortholog, \
        as well as remove unmapped genes.

        :param map_to_organism: organism name or NCBI taxon ID of the target species for ortholog mapping.
        :type map_to_organism: str or int
        :param map_from_organism: organism name or NCBI taxon ID of the input genes' source species.
        :type map_from_organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :param filter_percent_identity: if True (default), when encountering non-unique ortholog mappings, \
        *RNAlysis* will only keep the mappings with the highest percent_identity score.
        :type filter_percent_identity: bool (default=True)
        :param non_unique_mode: How to handle non-unique mappings. 'first' will keep the first mapping found for each \
        gene; 'last' will keep the last; 'random' will keep a random mapping; \
        and 'none' will discard all non-unique mappings.
        :type non_unique_mode: 'first', 'last', 'random', or 'none' (default='first')
        :param remove_unmapped_genes: if True, rows with gene names/IDs that could not be mapped to an ortholog \
        will be dropped from the table. \
        Otherwise, they will remain in the table with their original gene name/ID.
        :type remove_unmapped_genes: bool (default=False)
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: DataFrame describing all discovered mappings (unique and otherwise). If inplace=True, \
        returns a filtered instance of the Filter object as well.
        """
        taxon_id_to = io.map_taxon_id(map_to_organism)[0]
        (taxon_id_from, organism), gene_id_type = io.get_taxon_and_id_type(map_from_organism, gene_id_type,
                                                                           self.index_set)

        mapper = io.EnsemblOrthologMapper(taxon_id_to, taxon_id_from, gene_id_type)
        ids = parsing.data_to_tuple(self.index_set)
        ortho_dict_one2one, ortho_dict_one2many = mapper.get_orthologs(ids, non_unique_mode, filter_percent_identity)
        one2many_table = self._create_one2many_table(ortho_dict_one2many)

        new_df = self._apply_dict_map(ortho_dict_one2one, remove_unmapped_genes)
        suffix = f'_orthologsEnsembl{taxon_id_to}'
        return self._inplace(new_df, False, inplace, suffix, 'translate'), one2many_table

    @readable_name('Map genes to nearest orthologs (using PhylomeDB)')
    def map_orthologs_phylomedb(self, map_to_organism: Union[str, int, Literal[get_phylomedb_taxons()]],
                                map_from_organism: Union[
                                    Literal['auto'], str, int, Literal[get_phylomedb_taxons()]] = 'auto',
                                gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                consistency_score_threshold: Fraction = 0.5, filter_consistency_score: bool = True,
                                non_unique_mode: Literal[ORTHOLOG_NON_UNIQUE_MODES] = 'first',
                                remove_unmapped_genes: bool = False, inplace: bool = True):
        """
        Map genes to their nearest orthologs in a different species using the PhylomeDB database. \
        This function generates a table describing all matching discovered ortholog pairs (both unique and non-unique) \
        and returns it, and can also translate the genes in this data table into their nearest ortholog, \
        as well as remove unmapped genes.

        :param map_to_organism: organism name or NCBI taxon ID of the target species for ortholog mapping.
        :type map_to_organism: str or int
        :param map_from_organism: organism name or NCBI taxon ID of the input genes' source species.
        :type map_from_organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
    `   :param consistency_score_threshold: the minimum consistency score required for an ortholog mapping to be \
        considered valid. Consistency scores are calculated by PhylomeDB and represent the confidence of the \
        ortholog mapping. setting consistency_score_threshold to 0 will keep all mappings. You can read more about \
        PhylomeDB consistency score on the PhylomeDB website: orthology.phylomedb.org/help
        :type consistency_score_threshold: float between 0 and 1 (default=0.5)
        :param filter_consistency_score: if True (default), when encountering non-unique ortholog mappings, \
        *RNAlysis* will only keep the mappings with the highest consistency score.
        :type filter_consistency_score: bool (default=True)
        :param non_unique_mode: How to handle non-unique mappings. 'first' will keep the first mapping found for each \
        gene; 'last' will keep the last; 'random' will keep a random mapping; \
        and 'none' will discard all non-unique mappings.
        :type non_unique_mode: 'first', 'last', 'random', or 'none' (default='first')
        :param remove_unmapped_genes: if True, rows with gene names/IDs that could not be mapped to an ortholog \
        will be dropped from the table. \
        Otherwise, they will remain in the table with their original gene name/ID.
        :type remove_unmapped_genes: bool (default=False)
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: DataFrame describing all discovered mappings (unique and otherwise). If inplace=True, \
        returns a filtered instance of the Filter object as well.
        """
        taxon_id_to = io.map_taxon_id(map_to_organism)[0]
        (taxon_id_from, organism), gene_id_type = io.get_taxon_and_id_type(map_from_organism, gene_id_type,
                                                                           self.index_set)

        mapper = io.PhylomeDBOrthologMapper(taxon_id_to, taxon_id_from, gene_id_type)
        ids = parsing.data_to_tuple(self.index_set)
        ortho_dict_one2one, ortho_dict_one2many = mapper.get_orthologs(ids, non_unique_mode,
                                                                       consistency_score_threshold,
                                                                       filter_consistency_score)
        one2many_table = self._create_one2many_table(ortho_dict_one2many)
        new_df = self._apply_dict_map(ortho_dict_one2one, remove_unmapped_genes)

        suffix = f'_orthologsPhylomeDB{taxon_id_to}'
        return self._inplace(new_df, False, inplace, suffix, 'translate'), one2many_table

    @readable_name('Map genes to nearest orthologs (using OrthoInspector)')
    def map_orthologs_orthoinspector(self, map_to_organism: Union[str, int, Literal[get_panther_taxons()]],
                                     map_from_organism: Union[
                                         Literal['auto'], str, int, Literal[get_panther_taxons()]] = 'auto',
                                     gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                     non_unique_mode: Literal[ORTHOLOG_NON_UNIQUE_MODES] = 'first',
                                     remove_unmapped_genes: bool = False, inplace: bool = True):
        """
        Map genes to their nearest orthologs in a different species using the OrthoInspector database. \
        This function generates a table describing all matching discovered ortholog pairs (both unique and non-unique) \
        and returns it, and can also translate the genes in this data table into their nearest ortholog, \
        as well as remove unmapped genes.

        :param map_to_organism: organism name or NCBI taxon ID of the target species for ortholog mapping.
        :type map_to_organism: str or int
        :param map_from_organism: organism name or NCBI taxon ID of the input genes' source species.
        :type map_from_organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :param non_unique_mode: How to handle non-unique mappings. 'first' will keep the first mapping found for each \
        gene; 'last' will keep the last; 'random' will keep a random mapping; \
        and 'none' will discard all non-unique mappings.
        :type non_unique_mode: 'first', 'last', 'random', or 'none' (default='first')
        :param remove_unmapped_genes: if True, rows with gene names/IDs that could not be mapped to an ortholog \
        will be dropped from the table. \
        Otherwise, they will remain in the table with their original gene name/ID.
        :type remove_unmapped_genes: bool (default=False)
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: DataFrame describing all discovered mappings (unique and otherwise). If inplace=True, \
        returns a filtered instance of the Filter object as well.
        """
        taxon_id_to = io.map_taxon_id(map_to_organism)[0]
        (taxon_id_from, organism), gene_id_type = io.get_taxon_and_id_type(map_from_organism, gene_id_type,
                                                                           self.index_set)

        mapper = io.OrthoInspectorOrthologMapper(taxon_id_to, taxon_id_from, gene_id_type)
        ids = parsing.data_to_tuple(self.index_set)
        ortho_dict_one2one, ortho_dict_one2many = mapper.get_orthologs(ids, non_unique_mode)
        one2many_table = self._create_one2many_table(ortho_dict_one2many)

        new_df = self._apply_dict_map(ortho_dict_one2one, remove_unmapped_genes)
        suffix = f'_orthologsOrthoInspector{taxon_id_to}'
        return self._inplace(new_df, False, inplace, suffix, 'translate'), one2many_table

    @readable_name('Translate gene IDs')
    def translate_gene_ids(self, translate_to: Union[str, Literal[get_gene_id_types()]],
                           translate_from: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                           remove_unmapped_genes: bool = False, inplace: bool = True):
        """
        Translates gene names/IDs from one type to another. \
        Mapping is done using the UniProtKB Gene ID Mapping service. \
        You can choose to optionally drop from the table all rows that failed to be translated.

        :param translate_to: the gene ID type to translate gene names/IDs to. \
        For example: UniProtKB, Ensembl, Wormbase.
        :type translate_to: str
        :param translate_from: the gene ID type to translate gene names/IDs from. \
        For example: UniProtKB, Ensembl, Wormbase. If translate_from='auto', \
        *RNAlysis* will attempt to automatically determine the gene ID type of the features in the table.
        :type translate_from: str or 'auto' (default='auto')
        :param remove_unmapped_genes: if True, rows with gene names/IDs that could not be translated \
        will be dropped from the table. \
        Otherwise, they will remain in the table with their original gene name/ID.
        :type remove_unmapped_genes: bool (default=False)
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and translated instance of the Filter object.
        """
        gene_ids = parsing.data_to_tuple(self.df.select(pl.first()).to_series().to_list())
        if translate_from.lower() == 'auto':
            translator, translate_from, _ = io.find_best_gene_mapping(gene_ids, None, (translate_to,))
        else:
            translator = io.GeneIDTranslator(translate_from, translate_to).run(gene_ids)
        new_df = self._apply_dict_map(translator, remove_unmapped_genes)
        suffix = f'_translateFrom{translate_from.replace(" ", "").replace("/", "")}' \
                 f'to{translate_to.replace(" ", "").replace("/", "")}'
        return self._inplace(new_df, False, inplace, suffix, 'translate')

    @readable_name('Filter by percentile')
    def filter_percentile(self, percentile: param_typing.Fraction, column: param_typing.ColumnName,
                          interpolate: param_typing.QUANTILE_INTERPOLATION_METHODS = 'linear', opposite: bool = False,
                          inplace: bool = True):

        """
        Removes all entries above the specified percentile in the specified column. \
        For example, if the column were 'pvalue' and the percentile was 0.5, then all features whose pvalue is above \
        the median pvalue will be filtered out.

        :type percentile: float between 0 and 1
        :param percentile: The percentile that all features above it will be filtered out.
        :type column: str
        :param column: Name of the DataFrame column according to which the filtering will be performed.
        :param interpolate: interpolation method to use when the desired quantile lies between two data points.
        :type interpolate: 'nearest', 'higher', 'lower', 'midpoint' or 'linear' (default='linear')
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
        assert isinstance(column, str) and column in self.columns, "Invalid column name!"
        suffix = f'_below{percentile}percentile'
        new_df = self.df.filter(pl.col(column) <= pl.quantile(column, percentile, interpolate))
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Split by percentile')
    def split_by_percentile(self, percentile: param_typing.Fraction, column: param_typing.ColumnName,
                            interpolate: param_typing.QUANTILE_INTERPOLATION_METHODS = 'linear') -> tuple:

        """
        Splits the features in the Filter object into two non-overlapping Filter objects: \
        one containing features below the specified percentile in the specfieid column, \
        and the other containing features about the specified percentile in the specified column.

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
        return self.filter_percentile(percentile, column, interpolate, opposite=False,
                                      inplace=False), self.filter_percentile(percentile, column, interpolate,
                                                                             opposite=True, inplace=False)

    def _get_ref_srs_from_gtf(self, gtf_path: Union[str, Path],
                              attribute_name: Union[Literal[BIOTYPE_ATTRIBUTE_NAMES], str] = 'gene_biotype',
                              feature_type: Literal['gene', 'transcript'] = 'gene'):
        with tqdm(desc='Parsing GTF file', total=8) as pbar:
            for use_name in [False, True]:
                for use_version in [True, False]:
                    for split_ids in [True, False]:
                        try:
                            mapping = genome_annotation.map_gene_to_attr(gtf_path, attribute_name, feature_type,
                                                                         use_name, use_version, split_ids)
                        except KeyError:
                            continue
                        pbar.update(1)
                        if len(mapping) == 0:
                            continue
                        ref_df = pl.DataFrame([list(mapping.keys()), list(mapping.values())],
                                              schema=[feature_type, 'biotype'])
                        if len(parsing.data_to_set(ref_df.select(pl.first())).intersection(
                            parsing.data_to_set(self.df.select(pl.first())))) == 0:
                            continue

                        pbar.update(8)
                        return ref_df.drop_nulls()
            raise ValueError("Failed to map features to biotypes with the given GTF file and parameters. "
                             "Please try again with different parameters or a different GTF file. ")

    @readable_name('Filter by feature biotype (based on a GTF file)')
    def filter_biotype_from_gtf(self, gtf_path: Union[str, Path],
                                biotype: Union[Literal[BIOTYPES], str, List[str]] = 'protein_coding',
                                attribute_name: Union[Literal[BIOTYPE_ATTRIBUTE_NAMES], str] = 'gene_biotype',
                                feature_type: Literal['gene', 'transcript'] = 'gene',
                                opposite: bool = False, inplace: bool = True):
        """
        Filters out all features that do not match the indicated biotype/biotypes \
        (for example: 'protein_coding', 'ncRNA', etc). \
        The data about feature biotypes is drawn from a GTF (Gene transfer format) file supplied by the user.

        :param gtf_path: Path to your GTF (Gene transfer format) file. The file should match the type of \
        gene names/IDs you use in your table, and should contain an attribute describing biotype.
        :type gtf_path: str or Path
        :param biotype: the biotypes which will not be filtered out.
        :type biotype: str or list of strings
        :param attribute_name: name of the attribute in your GTF file that describes feature biotype.
        :type attribute_name: str (default='gene_biotype')
        :param feature_type: determined whether the features/rows in your data table describe \
        individual genes or transcripts.
        :type feature_type: 'gene' or 'transcript' (default='gene')
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        """
        biotype = parsing.data_to_set(biotype)
        # make sure 'biotype' is a list of strings
        assert validation.isinstanceiter(biotype, str), "biotype must be a string or a list of strings!"
        assert Path(gtf_path).exists(), "the given gtf path does not exist!"
        suffix = f"_{'_'.join(sorted(biotype))}"

        ref_srs = self._get_ref_srs_from_gtf(gtf_path, attribute_name, feature_type)
        # gene names which remain after filtering are True in the mask AND were previously in the df
        gene_names = parsing.data_to_set(ref_srs.filter(pl.last().is_in(biotype))).intersection(
            parsing.data_to_set(self.df.select(pl.first())))
        new_df = self.df.filter(pl.first().is_in(gene_names))
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by feature biotype (based on a reference table)')
    def filter_biotype_from_ref_table(self, biotype: Union[Literal[BIOTYPES], str, List[str]] = 'protein_coding',
                                      ref: Union[str, Path, Literal['predefined']] = 'predefined',
                                      opposite: bool = False, inplace: bool = True):

        """
        Filters out all features that do not match the indicated biotype/biotypes \
        (for example: 'protein_coding', 'ncRNA', etc). \
        The data about feature biotypes is drawn from a Biotype Reference Table supplied by the user.

        :type biotype: string or list of strings
        :param biotype: the biotypes which will not be filtered out.
        :param ref: Name of the biotype reference file used to determine biotypes. \
        Default is the path defined by the user in the settings.yaml file.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.

        :Examples:
            >>> from rnalysis import filtering
            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows whose biotype is 'protein_coding'
            >>> counts.filter_biotype_from_ref_table('protein_coding',ref='tests/biotype_ref_table_for_tests.csv')
            Filtered 9 features, leaving 13 of the original 22 features. Filtered inplace.

            >>> counts = filtering.Filter('tests/test_files/counted.csv')
            >>> # keep only rows whose biotype is 'protein_coding' or 'pseudogene'
            >>> counts.filter_biotype_from_ref_table(['protein_coding','pseudogene'],ref='tests/biotype_ref_table_for_tests.csv')
            Filtered 0 features, leaving 22 of the original 22 features. Filtered inplace.

        """
        biotype = parsing.data_to_set(biotype)
        # make sure 'biotype' is a list of strings
        assert validation.isinstanceiter(biotype, str), "biotype must be a string or a list of strings!"
        suffix = f"_{'_'.join(sorted(biotype))}"
        # load the biotype reference table
        ref = settings.get_biotype_ref_path(ref)
        ref_df = io.load_table(ref).drop_nulls()
        validation.validate_biotype_table(ref_df)
        # generate set of legal inputs
        legal_inputs = set(ref_df['biotype'].unique())
        for bio in biotype:
            assert bio in legal_inputs, f"biotype {bio} is not a legal string!"
        gene_names = parsing.data_to_set(ref_df.filter(pl.last().is_in(biotype))).intersection(
            parsing.data_to_set(self.df.select(pl.first())))
        new_df = self.df.filter(pl.first().is_in(gene_names))
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by Gene Ontology (GO) annotation')
    def filter_by_go_annotations(self, go_ids: Union[str, List[str]], mode: Literal['union', 'intersection'] = 'union',
                                 organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                                 gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                 propagate_annotations: bool = True,
                                 evidence_types: Union[Literal[('any',) + GO_EVIDENCE_TYPES], Iterable[
                                     Literal[GO_EVIDENCE_TYPES]]] = 'any',
                                 excluded_evidence_types: Union[Literal[GO_EVIDENCE_TYPES], Iterable[Literal[
                                     GO_EVIDENCE_TYPES]]] = (),
                                 databases: Union[str, Iterable[str]] = 'any',
                                 excluded_databases: Union[str, Iterable[str]] = (),
                                 qualifiers: Union[
                                     Literal[('any',) + GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'any',
                                 excluded_qualifiers: Union[
                                     Literal[GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'not',
                                 opposite: bool = False, inplace: bool = True):
        """
        Filters genes according to GO annotations, keeping only genes that are annotated with a specific GO term. \
        When multiple GO terms are given, filtering can be done in 'union' mode \
        (where genes that belong to at least one GO term are not filtered out), or in 'intersection' mode \
        (where only genes that belong to ALL GO terms are not filtered out).

        :param go_ids:
        :type go_ids: str or list of str
        :type mode: 'union' or 'intersection'.
        :param mode: If 'union', filters out every genomic feature that does not belong to one or more \
        of the indicated attributes. If 'intersection', \
        filters out every genomic feature that does not belong to ALL of the indicated attributes.
        param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :param propagate_annotations: determines the propagation method of GO annotations. \
        'no' does not propagate annotations at all; 'classic' propagates all annotations up to the DAG tree's root; \
        'elim' terminates propagation at nodes which show significant enrichment; 'weight' performs propagation in a \
        weighted manner based on the significance of children nodes relatively to their parents; and 'allm' uses a \
        combination of all proopagation methods. To read more about the propagation methods, \
        see Alexa et al: https://pubmed.ncbi.nlm.nih.gov/16606683/
        :type propagate_annotations: 'classic', 'elim', 'weight', 'all.m', or 'no' (default='elim')
        :param evidence_types: only annotations with the specified evidence types will be included in the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', 'author', \
        'curator', 'electronic', or 'any' (default='any')
        :param excluded_evidence_types: annotations with the specified evidence types will be \
        excluded from the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type excluded_evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', \
        'author', 'curator', 'electronic', or None (default=None)
        :param databases: only annotations from the specified databases will be included in the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type databases: str, Iterable of str, or 'any' (default)
        :param excluded_databases: annotations from the specified databases will be excluded from the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type excluded_databases: str, Iterable of str, or None (default)
        :param qualifiers: only annotations with the speficied qualifiers will be included in the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type qualifiers: str, Iterable of str, or 'any' (default)
        :param excluded_qualifiers: annotations with the speficied qualifiers will be excluded from the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type excluded_qualifiers: str, Iterable of str, or None (default='not')
        :type opposite: bool (default=False)
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current FeatureSet object. If False, \
        the function will return a new FeatureSet instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new, filtered instance of the FeatureSet object.
        """
        suffix = '_filtGO'
        go_ids = parsing.data_to_set(go_ids)
        # make sure 'mode' is legal
        assert mode in {'union', 'intersection'}, \
            f"Illegal mode '{mode}': mode must be either 'union' or 'intersection'"

        (taxon_id, organism), gene_id_type = io.get_taxon_and_id_type(organism, gene_id_type, self.index_set,
                                                                      'UniProtKB')

        dag_tree = ontology.fetch_go_basic()
        # make sure all GO IDs are valid
        for go_id in go_ids.copy():
            if go_id in dag_tree:
                continue
            assert isinstance(go_id, str), "'go_ids' must all be strings!"
            stripped_id = go_id.strip()
            if stripped_id in dag_tree:
                go_ids.remove(go_id)
                go_ids.add(stripped_id)
                continue
            raise ValueError(f"Invalid GO ID: '{go_id}'")
        # find the minimal set of GO aspects that need to be fetched
        aspects = {dag_tree[go_id].namespace for go_id in go_ids}

        annotations = io.GOlrAnnotationIterator(taxon_id, aspects,
                                                evidence_types, excluded_evidence_types,
                                                databases, excluded_databases,
                                                qualifiers, excluded_qualifiers)
        assert annotations.n_annotations > 0, "No GO annotations were found for the given parameters. " \
                                              "Please try again with a different set of parameters. "
        go_to_genes = {}
        source_to_genes = {}
        for annotation in tqdm(annotations, desc="Fetching GO annotations", total=annotations.n_annotations,
                               unit=' annotations'):
            # extract gene_id, go_id, source from the annotation. skip annotations that don't appear in the GO DAG

            if annotation['annotation_class'] not in dag_tree:
                continue
            go_id: str = dag_tree[annotation['annotation_class']].id  # replace alt_ids with their main id
            gene_id: str = annotation['bioentity_internal_id']
            source: str = annotation['source']

            # add annotation to annotation dictionary
            if go_id not in go_ids and not propagate_annotations:
                continue
            if go_id not in go_to_genes:
                go_to_genes[go_id] = set()
            go_to_genes[go_id].add(gene_id)

            # add gene id and source to source dict
            if source not in source_to_genes:
                source_to_genes[source] = set()
            source_to_genes[source].add(gene_id)
            # propagate annotations
            if propagate_annotations:
                parents = dag_tree.upper_induced_graph_iter(go_id)
                for parent in parents:
                    if parent not in go_ids:
                        continue
                    if parent not in go_to_genes:
                        go_to_genes[parent] = set()
                    go_to_genes[parent].add(gene_id)

        # translate gene IDs
        go_to_translated_genes = {go_id: set() for go_id in go_ids if go_id in go_to_genes}

        for source in source_to_genes:
            translator = io.GeneIDTranslator(source, gene_id_type).run(parsing.data_to_tuple(source_to_genes[source]))
            for go_id in go_to_genes:
                if go_id not in go_to_translated_genes:
                    continue
                for gene_id in go_to_genes[go_id].copy():
                    if gene_id in translator:
                        go_to_translated_genes[go_id].add(translator[gene_id])
                        go_to_genes[go_id].remove(gene_id)

        # generate list of the indices that are positive for each attribute
        chosen_genes = set()
        # if in union mode, calculate union between the indices of all attributes
        if mode == 'union':
            suffix += 'Union'
            for grp in go_to_translated_genes.values():
                chosen_genes = chosen_genes.union(grp)
            chosen_genes = chosen_genes.intersection(self.index_set)
        # if in intersection mode, calculate intersection between the indices of all attributes
        elif mode == 'intersection':
            suffix += 'Intersection'
            chosen_genes = self.index_set
            for grp in go_to_translated_genes.values():
                chosen_genes = chosen_genes.intersection(grp)

        new_df = self.df.filter(pl.first().is_in(chosen_genes))
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by KEGG Pathways annotations')
    def filter_by_kegg_annotations(self, kegg_ids: Union[str, List[str]],
                                   mode: Literal['union', 'intersection'] = 'union',
                                   organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                                   gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                   opposite: bool = False, inplace: bool = True):
        """
        Filters genes according to KEGG pathways, keeping only genes that belong to specific KEGG pathway. \
        When multiple KEGG IDs are given, filtering can be done in 'union' mode \
        (where genes that belong to at least one pathway are not filtered out), or in 'intersection' mode \
        (where only genes that belong to ALL pathways are not filtered out).

        :param kegg_ids: the KEGG pathway IDs according to which the table will be filtered. \
        An example for a legal KEGG pathway ID would be 'path:cel04020' for the C. elegans calcium signaling pathway.
        :type kegg_ids: str or list of str
        :type mode: 'union' or 'intersection'.
        :param mode: If 'union', filters out every genomic feature that does not belong to one or more \
        of the indicated attributes. If 'intersection', \
        filters out every genomic feature that does not belong to ALL of the indicated attributes.
        param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :type opposite: bool (default=False)
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new, filtered instance of the Filter object.
        """
        suffix = '_filtKEGG'
        kegg_ids = parsing.data_to_set(kegg_ids)
        # make sure 'mode' is legal
        assert mode in {'union', 'intersection'}, \
            f"Illegal mode '{mode}': mode must be either 'union' or 'intersection'"

        (taxon_id, organism), gene_id_type = io.get_taxon_and_id_type(organism, gene_id_type, self.index_set, 'KEGG')

        # fetch KEGG annotations for the specified KEGG IDs
        annotations = io.KEGGAnnotationIterator(taxon_id, list(kegg_ids))
        kegg_to_genes = {}
        for pathway_id, _, genes in tqdm(annotations, desc='Fetching KEGG annotations',
                                         total=annotations.n_annotations, unit=' annotations'):
            if pathway_id in kegg_ids:
                kegg_to_genes[pathway_id] = parsing.data_to_set(genes)
        if len(kegg_to_genes) < len(kegg_ids):
            not_annotated = kegg_ids.difference(parsing.data_to_set(kegg_to_genes))
            warnings.warn(f"Could not find {len(not_annotated)} KEGG IDs: {not_annotated}. \n"
                          f"Proceeding with the other {len(kegg_to_genes)} KEGG IDs. ")

        genes_to_translate = set()
        for grp in kegg_to_genes.values():
            genes_to_translate.update(grp)

        kegg_to_translated_genes = {}
        translator = io.GeneIDTranslator('KEGG', gene_id_type).run(parsing.data_to_tuple(genes_to_translate))
        for kegg_id in kegg_to_genes:
            kegg_to_translated_genes[kegg_id] = set()
            for gene_id in kegg_to_genes[kegg_id]:
                if gene_id in translator:
                    kegg_to_translated_genes[kegg_id].add(translator[gene_id])

        # generate list of the indices that are positive for each attribute
        chosen_genes = set()
        # if in union mode, calculate union between the indices of all attributes
        if mode == 'union':
            suffix += 'Union'
            for grp in kegg_to_translated_genes.values():
                chosen_genes = chosen_genes.union(grp)
            chosen_genes = chosen_genes.intersection(self.index_set)
        # if in intersection mode, calculate intersection between the indices of all attributes
        elif mode == 'intersection':
            suffix += 'Intersection'
            chosen_genes = self.index_set
            for grp in kegg_to_translated_genes.values():
                chosen_genes = chosen_genes.intersection(grp)

        new_df = self.df.filter(pl.first().is_in(chosen_genes))
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by user-defined attribute')
    def filter_by_attribute(self, attributes: Union[str, List[str]] = None,
                            mode: Literal['union', 'intersection'] = 'union',
                            ref: Union[str, Path, Literal['predefined']] = 'predefined',
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
        :type ref: str or pathlib.Path (default='predefined')
        :param ref: filename/path of the attribute reference table to be used as reference.
        :type opposite: bool (default=False)
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
        else:
            attributes = parsing.data_to_list(attributes)
        # make sure 'mode' is legal
        assert mode in {'union', 'intersection'}, \
            f"Illegal mode '{mode}': mode must be either 'union' or 'intersection'"
        # load the Attribute Reference Table
        attr_ref_table = io.load_table(settings.get_attr_ref_path(ref))
        validation.validate_attr_table(attr_ref_table)
        attr_ref_table.fill_nan(None)
        # generate list of the indices that are positive for each attribute
        attr_indices_list = [
            {i[0] for i in attr_ref_table.filter(pl.col(attr).is_not_null()).select(pl.first()).iter_rows()} for attr in
            attributes]
        chosen_genes = set()
        # if in union mode, calculate union between the indices of all attributes
        if mode == 'union':
            for idx in attr_indices_list:
                chosen_genes = chosen_genes.union(idx)
        # if in intersection mode, calculate intersection between the indices of all attributes
        elif mode == 'intersection':
            chosen_genes = attr_indices_list[0]
            for idx in attr_indices_list[1:]:
                chosen_genes = chosen_genes.intersection(idx)

        new_df = self.df.filter(pl.first().is_in(chosen_genes))

        if len(attributes) > 1:
            suffix += mode.capitalize()
        else:
            suffix += str(attributes[0])
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Split by user-defined attribute')
    def split_by_attribute(self, attributes: Union[str, List[str]],
                           ref: Union[str, Path, Literal['predefined']] = 'predefined') -> tuple:

        """
        Splits the features in the Filter object into multiple Filter objects, \
        each corresponding to one of the specified Attribute Reference Table attributes. \
        Each new Filter object will contain only features that belong to its Attribute Reference Table attribute.

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

    @readable_name('Table descriptive statistics')
    def describe(self, percentiles: Union[float, List[float]] = (0.01, 0.25, 0.5, 0.75, 0.99)) -> pl.DataFrame:

        """
        Generate descriptive statistics that summarize the central tendency, dispersion and shape \
        of the dataset's distribution, excluding NaN values. \
        For more information see the documentation of pandas.DataFrame.describe.

        :type percentiles: list-like of floats (default=(0.01, 0.25, 0.5, 0.75, 0.99))
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
        return self.df.describe(percentiles=parsing.data_to_tuple(percentiles))

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
        data = self.df.select(pl.first())
        if data.is_duplicated().any():
            warnings.warn(" this filter object contains multiple rows with the same gene name/ID. When "
                          "returning a set or string of features from this table, each gene ID will "
                          "appear ONLY ONCE!")
        return {item[0] for item in data.iter_rows()}

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
        return "\n".join((str(ind) for ind in self.df.select(pl.first()).to_series().to_list()))

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

    @readable_name('Summarize feature biotypes (based on a GTF file)')
    def biotypes_from_gtf(self, gtf_path: Union[str, Path],
                          attribute_name: Union[Literal[BIOTYPE_ATTRIBUTE_NAMES], str] = 'gene_biotype',
                          feature_type: Literal['gene', 'transcript'] = 'gene',
                          long_format: bool = False) -> pl.DataFrame:

        """
        Returns a DataFrame describing the biotypes in the table and their count. \
        The data about feature biotypes is drawn from a GTF (Gene transfer format) file supplied by the user.

        :param gtf_path: Path to your GTF (Gene transfer format) file. The file should match the type of \
        gene names/IDs you use in your table, and should contain an attribute describing biotype.
        :type gtf_path: str or Path
        :param attribute_name: name of the attribute in your GTF file that describes feature biotype.
        :type attribute_name: str (default='gene_biotype')
        :param feature_type: determined whether the features/rows in your data table describe \
        individual genes or transcripts.
        :type feature_type: 'gene' or 'transcript' (default='gene')
        :type long_format: bool (default=False)
        :param long_format:if True, returns a short-form DataFrame, which states the biotypes \
        in the Filter object and their count. Otherwise, returns a long-form DataFrame,
        which also provides descriptive statistics of each column per biotype.
        :rtype: pandas.DataFrame
        :returns: a pandas DataFrame showing the number of values belonging to each biotype, \
        as well as additional descriptive statistics of format=='long'.
        """
        # load the Biotype Reference Table
        ref_df = self._get_ref_srs_from_gtf(gtf_path, attribute_name, feature_type)
        return self._biotypes(ref_df, long_format)

    def _biotypes(self, ref_df: pl.DataFrame, long_format: bool):
        # find which genes from tne Filter object don't appear in the Biotype Reference Table
        not_in_ref = parsing.data_to_set(self.df.select(pl.first())).difference(
            parsing.data_to_set(ref_df.select(pl.first())))
        if len(not_in_ref) > 0:
            warnings.warn(
                f'{len(not_in_ref)} of the features in the table do not appear in the Biotype Reference Table')
            missing_df = pl.DataFrame(
                [parsing.data_to_list(not_in_ref), ['_missing_from_biotype_reference'] * len(not_in_ref)],
                schema=ref_df.columns)
            ref_df = pl.concat([ref_df, missing_df])

        if long_format:
            # additionally return descriptive statistics for each biotype
            df_with_biotype = self.df.join(ref_df, left_on=self.df.columns[0], right_on=ref_df.columns[0], how='left')
            res = df_with_biotype.to_pandas().groupby('biotype').describe()
            res.columns = ['_'.join(col) for col in res.columns]
            res = pl.from_pandas(res.reset_index())
            return res
        else:
            # return just the number of genes/indices belonging to each biotype
            return ref_df.filter(pl.first().is_in(self.df.select(pl.first()))).group_by(pl.last()).count()

    @readable_name('Summarize feature biotypes (based on a reference table)')
    def biotypes_from_ref_table(self, long_format: bool = False,
                                ref: Union[str, Path, Literal['predefined']] = 'predefined') -> pl.DataFrame:

        """
        Returns a DataFrame describing the biotypes in the table and their count. \
        The data about feature biotypes is drawn from a Biotype Reference Table supplied by the user.

        :type long_format: bool (default=False)
        :param long_format:if True, returns a short-form DataFrame, which states the biotypes \
        in the Filter object and their count. Otherwise, returns a long-form DataFrame,
        which also provides descriptive statistics of each column per biotype.
        :param ref: Name of the biotype reference table used to determine biotype. Default is ce11 (included in the package).
        :rtype: pandas.DataFrame
        :returns: a pandas DataFrame showing the number of values belonging to each biotype, \
        as well as additional descriptive statistics of format=='long'.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> # short-form view
            >>> d.biotypes_from_ref_table(ref='tests/biotype_ref_table_for_tests.csv')
                            gene
            biotype
            protein_coding    26
            pseudogene         1
            unknown            1

            >>> # long-form view
            >>> d.biotypes_from_ref_table(long_format=True,ref='tests/biotype_ref_table_for_tests.csv')
                           baseMean               ...           padj
                              count         mean  ...            75%            max
            biotype                               ...
            protein_coding     26.0  1823.089609  ...   1.005060e-90   9.290000e-68
            pseudogene          1.0  2688.043701  ...   1.800000e-94   1.800000e-94
            unknown             1.0  2085.995094  ...  3.070000e-152  3.070000e-152
            <BLANKLINE>
            [3 rows x 48 columns]

        """
        # load the Biotype Reference Table
        ref = settings.get_biotype_ref_path(ref)
        ref_df = io.load_table(ref).drop_nulls()
        validation.validate_biotype_table(ref_df)
        return self._biotypes(ref_df, long_format)

    @readable_name('Filter with a number filter')
    def number_filters(self, column: param_typing.ColumnName,
                       operator: Literal['greater than', 'equals', 'lesser than', 'abs greater than'], value: float,
                       opposite: bool = False, inplace: bool = True):

        """
        Applay a number filter (greater than, equal, lesser than) on a particular column in the Filter object.

        :type column: str
        :param column: name of the column to filter by
        :type operator: str: 'gt' / 'greater than' / '>', 'eq' / 'equals' / '=', 'lt' / 'lesser than' / '<'
        :param operator: the operator to filter the column by (greater than, equal or lesser than)
        :type value: float
        :param value: the value to filter by
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
                         'lesser than': 'lt', '<': 'lt', 'equal': 'eq', 'abs greater than': 'abs_gt',
                         'abs_gt': 'abs_gt', '|x|>': 'abs_gt'}
        operator = operator.lower()
        assert operator in operator_dict, f"Invalid operator {operator}"
        op = operator_dict[operator]
        # determine that 'value' is a number
        assert isinstance(value, (int, float)), "'value' must be a number!"
        # determine that the column is legal
        assert column in self.columns, f"column {column} not in DataFrame!"

        suffix = f"_{column}{op}{value}"
        # perform operation according to operator
        if op == 'eq':
            new_df = self.df.filter(pl.col(column) == value)
        elif op == 'gt':
            new_df = self.df.filter(pl.col(column) > value)
        elif op == 'lt':
            new_df = self.df.filter(pl.col(column) < value)
        elif op == 'abs_gt':
            new_df = self.df.filter(pl.col(column).abs() > value)

        # noinspection PyUnboundLocalVariable
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter with a text filter')
    def text_filters(self, column: param_typing.ColumnName,
                     operator: Literal['equals', 'contains', 'starts with', 'ends with'], value: str,
                     opposite: bool = False, inplace: bool = True):

        """
        Applay a text filter (equals, contains, starts with, ends with) on a particular column in the Filter object.

        :type column: str
        :param column: name of the column to filter by
        :type operator: str: 'eq' / 'equals' / '=', 'ct' / 'contains' / 'in', 'sw' / 'starts with', 'ew' / 'ends with'
        :param operator: the operator to filter the column by (equals, contains, starts with, ends with)
        :type value: number (int or float)
        :param value: the value to filter by
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
        assert isinstance(value, str), "'value' must be a string!"
        # determine that the column is legal
        assert column in self.columns, f"column {column} not in DataFrame!"

        suffix = f"_{column}{op}{value}"
        # perform operation according to operator
        if op == 'eq':
            new_df = self.df.filter(pl.col(column) == value)
        elif op == 'ct':
            new_df = self.df.filter(pl.col(column).str.contains(value))
        elif op == 'ew':
            new_df = self.df.filter(pl.col(column).str.ends_with(value))
        elif op == 'sw':
            new_df = self.df.filter(pl.col(column).str.starts_with(value))

        # noinspection PyUnboundLocalVariable
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Remove rows with missing values')
    def filter_missing_values(self, columns: Union[param_typing.ColumnNames, Literal['all']] = 'all',
                              opposite: bool = False, inplace: bool = True):
        """
        Remove all rows whose values in the specified columns are missing (NaN).

        :type columns: str or list of str (default='all')
        :param columns:name/names of the columns to check for missing values.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
            >>> filt_no_nan_basemean_pval = filt.filter_missing_values(columns=['baseMean','pvalue'], inplace=False)
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
            new_df = self.df.drop_nulls(subset=subset)
        else:
            new_df = self.df.drop_nulls()

        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Apply a transformation to the table')
    def transform(self, function: Union[Literal['Box-Cox', 'log2', 'log10', 'ln', 'Standardize'], Callable],
                  columns: Union[param_typing.ColumnNames, Literal['all']] = 'all',
                  inplace: bool = True, **function_kwargs):
        """
        Transform the values in the Filter object with the specified function.

        :param function: The function or function name to be applied.
        :type function: Callable or str ('logx' for base-x log of the data + 1, \
        'box-cox' for Box-Cox transform of the data + 1, 'standardize' for standardization)
        :param columns: The columns to which the transform should be applied.
        :type columns: str, list of str, or 'all' (default='all')
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :param function_kwargs: Any additional keyworded arguments taken by the supplied function.
        :return: If 'inplace' is False, returns a new instance of the Filter object.

        :Examples:
            >>> from rnalysis import filtering
            >>> filt = filtering.Filter('tests/test_files/counted.csv')
            >>> filt_log10 = filt.transform('log10', inplace=False)
            Transformed 22 features. Transformation result saved to new object.
            >>> filt.transform(lambda x: x+1, columns=['cond1','cond4'])
            Transformed 22 features. Transformed inplace.
        """

        predefined_funcs = {'ln': np.log1p, 'loge': np.log1p, 'standardize': StandardScaler.fit_transform,
                            'box-cox': lambda x: PowerTransformer(method='box-cox').fit_transform(x + 1)}

        suffix = "_customtransform" if callable(function) else f"_{str(function)}transform"

        if isinstance(columns, str) and columns.lower() == 'all':
            columns = parsing.data_to_list(self.columns)
        else:
            columns = parsing.data_to_list(columns)

        if isinstance(function, str):
            function = function.lower()
            log_match = re.findall(r"log[\d]+", function)
            if len(log_match) == 1:
                log_base = int(log_match[0][3:])
                new_df = self.df.with_columns(self.df.select(pl.col(columns).add(1).log(log_base)))

            elif function in predefined_funcs:
                function = predefined_funcs[function]
                new_df = self.df.with_columns(
                    pl.DataFrame(function(self.df.select(pl.col(columns)).to_numpy()), schema=columns))
            else:
                raise TypeError(f"Invalid value for 'function': '{function}'. ")
        elif callable(function):
            partial = functools.partial(function, **function_kwargs)
            new_df = self.df.with_columns(self.df.select(pl.col(columns).map_elements(partial)))
        else:
            raise TypeError(f"Invalid value for 'function': {function} (type {type(function)}). "
                            f"Please specify a function or a recognized function name. ")
        return self._inplace(new_df, False, inplace, suffix, 'transform')

    @readable_name('Sort table rows')
    def sort(self, by: Union[str, List[str]], ascending: Union[bool, List[bool]] = True, na_position: str = 'last',
             inplace: bool = True):

        """
        Sort the rows by the values of specified column or columns.

        :type by: str or list of str
        :param by: Names of the column or columns to sort by.
        :type ascending: bool or list of bool (default=True)
        :param ascending: Sort ascending vs. descending. Specify list for multiple sort orders. \
        If this is a list of bools, it must have the same length as 'by'.
        :type na_position: 'first' or 'last', default 'last'
        :param na_position: If 'first', puts NaNs at the beginning; if 'last', puts NaNs at the end.
        :type inplace: bool (default=True)
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
        new_df = self._sort(by=by, ascending=ascending, na_position=na_position)
        return self._inplace(new_df, False, inplace, suffix, 'sort')

    def _sort(self, by: Union[str, List[str]], ascending: Union[bool, List[bool]] = True, na_position: str = 'last'):
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
        if isinstance(ascending, bool):
            descending = not ascending
        else:
            descending = [not item for item in parsing.data_to_list(ascending)]
        return self.df.sort(by=by, descending=descending, nulls_last=na_position == 'last')

    @readable_name('Filter all but top N values')
    def filter_top_n(self, by: param_typing.ColumnNames, n: PositiveInt = 100,
                     ascending: Union[bool, List[bool]] = True, na_position: str = 'last', opposite: bool = False,
                     inplace: bool = True, ):

        """
        Sort the rows by the values of specified column or columns, then keep only the top 'n' rows.

        :type by: name of column/columns (str/List[str])
        :param by: Names of the column or columns to sort and then filter by.
        :type n: int
        :param n: How many features to keep in the Filter object.
        :type ascending: bool or list of bools (default=True)
        :param ascending: Sort ascending vs. descending. Specify list for multiple sort orders. \
        If this is a list of bools, it must have the same length as 'by'.
        :type na_position: 'first' or 'last', default 'last'
        :param na_position: If 'first', puts NaNs at the beginning; if 'last', puts NaNs at the end.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
        if n > self.shape[0]:
            warnings.warn(f'Current number of rows {self.shape[0]} is smaller than the specified n={n}. '
                          f'Therefore output Filter object will only have {self.shape[0]} rows. ')
            n = self.shape[0]
        new_df = self._sort(by=by, ascending=ascending, na_position=na_position)[0:n, :]
        return self._inplace(new_df, opposite, inplace, suffix)

    @staticmethod
    def _return_type(index_set: set, return_type: Literal['set', 'str']):
        assert isinstance(return_type, str), "'return_type' must be a string!!"
        if return_type == 'set':
            return index_set
        elif return_type == 'str':
            return "\n".join(sorted(index_set))
        else:
            raise ValueError(f"'return type' must be either 'set' or 'str', instead got '{return_type}'.")

    def _set_ops(self, others, return_type: Literal['set', 'str'], op: Any, **kwargs):
        """
        Apply the supplied set operation (union/intersection/difference/symmetric difference) to the supplied objects.

        :param others: the other objects to apply the set operation to
        :type others: Filter or set objects.
        :param return_type: the return type of the output
        :type return_type: 'set' or 'str'
        :param op: the set operation
        :type op: function (set.union, set.intersection, set.difference or set.symmetric_difference)
        :param kwargs: any additional keyworded arguments to be supplied to the set operation.
        :return: a set/string of indices resulting from the set operation
        :rtype: set or str
        """
        others = parsing.data_to_list(others).copy()

        for i, other in enumerate(others):
            if isinstance(other, Filter):
                others[i] = set(other.df.select(pl.col(other.df.columns[0])).to_series().to_list())
            elif isinstance(other, set):
                pass
            else:
                raise TypeError(f"'others' must contain only Filter objects or sets, "
                                f"instead got object {other} of type {type(other)}.")
        try:
            op_indices = op(set(self.df.select(pl.col(self.df.columns[0])).to_series().to_list()), *others, **kwargs)
        except TypeError as e:
            if op == set.symmetric_difference:
                raise TypeError(
                    f"Symmetric difference can only be calculated for two objects, {len(others) + 1} were given!")
            else:
                raise e
        return Filter._return_type(op_indices, return_type)

    def intersection(self, *others: Union['Filter', set], return_type: Literal['set', 'str'] = 'set',
                     inplace: bool = False):

        """
        Keep only the features that exist in ALL of the given Filter objects/sets. \
        Can be done either inplace on the first Filter object, or return a set/string of features.

        :type others: Filter or set objects.
        :param others: Objects to calculate intersection with.
        :type return_type: 'set' or 'str' (default='set')
        :param return_type: If 'set', returns a set of the intersecting features. If 'str', returns a string of \
        the intersecting features, delimited by a comma.
        :type inplace: bool (default=False)
        :param inplace: If True, the function will be applied in-place to the current Filter object. \
        If False (default), the function will return a set/str that contains the intersecting indices.
        :rtype: set or str
        :return: If inplace=False, returns a set/string of the features that \
        intersect between the given Filter objects/sets.


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
            new_set = self._set_ops(others, 'set', set.intersection)

            suffix = "_intersection"
            return self._inplace(self.df.filter(pl.first().is_in(new_set)), opposite=False, inplace=inplace,
                                 suffix=suffix)
        # if intersection is not performed inplace, return a set/string according to user's request
        new_set = self._set_ops(others, return_type, set.intersection)
        return new_set

    def majority_vote_intersection(self, *others: Union['Filter', set], majority_threshold: float = 0.5,
                                   return_type: Literal['set', 'str'] = 'set'):

        """
        Returns a set/string of the features that appear in at least \
        (majority_threhold * 100)% of the given Filter objects/sets. \
        Majority-vote intersection with majority_threshold=0 is equivalent to Union. \
        Majority-vote intersection with majority_threshold=1 is equivalent to Intersection.

        :type others: Filter or set objects.
        :param others: Objects to calculate intersection with.
        :type return_type: 'set' or 'str' (default='set')
        :param return_type: If 'set', returns a set of the intersecting WBGene indices. If 'str', returns a string of \
        the intersecting indices, delimited by a comma.
        :type majority_threshold: float (default=0.5)
        :param majority_threshold: The threshold that determines what counts as majority. Features will be returned \
        only if they appear in at least (majority_threshold * 100)% of the given Filter objects/sets.
        :rtype: set or str
        :return: If inplace=False, returns a set/string of the features that \
        uphold majority vote intersection between two given Filter objects/sets.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> a_set = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> b_set = {'WBGene00000002','WBGene00000004'}
            >>> # calculate majority-vote intersection and return a set
            >>> d.majority_vote_intersection(a_set, b_set, majority_threshold=2/3)
            {'WBGene00000002', 'WBGene00000003', 'WBGene00000004'}

        """
        new_set = self._set_ops(others, return_type, generic.SetWithMajorityVote.majority_vote_intersection,
                                majority_threshold=majority_threshold)
        return new_set

    def union(self, *others: Union['Filter', set], return_type: Literal['set', 'str'] = 'set'):

        """
        Returns a set/string of the union of features between multiple Filter objects/sets \
        (the features that exist in at least one of the Filter objects/sets).

        :type others: Filter or set objects.
        :param others: Objects to calculate union with.
        :type return_type: 'set' or 'str' (default='set')
        :param return_type: If 'set', returns a set of the union features. If 'str', returns a string of \
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

    def difference(self, *others: Union['Filter', set], return_type: Literal['set', 'str'] = 'set',
                   inplace: bool = False):

        """
        Keep only the features that exist in the first Filter object/set but NOT in the others. \
        Can be done inplace on the first Filter object, or return a set/string of features.

        :type others: Filter or set objects.
        :param others: Objects to calculate difference with.
        :type return_type: 'set' or 'str' (default='set')
        :param return_type: If 'set', returns a set of the features that exist only in the first Filter object. \
        If 'str', returns a string of the WBGene indices that exist only in the first Filter object, \
        delimited by a comma.
        :type inplace: bool, default False
        :param inplace: If True, filtering will be applied to the current Filter object. If False (default), \
        the function will return a set/str that contains the intersecting features.
        :rtype: set or str
        :return: If inplace=False, returns a set/string of the features\
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
            new_set = self._set_ops(others, 'set', set.difference)
            suffix = "_difference"
            return self._inplace(self.df.filter(pl.first().is_in(new_set)), opposite=False, inplace=inplace,
                                 suffix=suffix)
            # if difference is not performed inplace, return a set/string according to user's request
        new_set = self._set_ops(others, return_type, set.difference)
        return new_set

    def symmetric_difference(self, other: Union['Filter', set], return_type: Literal['set', 'str'] = 'set'):

        """
        Returns a set/string of the WBGene indices that exist either in the first Filter object/set OR the second, \
        but NOT in both (set symmetric difference).

        :type other: Filter or set.
        :param other: a second Filter object/set to calculate symmetric difference with.
        :type return_type: 'set' or 'str' (default='set')
        :param return_type: If 'set', returns a set of the features that exist in exactly one Filter object. \
        If 'str', returns a string of the features that exist in exactly one Filter object, delimited by a comma.
        :rtype: set or str
        :return: a set/string of the features that that exist t in exactly one Filter. (set symmetric difference).


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


@readable_name('Fold change table')
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

    def __init__(self, fname: Union[str, Path, tuple], numerator_name: str, denominator_name: str,
                 suppress_warnings: bool = False):
        """
        Load a fold-change table. Valid fold-change tables should contain exactly two columns: \
        the first column containing gene names/indices, and the second column containing log2(fold change) values.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param numerator_name: name of the numerator condition in the fold-change table
        :type numerator_name: str
        :param denominator_name: name of the denominator condition in the fold-change table
        :type denominator_name: str
        :param suppress_warnings: if True, RNAlysis will not issue warnings about the loaded table's \
        structure or content.
        :type suppress_warnings: bool (default=False)
        """
        super().__init__(fname, suppress_warnings=True)
        self.numerator = numerator_name
        self.denominator = denominator_name
        if not suppress_warnings:
            self._init_warnings()

    def _init_warnings(self):
        super()._init_warnings()
        # inf/0 can be problematic for functions down the line (like randomization test)
        if np.inf in self.df or 0 in self.df:
            warnings.warn(
                " FoldChangeFilter does not support 'inf' or '0' values! "
                "Unexpected results may occur during filtering or statistical analyses. ")

    def __repr__(self):
        return f"{type(self).__name__}('{self.fname.as_posix()}', '{self.numerator}', '{self.denominator}')"

    def __str__(self):
        return f"{type(self).__name__} (numerator: '{self.numerator}', denominator: '{self.denominator}') " \
               f"of file {self.fname.name}"

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame, name: Union[str, Path], numerator_name: str = 'numerator',
                       denominator_name: str = 'denominator', suppress_warnings: bool = False) -> 'FoldChangeFilter':
        obj = cls.__new__(cls)
        obj.df = df.clone()
        obj.fname = Path(name)
        obj.numerator = numerator_name
        obj.denominator = denominator_name
        if not suppress_warnings:
            obj._init_warnings()
        return obj

    def __copy__(self):
        return type(self).from_dataframe(self.df, self.fname, numerator_name=self.numerator,
                                         denominator_name=self.denominator)

    @readable_name('Perform randomization test')
    def randomization_test(self, ref, alpha: param_typing.Fraction = 0.05, reps: PositiveInt = 10000,
                           save_csv: bool = False,
                           fname: Union[str, None] = None, random_seed: Union[int, None] = None) -> pl.DataFrame:

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
        :type random_seed: non-negative integer (default=None)
        :type random_seed: The random seed used to initialize the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs.
        :rtype: pandas DataFrame
        :return: A Dataframe with the number of given genes, the observed fold change for the given group of genes, \
        the expected fold change for a group of genes of that size and the p value for the comparison.


        :Examples:
            >>> from rnalysis import filtering
            >>> f = filtering.FoldChangeFilter('tests/test_files/fc_1.csv' , 'numerator' , 'denominator')
            >>> f_background = f.filter_biotype_from_ref_table('protein_coding', ref='tests/biotype_ref_table_for_tests.csv', inplace=False) #keep only protein-coding genes as reference
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
        obs_fc = self.df['Fold Change'].mean()
        exp_fc = ref.df['Fold Change'].mean()
        n = self.shape[0]
        # set random seed if requested
        if random_seed is not None:
            assert isinstance(random_seed, int) and random_seed >= 0, f"random_seed must be a non-negative integer. " \
                                                                      f"Value {random_seed} invalid."
            np.random.seed(random_seed)
        # run randomization test
        print('Calculating...')
        res = self._foldchange_randomization(ref.df['Fold Change'].to_numpy(), reps, obs_fc, exp_fc, n)
        # format the output DataFrame
        res_df = pl.DataFrame(res, schema=['group size', 'observed fold change', 'expected fold change', 'pval'])
        res_df = res_df.with_columns(significant=res_df.select('pval') <= alpha)
        # save the output DataFrame if requested
        if save_csv:
            io.save_table(res_df, fname)
        print(res_df)

        return res_df

    @staticmethod
    @generic.numba.jit(nopython=True)
    def _foldchange_randomization(vals: np.ndarray, reps: PositiveInt, obs_fc: float, exp_fc: float,
                                  n: int):  # pragma: no cover
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

    @readable_name('Remove rows with missing values')
    def filter_missing_values(self, opposite: bool = False, inplace: bool = True):
        """
        Remove all rows with missing values.

        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of the Filter object.
        """
        return super().filter_missing_values(columns=self.columns[-1], opposite=opposite, inplace=inplace)

    @readable_name('Filter by absolute log2 fold-change magnitude')
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
        :type inplace: bool (default=True)
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
        new_df = self.df.filter(pl.last().log(2).abs() >= abslog2fc)
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by fold-change direction')
    def filter_fold_change_direction(self, direction: Literal['pos', 'neg'] = 'pos', opposite: bool = False,
                                     inplace: bool = True):
        """
        Filters out features according to the direction in which they changed between the two conditions.

        :param direction: 'pos' or 'neg'. If 'pos', will keep only features that have positive log2foldchange. \
        If 'neg', will keep only features that have negative log2foldchange.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
            new_df = self.df.filter(pl.last() > 1)
            suffix = '_PositiveLog2FC'
        elif direction == 'neg':
            new_df = self.df.filter(pl.last() < 1)
            suffix = '_NegativeLog2FC'
        else:
            raise ValueError(
                "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. ")
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Split by fold-change direction')
    def split_fold_change_direction(self) -> tuple:
        """
        Splits the features in the FoldChangeFilter object into two non-overlapping \
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


@readable_name('Differential expression table')
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
    __slots__ = {'log2fc_col': 'name of the log2 fold change column', 'padj_col': 'name of the adjusted p-value column',
                 'pval_col': 'name of the p-value column'}

    def __init__(self, fname: Union[str, Path, tuple], drop_columns: Union[str, List[str]] = None,
                 log2fc_col: Union[str, Literal['log2FoldChange', 'logFC']] = 'log2FoldChange',
                 padj_col: Union[str, Literal['padj', 'adj.P.Val']] = 'padj',
                 pval_col: Union[str, Literal['pvalue', 'P.Value']] = 'pvalue',
                 suppress_warnings: bool = False):

        """
        Load a differential expression table. A valid differential expression table should have \
        a column containing log2(fold change) values for each gene, and another column containing \
        adjusted p-values for each gene.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded table.
        :type drop_columns: str, list of str, or None (default=None)
        :param log2fc_col: name of the table column containing log2(fold change) values.
        :type log2fc_col: str (default='Log2FoldChange')
        :param padj_col: name of the table column containing adjusted p-values.
        :type padj_col: str (default='padj')
        :param suppress_warnings: if True, RNAlysis will not issue warnings about the loaded table's \
        structure or content.
        :type suppress_warnings: bool (default=False)
        """
        super().__init__(fname, drop_columns, suppress_warnings=True)
        self.log2fc_col = log2fc_col
        self.padj_col = padj_col
        self.pval_col = pval_col
        if not suppress_warnings:
            self._init_warnings()

    def _init_warnings(self):
        super()._init_warnings()
        if self.log2fc_col not in self.columns:
            warnings.warn(f"The specified log2fc_col '{self.log2fc_col}' does not appear in the DESeqFilter's columns: "
                          f"{self.columns}. DESeqFilter-specific functions that depend on "
                          f"log2(fold change) may fail to run. ")
        if self.padj_col not in self.columns:
            warnings.warn(f"The specified padj_col '{self.padj_col}' does not appear in the DESeqFilter's columns: "
                          f"{self.columns}. DESeqFilter-specific functions that depend on p-values may fail to run. ")
        if self.pval_col not in self.columns:
            warnings.warn(f"The specified pval_col '{self.pval_col}' does not appear in the DESeqFilter's columns: "
                          f"{self.columns}. DESeqFilter-specific functions that depend on p-values may fail to run. ")

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame, name: Union[str, Path],
                       log2fc_col: Union[str, Literal['log2FoldChange', 'logFC']] = 'log2FoldChange',
                       padj_col: Union[str, Literal['padj', 'adj.P.Val']] = 'padj',
                       pval_col: Union[str, Literal['pvalue', 'P.Value']] = 'pvalue',
                       suppress_warnings: bool = False) -> 'DESeqFilter':
        obj = cls.__new__(cls)
        obj.df = df.clone()
        obj.fname = Path(name)
        obj.log2fc_col = log2fc_col
        obj.padj_col = padj_col
        obj.pval_col = pval_col
        if not suppress_warnings:
            obj._init_warnings()
        return obj

    def __copy__(self):
        return type(self).from_dataframe(self.df, self.fname, self.log2fc_col, self.padj_col, self.pval_col,
                                         suppress_warnings=True)

    def _assert_padj_col(self):
        if self.padj_col not in self.columns:
            raise KeyError(f"A column with adjusted p-values under the name padj_col='{self.padj_col}' "
                           f"could not be found. Try setting a different value for the parameter 'padj_col' "
                           f"when creating the DESeqFilter object.")

    def _assert_pval_col(self):
        if self.pval_col not in self.columns:
            raise KeyError(f"A column with p-values under the name pval_col='{self.pval_col}' "
                           f"could not be found. Try setting a different value for the parameter 'pval_col' "
                           f"when creating the DESeqFilter object.")

    def _assert_log2fc_col(self):
        if self.log2fc_col not in self.columns:
            raise KeyError(f"A column with log2 fold change values under the name log2fc_col='{self.log2fc_col}' "
                           f"could not be found. Try setting a different value for the parameter 'log2fc_col' "
                           f"when creating the DESeqFilter object.")

    @readable_name('Filter by statistical significance')
    def filter_significant(self, alpha: param_typing.Fraction = 0.1, opposite: bool = False, inplace: bool = True):

        """
        Removes all features which did not change significantly, according to the provided alpha.

        :param alpha: the significance threshold to determine which genes will be filtered. between 0 and 1.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
        self._assert_padj_col()
        new_df = self.df.filter(pl.col(self.padj_col) <= alpha)
        suffix = f"_sig{alpha}"
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Plot histogram of p-values')
    def pval_histogram(self, adjusted_pvals: bool = False, bin_size: Fraction = 0.05,
                       title: Union[str, Literal['auto']] = 'auto') -> plt.Figure:
        """
        Plots a histogram of the p-values in the DESeqFilter object. \
        This is often used to troubleshoot the results of a differential expression analysis. \
        For more information about interpreting p-value histograms, see the following blog post by /
        David Robinson: https://varianceexplained.org/statistics/interpreting-pvalue-histogram/

        :param adjusted_pvals: if True, will plot a histogram of the adjusted p-values instead of the raw p-values.
        :type adjusted_pvals: bool (default=False)
        :param bin_size: determines the size of the bins in the histogram.
        :type bin_size: float between 0 and 1 (default=0.05)
        :param title: the title of the histogram. If 'auto', will be set automatically.
        :type title: str or 'auto' (default='auto')
        """
        if adjusted_pvals:
            self._assert_padj_col()
            data = self.df.select(pl.col(self.padj_col))
            title = 'Histogram of adjusted p-values' if title == 'auto' else title
            x_label = 'adjusted p-value'
        else:
            self._assert_pval_col()
            data = self.df[self.pval_col]
            title = 'Histogram of p-values' if title == 'auto' else title
            x_label = 'p-value'
        fig, ax = plt.subplots()
        ax.hist(data, bins=int(1 / bin_size), color='skyblue', edgecolor='black', range=(0, 1))
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel('Frequency')
        ax.set_xticks(np.arange(0, 1.1, 0.1))
        plt.show()
        return fig

    @readable_name('Filter by absolute log2 fold-change magnitude')
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
        :type inplace: bool (default=True)
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
        self._assert_log2fc_col()

        suffix = f"_{abslog2fc}abslog2foldchange"
        new_df = self.df.filter(pl.col(self.log2fc_col).abs() >= abslog2fc)
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by log2 fold-change direction')
    def filter_fold_change_direction(self, direction: Literal['pos', 'neg'] = 'pos', opposite: bool = False,
                                     inplace: bool = True):

        """
        Filters out features according to the direction in which they changed between the two conditions.

        :param direction: 'pos' or 'neg'. If 'pos', will keep only features that have positive log2foldchange. \
        If 'neg', will keep only features that have negative log2foldchange.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
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
        self._assert_log2fc_col()

        if direction == 'pos':
            new_df = self.df.filter(pl.col(self.log2fc_col) > 0)
            suffix = '_PositiveLog2FC'
        elif direction == 'neg':
            new_df = self.df.filter(pl.col(self.log2fc_col) < 0)
            suffix = '_NegativeLog2FC'
        else:
            raise ValueError(
                "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. ")
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Split by log2 fold-change direction')
    def split_fold_change_direction(self) -> tuple:

        """
        Splits the features in the DESeqFilter object into two non-overlapping DESeqFilter \
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

    @readable_name('Volcano plot')
    def volcano_plot(self, alpha: param_typing.Fraction = 0.1, log2fc_threshold: Union[float, None] = 1,
                     title: Union[str, Literal['auto']] = 'auto', title_fontsize: float = 16,
                     label_fontsize: float = 16, tick_fontsize: float = 12,
                     annotation_fontsize: float = 10, point_size: float = 10, opacity: param_typing.Fraction = 0.65,
                     interactive: bool = True, show_cursor: bool = False) -> plt.Figure:

        """
        Plots a volcano plot (log2(fold change) vs -log10(adj. p-value)) of the DESeqFilter object. \
        Significantly upregulated features are colored in red, \
        and significantly downregulated features are colored in blue. \
        If the plot is generated in interactive mode, data points can be labeled by clicking on them.

        :type alpha: float between 0 and 1
        :param alpha: the significance threshold to paint data points as significantly up/down-regulated.
        :param log2fc_threshold: the absolute log2(fold change) threshold to paint data as \
        significantly up/down-regulated. if log2fc_threshold is None, no threshold will be used.
        :type log2fc_threshold: non-negative float or None (default=1)
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
         :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :param annotation_fontsize: determines the font size of the point annotations created in interactive mode.
        :type annotation_fontsize: float (default=10)
        :param opacity: float between 0 and 1 (default=0.65)
        :type opacity: determines the opacity of the points in the scatter plot. 0 indicates completely transparent, \
        while 1 indicates completely opaque.
        :param point_size: determines the size of the points in the scatter plot
        :type point_size: float (default=10)
        :param interactive: if True, turns on interactive mode. While in interactive mode, you can click on a data \
        point to label it with its gene name/ID, or click on a labeled data point to unlabel it.
        :type interactive: bool (default=True)
        :param show_cursor: if True, show the cursor position on the plot during interactive mode
        :type show_cursor: bool (default=False)
        :rtype: A matplotlib Figure

        .. figure:: /figures/volcano.png
           :align:   center
           :scale: 70 %

           Example plot of volcano_plot()

        """
        self._assert_padj_col()
        self._assert_log2fc_col()

        if log2fc_threshold is None:
            log2fc_threshold = 0
        else:
            assert isinstance(log2fc_threshold, (int, float)) and log2fc_threshold >= 0, \
                "'log2fc_threshold' must be a non-negative number!"

        if interactive:
            fig = plt.figure(constrained_layout=True, FigureClass=generic.InteractiveScatterFigure,
                             labels=parsing.data_to_tuple(self.df.select(pl.first())),
                             annotation_fontsize=annotation_fontsize,
                             show_cursor=show_cursor)
            ax = fig.axes[0]
            kwargs = {'picker': 5}
        else:
            fig = plt.figure(constrained_layout=True)
            ax = fig.add_subplot(111)
            kwargs = {}
        colors = self.df.with_columns(
            pl.when((pl.col(self.padj_col) <= alpha).and_(pl.col(self.log2fc_col) >= log2fc_threshold)).then(
                pl.lit('tab:red')).when(
                ((pl.col(self.padj_col) <= alpha).and_(pl.col(self.log2fc_col) <= -log2fc_threshold))).then(
                pl.lit('tab:blue')).otherwise(pl.lit('grey')).alias('color')).select(
            pl.col('color')).to_series().to_list()
        ax.scatter(self.df[self.log2fc_col], -np.log10(self.df[self.padj_col]),
                   c=colors, s=point_size, alpha=opacity, **kwargs)
        if title == 'auto':
            title = f"Volcano plot of {self.fname.stem}"
        ax.set_title(title, fontsize=title_fontsize)
        ax.set_xlabel(r"$\log_2$(Fold Change)", fontsize=label_fontsize)
        ax.set_ylabel(r'-$\log_{10}$(Adj. P-Value)', fontsize=label_fontsize)
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)

        if log2fc_threshold > 0:
            ax.axvline(log2fc_threshold, linestyle='--', color='black')
            ax.axvline(-log2fc_threshold, linestyle='--', color='black')
        generic.despine(ax)

        plt.show()
        return fig


@readable_name('Count matrix')
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
    _precomputed_metrics = clustering.ClusteringRunner.precomputed_metrics
    _transforms = {True: generic.standard_box_cox, False: generic.standardize}
    __slots__ = {'_is_normalized': 'indicates whether the values in this CountFilter were normalized'}

    def __init__(self, fname: Union[str, Path, tuple], drop_columns: Union[str, List[str]] = None,
                 is_normalized: bool = False, suppress_warnings: bool = False):
        """
        Load a count matrix. A valid count matrix should have one row per gene/genomic feature \
        and one column per condition/RNA library. The contents of the count matrix can be raw or pre-normalized.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded table.
        :type drop_columns: str, list of str, or None (default=None)
        :param is_normalized: indicates whether this count table is pre-normalized. \
        RNAlysis issues a warning when a function meant for normalized tables is applied to a \
        table that was not already normalized.
        :type is_normalized: bool (default=False)
        """
        super().__init__(fname, drop_columns, suppress_warnings)
        self._is_normalized = is_normalized

    def _init_warnings(self):
        super()._init_warnings()
        if len(self._numeric_columns) < len(self.columns):
            warnings.warn(f"The following columns in the CountFilter are not numeric, and will therefore be ignored "
                          f"when running some CountFilter-specific functions: "
                          f"{set(self.columns).difference(self._numeric_columns)}")

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame, name: Union[str, Path], is_normalized: bool = False,
                       suppress_warnings: bool = False) -> 'CountFilter':
        obj = cls.__new__(cls)
        obj.df = df.clone()
        obj.fname = Path(name)
        obj._is_normalized = is_normalized
        if not suppress_warnings:
            obj._init_warnings()
        return obj

    def __copy__(self):
        return type(self).from_dataframe(self.df, self.fname, self.is_normalized, suppress_warnings=True)

    @property
    def is_normalized(self) -> bool:
        return self._is_normalized

    @property
    def _numeric_columns(self) -> list:
        """
        Returns a list of the numeric (int/float) columns in the DataFrame.
        """
        return [self.df.columns[i] for i, dtype in enumerate(self.df.dtypes) if dtype in pl.NUMERIC_DTYPES]

    @property
    def triplicates(self):

        """
        Returns a nested list of the column names in the CountFilter, grouped by alphabetical order into triplicates. \
        For example, if counts.columns is ['A_rep1','A_rep2','A_rep3','B_rep1','B_rep2',_B_rep3'], then \
        counts.triplicates will be  [['A_rep1','A_rep2','A_rep3'],['B_rep1','B_rep2',_B_rep3']]

        """

        multiplier = 3
        numeric_cols = sorted(self._numeric_columns)
        n_cols = len(numeric_cols)
        triplicate = [numeric_cols[i * multiplier:(1 + i) * multiplier] for i in range(n_cols // multiplier)]
        if len(numeric_cols[(n_cols // multiplier) * multiplier::]) > 0:
            triplicate.append(numeric_cols[(n_cols // multiplier) * multiplier::])
            warnings.warn(f'Number of samples {n_cols} is not divisible by 3. '
                          f'Appending the remaining {n_cols % multiplier} as an inncomplete triplicate.')
        return triplicate

    def _diff_exp_assertions(self, design_mat_df: pl.DataFrame):
        assert design_mat_df.shape[0] == self.shape[1], f"The number of items in the design matrix " \
                                                        f"({design_mat_df.shape[0]}) does not match the number of " \
                                                        f"columns in the count matrix ({self.shape[1]})."
        design_mat_samples = sorted(design_mat_df.select(pl.first()).to_series().to_list())
        assert design_mat_samples == sorted(self.columns), f"The sample names in the design matrix do not " \
                                                           f"match the sample names in the count matrix: " \
                                                           f"{design_mat_samples} != " \
                                                           f"{sorted(self.columns)}"
        assert len(design_mat_df.columns) == len(set(design_mat_df.columns)), "The design matrix contains " \
                                                                              "duplicate factor names."
        for factor in design_mat_df.columns:
            assert generic.sanitize_variable_name(
                factor) == factor, f"Invalid factor name '{factor}': contains invalid characters." \
                                   f" \nSuggested alternative name: '{generic.sanitize_variable_name(factor)}'. "

    @readable_name('Run Limma-Voom differential expression')
    def differential_expression_limma_voom(self, design_matrix: Union[str, Path],
                                           comparisons: Iterable[Tuple[str, str, str]],
                                           covariates: Iterable[str] = tuple(), lrt_factors: Iterable[str] = tuple(),
                                           model_factors: Union[Literal['auto'], Iterable[str]] = 'auto',
                                           r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                           output_folder: Union[str, Path, None] = None,
                                           random_effect: Union[str, None] = None, quality_weights: bool = False,
                                           return_design_matrix: bool = False, return_code: bool = False,
                                           return_log: bool = False) -> Tuple['DESeqFilter', ...]:
        """
        Run differential expression analysis on the count matrix using the \
        `Limma-Voom <https://doi.org/10.1186/gb-2014-15-2-r29>`_ pipeline. \
        The count matrix you are analyzing should be normalized (typically to Reads Per Million). \
        The analysis will be based on a design matrix supplied by the user. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix. \
        To run this function, a version of R must be installed.

        :param design_matrix: path to a csv file containing the experiment's design matrix. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix.
        :type design_matrix: str or Path
       :param comparisons: specifies what comparisons to build results tables out of. \
        each individual comparison should be a tuple with exactly three elements: \
        the name of a factor in the design formula, the name of the numerator level for the fold change, \
        and the name of the denominator level for the fold change.
        :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
        :param lrt_factors: optionally, specify factors to be tested using the likelihood ratio test (LRT). \
        If the factors are a continuous variable, you can also specify the number of polynomial degree to fit.
        :type lrt_factors: Iterable of tuple(factor, polynomial_degree) or None (default=None)
        :param covariates: optionally, specify a list of continuous covariates to include in the analysis. \
        The covariates should be column names in the design matrix. \
        The reported fold change values correspond to the expected fold change \
        for every increase of 1 unit in the covariate.
        :param model_factors: optionally, specify a list of factors to include in the differential expression model. \
        If 'auto', all factors in the design matrix will be included.
        :type model_factors: Iterable of factor names or 'auto' (default='auto')
        :param r_installation_folder: Path to the installation folder of R. For example: \
        'C:/Program Files/R/R-4.2.1'
        :type r_installation_folder: str, Path, or 'auto' (default='auto')
        :param output_folder: Path to a folder in which the analysis results, \
        as well as the log files and R script used to generate them, will be saved. \
        if output_folder is None, the results will not be saved to a specified directory.
        :type output_folder: str, Path, or None
        :param random_effect: optionally, specify a single factor to model as a random effect. \
        This is useful when your experimental design is nested. limma-voom can only treat one factor as a random effect.
        :type random_effect: str or None
        :param quality_weights: if True, the analysis will use estimate sample-specific quality weights \
        using the 'arrayWeights' function im limma. This is useful when lower quality samples are present in the data.
        :type quality_weights: bool (default=False)
        :param return_design_matrix: if True, the function will return the sanitized design matrix used in the analysis.
        :type return_design_matrix: bool (default=False)
        :param return_code: if True, the function will return the path to the R script \
        used to generate the analysis results.
        :type return_code: bool (default=False)
        :param return_log: if True, the function will return the path to the analysis logfile, \
        which includes session info.
        :type return_log: bool (default=False)
        :return: a tuple of DESeqFilter objects, one for each comparison
        """
        if output_folder is not None:
            output_folder = Path(output_folder)
            assert output_folder.exists(), 'Output folder does not exist!'

        self._validate_is_normalized(expect_normalized=False)
        data_path = io.get_todays_cache_dir().joinpath(f'{self.fname.name}.csv')
        design_mat_path = None
        i = 0
        while design_mat_path is None or design_mat_path.exists():
            design_mat_path = io.get_todays_cache_dir().joinpath(f'design_mat_{i}.csv')
            i += 1

        io.save_table(self.df.select(pl.first()).with_columns(
            self.df.select(pl.col(self._numeric_columns))), data_path)
        # use RNAlysis to automatically detect file delimiter type, then export it as a CSV file.
        design_mat_df = io.load_table(design_matrix)
        self._diff_exp_assertions(design_mat_df)
        samples = design_mat_df.select(pl.first()).to_series().to_list()
        design_mat_df = design_mat_df.lazy().with_columns(pl.Series([self.df.columns.index(i) for i in samples])).sort(
            pl.last()).drop(cs.last()).collect()
        assert parsing.data_to_list(design_mat_df.select(pl.first())) == list(
            self.columns), "The sample names in the design matrix do not match!"
        io.save_table(design_mat_df, design_mat_path)

        r_output_dir = differential_expression.LimmaVoomRunner(data_path, design_mat_path, comparisons, covariates,
                                                               lrt_factors, model_factors, r_installation_folder,
                                                               random_effect, quality_weights).run()
        code_path = None
        log_path = None
        outputs = []
        for item in r_output_dir.iterdir():
            if not item.is_file():
                continue
            if item.suffix == '.csv':
                outputs.append(DESeqFilter(item, log2fc_col='logFC', padj_col='adj.P.Val', pval_col='P.Value'))
            elif item.suffix == '.R':
                code_path = item
            elif item.suffix == '.log':
                log_path = item
            if output_folder is not None:
                with open(item) as infile, open(output_folder.joinpath(item.name), 'w') as outfile:
                    outfile.write(infile.read())
        return_val = parsing.data_to_tuple(outputs)
        if return_design_matrix:
            return_val = [return_val, design_mat_df]
        if return_code:
            if code_path is None:
                warnings.warn("No R script was generated during the analysis")
            else:
                return_val = [return_val, code_path]
        if return_log:
            if log_path is None:
                warnings.warn("No log file was generated during the analysis")
            else:
                return_val = [return_val, log_path]
        return return_val

    @readable_name('Run Limma-Voom differential expression (simplified mode)')
    def differential_expression_limma_voom_simplified(self, design_matrix: Union[str, Path],
                                                      comparisons: Iterable[Tuple[str, str, str]],
                                                      r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                                      output_folder: Union[str, Path, None] = None,
                                                      random_effect: Union[str, None] = None,
                                                      quality_weights: bool = False,
                                                      return_design_matrix: bool = False, return_code: bool = False,
                                                      return_log: bool = False) -> Tuple['DESeqFilter', ...]:
        """
       Run differential expression analysis on the count matrix using the \
       `Limma-Voom <https://doi.org/10.1186/gb-2014-15-2-r29>`_ pipeline. \
       The simplified mode supports only pairwise comparisons. \
       The count matrix you are analyzing should be normalized (typically to Reads Per Million). \
       The analysis will be based on a design matrix supplied by the user. \
       The design matrix should contain at least two columns: the first column contains all the sample names, \
       and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
       (see the User Guide and Tutorial for a complete example). \
       The analysis formula will contain all the factors in the design matrix. \
       To run this function, a version of R must be installed.

       :param design_matrix: path to a csv file containing the experiment's design matrix. \
       The design matrix should contain at least two columns: the first column contains all the sample names, \
       and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
       (see the User Guide and Tutorial for a complete example). \
       The analysis formula will contain all the factors in the design matrix.
       :type design_matrix: str or Path
      :param comparisons: specifies what comparisons to build results tables out of. \
       each individual comparison should be a tuple with exactly three elements: \
       the name of a factor in the design formula, the name of the numerator level for the fold change, \
       and the name of the denominator level for the fold change.
       :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
       :param r_installation_folder: Path to the installation folder of R. For example: \
       'C:/Program Files/R/R-4.2.1'
       :type r_installation_folder: str, Path, or 'auto' (default='auto')
       :param output_folder: Path to a folder in which the analysis results, \
       as well as the log files and R script used to generate them, will be saved. \
       if output_folder is None, the results will not be saved to a specified directory.
       :type output_folder: str, Path, or None
       :param random_effect: optionally, specify a single factor to model as a random effect. \
       This is useful when your experimental design is nested. limma-voom can only treat one factor as a random effect.
       :type random_effect: str or None
       :param quality_weights: if True, the analysis will use estimate sample-specific quality weights \
        using the 'arrayWeights' function im limma. This is useful when lower quality samples are present in the data.
        :type quality_weights: bool (default=False)
        :param return_design_matrix: if True, the function will return the sanitized design matrix used in the analysis.
        :type return_design_matrix: bool (default=False)
        :param return_code: if True, the function will return the R script used to generate the analysis results.
        :type return_code: bool (default=False)
       :return: a tuple of DESeqFilter objects, one for each comparison
        """
        return self.differential_expression_limma_voom(design_matrix, comparisons,
                                                       r_installation_folder=r_installation_folder,
                                                       output_folder=output_folder, random_effect=random_effect,
                                                       quality_weights=quality_weights,
                                                       return_design_matrix=return_design_matrix,
                                                       return_code=return_code, return_log=return_log)

    @readable_name('Run DESeq2 differential expression')
    def differential_expression_deseq2(self, design_matrix: Union[str, Path],
                                       comparisons: Iterable[Tuple[str, str, str]],
                                       covariates: Iterable[str] = tuple(),
                                       lrt_factors: Iterable[str] = tuple(),
                                       model_factors: Union[Literal['auto'], Iterable[str]] = 'auto',
                                       r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                       output_folder: Union[str, Path, None] = None, return_design_matrix: bool = False,
                                       scaling_factors: Union[str, Path, None] = None,
                                       cooks_cutoff: bool = True,
                                       return_code: bool = False, return_log: bool = False
                                       ) -> Tuple['DESeqFilter', ...]:
        """
        Run differential expression analysis on the count matrix using the \
        `DESeq2 <https://doi.org/10.1186/s13059-014-0550-8>`_ algorithm. \
        The count matrix you are analyzing should be unnormalized (meaning, raw read counts). \
        The analysis will be based on a design matrix supplied by the user. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix. \
        To run this function, a version of R must be installed.

        :param design_matrix: path to a csv file containing the experiment's design matrix. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix.
        :type design_matrix: str or Path
        :param comparisons: specifies what comparisons to build results tables out of. \
        each individual comparison should be a tuple with exactly three elements: \
        the name of a factor in the design formula, the name of the numerator level for the fold change, \
        and the name of the denominator level for the fold change.
        :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
        :param lrt_factors: optionally, specify factors to be tested using the likelihood ratio test (LRT). \
        If the factors are a continuous variable, you can also specify the number of polynomial degree to fit.
        :type lrt_factors: Iterable of factor names (default=tuple())
        :param covariates: optionally, specify a list of continuous covariates to include in the analysis. \
        The covariates should be column names in the design matrix. \
        The reported fold change values correspond to the expected fold change \
        for every increase of 1 unit in the covariate.
        :type covariates: Iterable of covariate names (default=tuple())
        :param model_factors: optionally, specify a list of factors to include in the differential expression model. \
        If 'auto', all factors in the design matrix will be included.
        :type model_factors: Iterable of factor names or 'auto' (default='auto')
        :param r_installation_folder: Path to the installation folder of R. For example: \
        'C:/Program Files/R/R-4.2.1'
        :type r_installation_folder: str, Path, or 'auto' (default='auto')
        :param output_folder: Path to a folder in which the analysis results, \
        as well as the log files and R script used to generate them, will be saved. \
        if output_folder is None, the results will not be saved to a specified directory.
        :type output_folder: str, Path, or None
        :param return_design_matrix: if True, the function will return the sanitized design matrix used in the analysis.
        :type return_design_matrix: bool (default=False)
        :param return_code: if True, the function will return the path to the R script \
        used to generate the analysis results.
        :type return_code: bool (default=False)
        :param return_log: if True, the function will return the path to the analysis logfile, \
        which includes session info.
        :type return_log: bool (default=False)
        :return: a tuple of DESeqFilter objects, one for each comparison
        """
        if output_folder is not None:
            output_folder = Path(output_folder)
            assert output_folder.exists(), 'Output folder does not exist!'

        self._validate_is_normalized(expect_normalized=False)
        data_path = io.get_todays_cache_dir().joinpath(f'{self.fname.stem}.csv')
        design_mat_path = None
        i = 0
        while design_mat_path is None or design_mat_path.exists():
            design_mat_path = io.get_todays_cache_dir().joinpath(f'design_mat_{i}.csv')
            i += 1

        counts = self.df.lazy().select(pl.first()).with_columns(
            self.df.lazy().select(pl.col(self._numeric_columns)).with_columns(cs.float().round()).collect()).sort(
            pl.first()).collect()
        io.save_table(counts, data_path)
        # use RNAlysis to automatically detect file delimiter type, then export it as a CSV file.
        design_mat_df = io.load_table(design_matrix)
        self._diff_exp_assertions(design_mat_df)
        samples = design_mat_df.select(pl.first()).to_series().to_list()
        design_mat_df = design_mat_df.lazy().with_columns(pl.Series([self.df.columns.index(i) for i in samples])).sort(
            pl.last()).drop(cs.last()).collect()
        assert parsing.data_to_list(design_mat_df.select(pl.first())) == list(
            self.columns), "The sample names in the design matrix do not match!"
        io.save_table(design_mat_df, design_mat_path)

        if scaling_factors is None:
            scale_factor_path = None
            scale_factor_ndims = 1
        else:
            scaling_factors_df = io.load_table(scaling_factors)
            if len(scaling_factors_df) == 1:
                scale_factor_ndims = 1
                scaling_factors_df = scaling_factors_df.select(self.df.columns[1:])
                assert parsing.data_to_list(scaling_factors_df.select(pl.first())) == list(self.columns), \
                    "The sample names in the scaling factors table do not match the sample names in the count matrix!"
            else:
                scale_factor_ndims = 2
                scaling_factors_df = scaling_factors_df.lazy().with_columns(pl.col(counts.columns[1:]))
                scaling_factors_df = scaling_factors_df.sort(pl.first()).collect()
                assert sorted(scaling_factors_df.columns) == sorted(self.columns), \
                    "The sample names in the scaling factors table do not match the sample names in the count matrix!"
                assert scaling_factors_df.select(pl.first()).equals(counts.select(pl.first())), \
                    "The gene names in the scaling factors table do not match the gene names in the count matrix!"
            scale_factor_path = None
            i = 0
            while scale_factor_path is None or scale_factor_path.exists():
                scale_factor_path = io.get_todays_cache_dir().joinpath(f'scaling_factors_{i}.csv')
                i += 1
            io.save_table(scaling_factors_df, scale_factor_path)

        r_output_dir = differential_expression.DESeqRunner(data_path, design_mat_path, comparisons, covariates,
                                                           lrt_factors, model_factors, r_installation_folder,
                                                           scale_factor_path, cooks_cutoff, scale_factor_ndims).run()
        outputs = []
        code_path = None
        log_path = None
        for item in r_output_dir.iterdir():
            if not item.is_file():
                continue
            if item.suffix == '.csv':
                outputs.append(DESeqFilter(item))
            if item.suffix == '.R':
                code_path = item
            elif item.suffix == '.log':
                log_path = item
            if output_folder is not None:
                with open(item) as infile, open(output_folder.joinpath(item.name), 'w') as outfile:
                    outfile.write(infile.read())
        return_val = parsing.data_to_tuple(outputs)
        if return_design_matrix:
            return_val = [return_val, design_mat_df]
        if return_code:
            if code_path is None:
                warnings.warn("No R script was generated during the analysis")
            else:
                return_val = [return_val, code_path]
        if return_log:
            if log_path is None:
                warnings.warn("No log file was generated during the analysis")
            else:
                return_val = [return_val, log_path]
        return return_val

    @readable_name('Run DESeq2 differential expression (simplified mode)')
    def differential_expression_deseq2_simplified(self, design_matrix: Union[str, Path],
                                                  comparisons: Iterable[Tuple[str, str, str]],
                                                  r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                                  output_folder: Union[str, Path, None] = None,
                                                  return_design_matrix: bool = False, return_code: bool = False,
                                                  return_log: bool = False) -> Tuple['DESeqFilter', ...]:
        """
        Run differential expression analysis on the count matrix using the \
        `DESeq2 <https://doi.org/10.1186/s13059-014-0550-8>`_ algorithm. \
        The simplified mode supports only pairwise comparisons. \
        The count matrix you are analyzing should be unnormalized (meaning, raw read counts). \
        The analysis will be based on a design matrix supplied by the user. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix. \
        To run this function, a version of R must be installed.

        :param design_matrix: path to a csv file containing the experiment's design matrix. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix.
        :type design_matrix: str or Path
        :param comparisons: specifies what comparisons to build results tables out of. \
        each individual comparison should be a tuple with exactly three elements: \
        the name of a factor in the design formula, the name of the numerator level for the fold change, \
        and the name of the denominator level for the fold change.
        :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
        :param r_installation_folder: Path to the installation folder of R. For example: \
        'C:/Program Files/R/R-4.2.1'
        :type r_installation_folder: str, Path, or 'auto' (default='auto')
        :param output_folder: Path to a folder in which the analysis results, \
        as well as the log files and R script used to generate them, will be saved. \
        if output_folder is None, the results will not be saved to a specified directory.
        :type output_folder: str, Path, or None
        :param return_design_matrix: if True, the function will return the sanitized design matrix used in the analysis.
        :type return_design_matrix: bool (default=False)
        :param return_code: if True, the function will return the R script used to generate the analysis results.
        :type return_code: bool (default=False)
        :return: a tuple of DESeqFilter objects, one for each comparison
        """
        return self.differential_expression_deseq2(design_matrix, comparisons,
                                                   r_installation_folder=r_installation_folder,
                                                   output_folder=output_folder,
                                                   return_design_matrix=return_design_matrix, return_code=return_code,
                                                   return_log=return_log)

    @readable_name('Calculate fold change')
    def fold_change(self, numerator: param_typing.ColumnNames, denominator: param_typing.ColumnNames,
                    numer_name: str = 'default', denom_name: str = 'default') -> 'FoldChangeFilter':
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

        numerator = parsing.data_to_list(numerator)
        denominator = parsing.data_to_list(denominator)
        assert validation.isinstanceiter(numerator, str), "numerator must be str or a list of str!"
        assert validation.isinstanceiter(denominator, str), "denominator must be str or a list of str"

        numeric_cols = self._numeric_columns
        for num in numerator:
            assert num in self.columns, f"'{num}' is not a column in the CountFilter object!"
            assert num in numeric_cols, f"Invalid dtype for column '{num}': {self.df.dtypes[num]}"
        for den in denominator:
            assert den in self.columns, f"'{den}' is not a column in the CountFilter object!"
            assert den in numeric_cols, f"Invalid dtype for column '{den}': {self.df.dtypes[den]}"

        fc_df = self.df.select(pl.first()).with_columns(
            ((self.df.select(numerator).mean_horizontal() + 1) / (
                self.df.select(denominator).mean_horizontal() + 1)).alias('Fold Change'))
        new_fname = Path(f"{str(self.fname.parent)}/{self.fname.stem}'_fold_change_'"
                         f"{numer_name}_over_{denom_name}_{self.fname.suffix}")

        fcfilt = FoldChangeFilter.from_dataframe(fc_df, new_fname, numerator_name=numer_name,
                                                 denominator_name=denom_name)

        return fcfilt

    @readable_name('Pair-plot')
    def pairplot(self, samples: Union[param_typing.GroupedColumns, Literal['all']] = 'all',
                 log2: bool = True, show_corr: bool = True, title: Union[str, Literal['auto']] = 'auto',
                 title_fontsize: float = 30,
                 label_fontsize: float = 16, tick_fontsize: float = 12) -> plt.Figure:

        """
        Plot pairwise relationships in the dataset. \
        Can plot both single samples and average multiple replicates. \
        For more information see the documentation of seaborn.pairplot.

        :type samples: 'all', list, or nested list.
        :param samples: A list of the sample names and/or grouped sample names to be included in the pairplot. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :type log2: bool (default=True)
        :param log2: if True, the pairplot will be calculated with log2 of the DataFrame (pseudocount+1 added), \
        and not with the raw data. If False, the pairplot will be calculated with the raw data.
        :type show_corr: bool (default=True)
        :param show_corr: if True, shows the Spearman correlation coefficient (R) between each pair of samples/groups.
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
         :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :return: A matplotlib Figure.

        .. figure:: /figures/pairplot.png
           :align:   center
           :scale: 40 %

           Example plot of pairplot()

        """
        if samples == 'all':
            sample_df = self.df
        else:
            sample_df = self._avg_subsamples(samples)

        if log2:
            # sample_df = np.log2(sample_df + 1)
            sample_df = sample_df.with_columns(pl.col(sample_df.columns[1:]).add(1).log(2))

        pairplt = sns.pairplot(sample_df.to_pandas(), corner=True,
                               plot_kws=dict(alpha=0.25, edgecolors=(0.1, 0.5, 0.15),
                                             facecolors=(0.1, 0.5, 0.15), s=3.5),
                               diag_kws=dict(edgecolor=(0, 0, 0), facecolor=(0.1, 0.5, 0.15)))

        if title == 'auto':
            title = 'Pairplot' + log2 * ' (logarithmic scale)'
        plt.suptitle(title, fontsize=title_fontsize)

        for i, row in enumerate(pairplt.axes):
            for j, ax in enumerate(row[0:i + 1]):
                ax.set_xlabel(ax.get_xlabel(), fontsize=label_fontsize)
                ax.set_ylabel(ax.get_ylabel(), fontsize=label_fontsize)
                ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)

                if i == j:
                    continue
                if show_corr:
                    spearman_corr = spearmanr(sample_df[:, [i, j]])[0]
                    ax.text(0.05, 0.9, f"Spearman \u03C1={spearman_corr:.2f}", transform=ax.transAxes)

        pairplt.figure.set_size_inches(9, 9)
        try:
            pairplt.figure.set_layout_engine("compressed")
        except AttributeError:
            pairplt.figure.set_constrained_layout(True)
        plt.show()
        return pairplt.figure

    @readable_name('Average the expression values of replicate columns')
    def average_replicate_samples(self, sample_grouping: param_typing.GroupedColumns,
                                  new_column_names: Union[Literal['auto'], List[str]] = 'auto',
                                  function: Literal['mean', 'median', 'geometric_mean'] = 'mean',
                                  inplace: bool = True) -> 'CountFilter':
        """
        Average the expression values of gene expression for each group of replicate samples. \
        Each group of samples (e.g. biological/technical replicates)

        :param sample_grouping: grouping of the samples into conditions. \
        Each grouping should containg all replicates of the same condition. \
        Each condition will be averaged separately.
        :type sample_grouping: nested list of column names
        :param new_column_names: names to be given to the columns in the new count matrix. \
        Each new name should match a group of samples to be averaged. \
        If `new_column_names`='auto', names will be generated automatically.
        :type new_column_names: list of str or 'auto' (default='auto'
        :param function: the function which will be used to average the values within each group.
        :type function: 'mean', 'median', or 'geometric_mean' (default='mean')
        :type inplace: bool (default=True)
        :param inplace: If True (default), averaging will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        """
        suffix = f'_{function}'
        avg_df = self._avg_subsamples(sample_grouping, function, new_column_names)
        return self._inplace(avg_df, False, inplace, suffix, 'transform')

    def _avg_subsamples(self, sample_grouping: param_typing.GroupedColumns,
                        function: Literal['mean', 'median', 'geometric_mean'] = 'mean',
                        new_column_names: Union[Literal['auto'], Literal['display'], List[str]] = 'display'):

        """
        Avarages subsamples/replicates according to the specified sample list. \
        Every member in the sample list should be either a name of a single sample (str), \
        or a list of multiple sample names to be averaged (list).

        :param sample_grouping: A list of the sample names and/or grouped sample names passed by the user. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6'] \
        and the resulting output will be a DataFrame containing the following columns: \
        ['SAMPLE1', 'SAMPLE2', 'SAMPLE3', 'SAMPLE6']
        :return: a pandas DataFrame containing samples/averaged subsamples according to the specified sample_list.

        """
        assert isinstance(function, str), "'function' must be a string!"
        function = function.lower()
        assert function in {'mean', 'median', 'geometric_mean'}, \
            "'function' must be 'mean', 'median', or 'geometric_mean'!"

        if new_column_names in ('auto', 'display'):
            new_column_names = parsing.make_group_names(sample_grouping, new_column_names)
        else:
            assert validation.isiterable(new_column_names) and validation.isinstanceiter(new_column_names, str), \
                "'new_column_names' must be either 'auto' or a list of strings!"
            assert len(new_column_names) == len(sample_grouping), \
                f"The number of new column names {len(new_column_names)} " \
                f"does not match the number of sample groups {len(sample_grouping)}!"

        averaged_df = self.df.select(pl.first())

        for group, new_name in zip(sample_grouping, new_column_names):
            if isinstance(group, str):
                assert group in self.columns, f"Column '{group}' does not exist in the original table!"
                averaged_df = averaged_df.with_columns(self.df[group].alias(new_name))
            elif isinstance(group, (list, tuple, set)):
                for item in group:
                    assert item in self.columns, f"Column '{item}' does not exist in the original table!"
                if function == 'mean':
                    averaged_df = averaged_df.with_columns(
                        self.df.select(pl.col(group)).mean_horizontal().alias(new_name))
                elif function == 'median':
                    averaged_df = averaged_df.with_columns(
                        self.df.select(pl.col(group)).median_horizontal().alias(new_name))
                else:
                    averaged_df = averaged_df.with_columns(
                        self.df.select(pl.col(group).apply(lambda x: gmean(x)).alias(new_name)))

            else:
                raise TypeError(f"'sample_list' cannot contain objects of type {type(group)}.")

        return averaged_df

    def _validate_is_normalized(self, expect_normalized: bool = True):
        if not self.is_normalized and expect_normalized:
            warnings.warn("This function is meant for normalized values, and your count matrix may not be normalized. ")
        elif self.is_normalized and not expect_normalized:
            warnings.warn(
                "This function is meant for raw, unnormalize counts, and your count matrix appears to be normalized. ")

    def _norm_scaling_factors(self, scaling_factors: pl.DataFrame):
        numeric_cols = self._numeric_columns
        new_df = pl.DataFrame().lazy()

        if scaling_factors.shape[0] == 1:
            assert scaling_factors.shape[1] == len(numeric_cols), \
                f"Number of scaling factors ({scaling_factors.shape[1]}) does not match " \
                f"number of numeric columns in your data table ({len(numeric_cols)})!"

            for column in self.df.columns:
                if column in numeric_cols:
                    norm_factor = scaling_factors[column]
                    new_df = new_df.with_columns((self.df.select(pl.col(column).truediv(norm_factor))))
                else:
                    new_df = new_df.with_columns(self.df[column].alias(column))
        else:
            assert scaling_factors.shape[0] >= self.shape[0] and scaling_factors.shape[1] == len(numeric_cols) + 1, \
                f"Dimensions of scaling factors table ({scaling_factors.shape}) does not match the " \
                f"dimensions of your data table ({(self.shape[0], len(numeric_cols))} - numeric columns only)!"
            for column in self.df.columns:
                if column in numeric_cols:
                    merged = self.df.select(cs.first() | cs.by_name(column)).join(
                        scaling_factors.select(cs.first() | cs.by_name(column)), left_on=self.df.columns[0],
                        right_on=scaling_factors.columns[0], how='left')
                    merged_div = merged.with_columns((pl.nth(-2).truediv(pl.nth(-1))).alias('div'))
                    new_df = new_df.with_columns((merged_div.select(pl.col('div').alias(column))))
                else:
                    new_df = new_df.with_columns(self.df[column].alias(column))

        return new_df.collect()

    @readable_name('Normalize to reads-per-million (RPM) - HTSeq-count output')
    def normalize_to_rpm_htseqcount(self, special_counter_fname: Union[str, Path], inplace: bool = True,
                                    return_scaling_factors: bool = False):

        """
        Normalizes the count matrix to Reads Per Million (RPM). \
        Uses a table of feature counts (ambiguous, no feature, not aligned, etc) from HTSeq-count's output. \
        Divides each column in the CountFilter object by (total reads + ambiguous + no feature)*10^-6 .

        :param special_counter_fname: the .csv file which contains feature information about the RNA library \
        (ambiguous, no feature, not aligned, etc).
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_to_rpm_htseqcount("tests/test_files/uncounted.csv")
           Normalized 22 features. Normalized inplace.

        """
        self._validate_is_normalized(expect_normalized=False)

        suffix = '_normtoRPMhtseqcount'
        scaling_factors = {}

        if isinstance(special_counter_fname, (str, Path)):
            features = io.load_table(special_counter_fname)
        elif isinstance(special_counter_fname, pl.DataFrame):
            features = special_counter_fname
        else:
            raise TypeError("Invalid type for 'special_counter_fname'!")

        for column in self.df.columns:
            if column in self._numeric_columns:
                norm_factor = (self.df.select(pl.col(column)).sum() + features.filter(
                    pl.first().is_in(['__no_feature', '__alignment_not_unique', '__ambiguous'])).select(
                    pl.col(column)).sum()) / (10 ** 6)
                scaling_factors[column] = norm_factor

        scaling_factors = pl.DataFrame(scaling_factors)
        new_df = self._norm_scaling_factors(scaling_factors)

        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @readable_name('Normalize to reads-per-million (RPM)')
    def normalize_to_rpm(self, inplace: bool = True, return_scaling_factors: bool = False):
        """
        Normalizes the count matrix to Reads Per Million (RPM). \
        Divides each column in the count matrix by (total reads)*10^-6 .

        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_to_rpm()
           Normalized 22 features. Normalized inplace.

        """
        self._validate_is_normalized(expect_normalized=False)

        suffix = '_normtoRPM'
        scaling_factors = self.df.select(pl.col(self._numeric_columns)).sum() / (10 ** 6)
        new_df = self._norm_scaling_factors(scaling_factors)

        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @staticmethod
    def _gene_len_kbp(gtf_file: Union[str, Path], feature_type: Literal['gene', 'transcript'],
                      method: Literal[LEGAL_GENE_LENGTH_METHODS]):
        gene_lengths_dict = {key: val / 1000 for key, val in
                             genome_annotation.get_genomic_feature_lengths(gtf_file, feature_type, method).items()}
        gene_lengths_kbp = pl.from_records([list(gene_lengths_dict.keys()), list(gene_lengths_dict.values())],
                                           schema=['gene', 'length'])
        return gene_lengths_kbp

    @readable_name('Normalize to reads-per-kilobase-million (RPKM)')
    def normalize_to_rpkm(self, gtf_file: Union[str, Path], feature_type: Literal['gene', 'transcript'] = 'gene',
                          method: Literal[LEGAL_GENE_LENGTH_METHODS] = 'mean', inplace: bool = True,
                          return_scaling_factors: bool = False):
        """
        Normalizes the count matrix to Reads Per Kilobase Million (RPKM). \
        Divides each column in the count matrix by (total reads)*(gene length / 1000)*10^-6. \


       :param gtf_file: Path to a GTF/GFF3 annotation file. This file will be used to determine the length of each \
        gene/transcript. The gene/transcript names in this annotation file should match the ones in count matrix.
        :type gtf_file: str or Path
        :param feature_type: the type of features in your count matrix. if feature_type is 'transcript', \
        lengths will be calculated per-transcript, and the 'method' parameter is ignored. \
        Otherwise, lengths will be aggregated per gene according to the method specified in the 'method' parameter.
        :type feature_type: 'gene' or 'transcript' (default='gene')
        :param method: if feature_type='gene', this determines the aggregation method to calculate gene lengths. \
        'mean', 'median', 'min', and 'max' will calculate the \
        mean/median/min/max of all transcripts' lengths of the given gene. \
        'geometric_mean' will calculate the goemetric mean of all transcripts' lengths of the given gene. \
        'merged_exons' will calculate the total lengths of all exons of a gene across all of its transcripts, \
        while counting overlapping exons/regions exactly once.
        :type method: 'mean', 'median', 'min', 'max', 'geometric_mean', or 'merged_exons' (deafult='mean')
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.
        """
        self._validate_is_normalized(expect_normalized=False)

        suffix = f"_normtoRPKM{method.replace('_', '')}"
        numeric_cols = self._numeric_columns

        lib_coefs = self.df.drop(cs.first()).sum() * (10 ** -6)
        scaling_factors = self.df.select(pl.first()).with_columns(pl.concat([lib_coefs] * len(self)))

        gene_lengths_kbp = self._gene_len_kbp(gtf_file, feature_type, method)
        scaling_factors = scaling_factors.join(gene_lengths_kbp, left_on=self.df.columns[0], right_on='gene',
                                               how='left')
        for this_col in numeric_cols:
            scaling_factors = scaling_factors.with_columns(
                scaling_factors.select(pl.col(this_col).mul(pl.col('length')).alias(this_col)))
        scaling_factors = scaling_factors.drop('length')

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @readable_name('Normalize to transcripts-per-million (TPM)')
    def normalize_to_tpm(self, gtf_file: Union[str, Path], feature_type: Literal['gene', 'transcript'] = 'gene',
                         method: Literal[LEGAL_GENE_LENGTH_METHODS] = 'mean', inplace: bool = True,
                         return_scaling_factors: bool = False):
        """
        Normalizes the count matrix to Transcripts Per Million (TPM). \
        First, normalizes each gene to Reads Per Kilobase (RPK) \
        by dividing each gene in the count matrix by its length in Kbp (gene length / 1000). \
        Then, divides each column in the RPK matrix by (total RPK in column)*10^-6. \
        This calculation is similar to that of Reads Per Kilobase Million (RPKM), but in the opposite order: \
        the "per million" normalization factors are calculated **after** normalizing to gene lengths, not before.

       :param gtf_file: Path to a GTF/GFF3 annotation file. This file will be used to determine the length of each \
        gene/transcript. The gene/transcript names in this annotation file should match the ones in count matrix.
        :type gtf_file: str or Path
        :param feature_type: the type of features in your count matrix. if feature_type is 'transcript', \
        lengths will be calculated per-transcript, and the 'method' parameter is ignored. \
        Otherwise, lengths will be aggregated per gene according to the method specified in the 'method' parameter.
        :type feature_type: 'gene' or 'transcript' (default='gene')
        :param method: if feature_type='gene', this determines the aggregation method to calculate gene lengths. \
        'mean', 'median', 'min', and 'max' will calculate the \
        mean/median/min/max of all transcripts' lengths of the given gene. \
        'geometric_mean' will calculate the goemetric mean of all transcripts' lengths of the given gene. \
        'merged_exons' will calculate the total lengths of all exons of a gene across all of its transcripts, \
        while counting overlapping exons/regions exactly once.
        :type method: 'mean', 'median', 'min', 'max', 'geometric_mean', or 'merged_exons' (deafult='mean')
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.
        """
        self._validate_is_normalized(expect_normalized=False)

        suffix = f"_normtoTPM{method.replace('_', '')}"
        gene_lengths_kbp = self._gene_len_kbp(gtf_file, feature_type, method)
        numeric_cols = self._numeric_columns

        tmp_df = self.df.join(gene_lengths_kbp, left_on=self.df.columns[0], right_on='gene', how='left')
        for this_col in numeric_cols:
            tmp_df = tmp_df.with_columns(
                tmp_df.select(pl.col(this_col).truediv(pl.nth(-1)).alias(this_col)))
        tmp_df = tmp_df.drop('length')

        lib_coefs = tmp_df.drop(cs.first()).sum() / (10 ** 6)
        scaling_factors = tmp_df.select(pl.first()).with_columns(
            pl.DataFrame({key: gene_lengths_kbp.select(pl.nth(-1)) for key in numeric_cols}))

        for this_col in numeric_cols:
            scaling_factors = scaling_factors.with_columns(
                scaling_factors.select(pl.col(this_col).mul(lib_coefs.select(pl.col(this_col))).alias(this_col)))

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @readable_name('Normalize with the Quantile method')
    def normalize_to_quantile(self, quantile: param_typing.Fraction = 0.75, inplace: bool = True,
                              return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the quantile method, generalized from \
        `Bullard et al 2010 <https://doi.org/10.1186/1471-2105-11-94>`_. \
        This is the default normalization method used by R's Limma. \
        To calculate the Quantile Method scaling factors, you first calculate the given quantile of gene expression \
        within each sample, excluding genes that have 0 reads in all samples. \
        You then divide those quantile values by the total number of reads in each sample, \
        which yields the scaling factors for each sample.

        :param quantile: the quantile from which scaling factors will be calculated.
        :type quantile: float between 0 and 1 (default=0.75)
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_to_quantile(0.75)
           Normalized 22 features. Normalized inplace.
        """
        self._validate_is_normalized(expect_normalized=False)

        suffix = f'_normto{int(quantile * 100)}quantile'
        data = self.df.select(pl.col(self._numeric_columns))
        expressed_genes = data.filter(data.sum_horizontal() != 0)
        quantiles = expressed_genes.quantile(quantile, interpolation='linear')
        if (quantiles.min_horizontal() == 0).any():
            warnings.warn("One or more quantiles are zero")

        scaling_factors = quantiles / quantiles.mean_horizontal()
        new_df = self._norm_scaling_factors(scaling_factors)

        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @readable_name('Normalize with the "Trimmed Mean of M-values" (TMM) method')
    def normalize_tmm(self, log_ratio_trim: float = 0.3, sum_trim: float = 0.05,
                      a_cutoff: Union[float, None] = -1 * 10 ** 10,
                      ref_column: Union[Literal['auto'], param_typing.ColumnName] = 'auto',
                      inplace: bool = True, return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the 'trimmed mean of M values' (TMM) method \
        `(Robinson and Oshlack 2010) <https://doi.org/10.1186/gb-2010-11-3-r25>`_. \
        This is the default normalization method used by R's edgeR. \
        To calculate the Trimmed Mean of M Values scaling factors, you first calculate the M-values of each gene \
        between each sample and the reference sample (log2 of each sample Minus log2 of the reference sample), \
        and the A-values of each gene between each sample and the reference sample \
        (log2 of each sample Added to log2 of the reference sample). \
        You then trim out genes with extreme values that are likely to be differentially expressed or non-indicative, \
        by trimming the top and bottom X% of M-values, the top and bottom Y% of A-values, all A-values which are \
        smaller than the specified cutuff, and all genes with 0 reads (to avoid log2 values of inf or -inf). \
        Next, a weighted mean is calculated on the filtered M-values, with the weights being an inverse of \
        an approximation of variance of each gene, which gives out the scaling factors for each sample. \
        Finally, the scaling factors are adjusted, for symmetry, so that they multiply to 1.

        :param log_ratio_trim: the fraction of M-values that should be trimmed from each direction (top and bottom X%).
        :type log_ratio_trim:  float between 0 and 0.5 (default=0.3)
        :param sum_trim: the fraction of A-values that should be trimmed from each direction (top and bottom Y%).
        :type sum_trim:  float between 0 and 0.5 (default=0.05)
        :param a_cutoff: a lower bound on the A-values that should be included in the trimmed mean. \
        If set to None, no lower bound will be used.
        :type a_cutoff: float or None (default = -1e10)
        :param ref_column: the column to be used as reference for normalization. If 'auto', \
        then the reference column will be chosen automatically to be \
        the column whose upper quartile is closest to the mean upper quartile.
        :type ref_column: name of a column or 'auto' (default='auto')
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_tmm()
           Normalized 22 features. Normalized inplace.
        """
        self._validate_is_normalized(expect_normalized=False)
        assert 0 <= log_ratio_trim < 0.5, "'log_ratio_trim' must be a value in the range 0 <= log_ratio_trim < 5"
        assert 0 <= sum_trim < 0.5, "'sum_trim' must be a value in the range 0 <= sum_trim < 5"
        suffix = '_normTMM'
        ref_column = self._get_ma_ref_column(ref_column)
        columns = self._numeric_columns
        m_data, a_data, weights = self._calculate_ma(ref_column, columns)
        scaling_factors = {}
        for col in columns:
            # remove genes with 0 expression in the current column
            zero_genes = parsing.data_to_set(
                self.df.lazy().select(pl.first(), col).filter(pl.last() == 0).select(pl.first()).collect())
            # apply cutoff to A values
            a_post_cutoff = a_data.lazy().select(pl.first(), pl.col(col)).filter(pl.last().is_not_nan())
            if a_cutoff is not None:
                a_post_cutoff = a_post_cutoff.filter(pl.col(col) > a_cutoff)
            a_post_cutoff = a_post_cutoff.filter(~pl.first().is_in(zero_genes)).collect()

            this_m = m_data.lazy().select(pl.first(), col).filter(pl.last().is_not_nan()).filter(
                ~pl.first().is_in(zero_genes)).filter(pl.first().is_in(a_post_cutoff.select(pl.first()))).collect()
            # Trim A and M values
            m_trim_idx = int(np.ceil(log_ratio_trim * len(this_m)))
            a_trim_idx = int(np.ceil(sum_trim * len(a_post_cutoff)))
            trimmed_m = this_m.sort(pl.last())[m_trim_idx:len(this_m) - m_trim_idx]
            trimmed_a = a_post_cutoff.sort(pl.last())[a_trim_idx:len(a_post_cutoff) - a_trim_idx]
            # Get list of genes post trimming
            genes_post_trimming = parsing.data_to_set(trimmed_m.select(pl.first())).intersection(
                parsing.data_to_set(trimmed_a.select(pl.first())))
            # weighted average of trimmed M-values
            if len(genes_post_trimming) == 0:
                tmm = 0
            else:
                tmm = np.average(
                    trimmed_m.lazy().filter(pl.first().is_in(genes_post_trimming)).sort(pl.first()).select(
                        pl.col(col)).collect(),
                    weights=weights.lazy().filter(pl.first().is_in(genes_post_trimming)).sort(pl.first()).select(
                        pl.col(col)).collect())
            scaling_factors[col] = 2 ** tmm

        # adjust scaling factors to multiply, for symmetry, to 1
        scaling_factors = pl.DataFrame(scaling_factors)
        scaling_factors = scaling_factors / np.exp(np.mean(np.log(scaling_factors)))

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @readable_name('Normalize with the "Relative Log Expression" (RLE) method')
    def normalize_rle(self, inplace: bool = True, return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the 'Relative Log Expression' (RLE) method \
        `(Anders and Huber 2010) <https://doi.org/10.1186/gb-2010-11-10-r106>`_. \
        This is the default normalization method used by R's DESeq2. \
        To calculate the Relative Log Expression scaling factors, you first generate a pseudo-sample by calculating \
        the geometric mean expression of each gene across samples. You then calculate the gene-wise ratio \
        of expression between each sample and the pseudo-sample. You then pick the median ratio within each \
        sample as the scaling factor for that sample.

        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_rle()
           Normalized 22 features. Normalized inplace.
        """
        self._validate_is_normalized(expect_normalized=False)

        suffix = '_normRLE'
        data = self.df.select(pl.col(self._numeric_columns)).drop_nulls()
        pseudo_sample = pl.Series(gmean(data, axis=1))
        ratios = (data / pseudo_sample).fill_nan(None).drop_nulls()
        if ratios.shape[0] == 0:
            raise ArithmeticError("Every gene in the count matrix contains at least one zero, "
                                  "cannot calculate geometric means. ")
        scaling_factors = ratios.median()

        # TODO: implement a 'control genes' parameter that calculates ratios only for the given control genes

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @readable_name('Normalize with the "Median of Ratios" (MRN) method')
    def normalize_median_of_ratios(self, sample_grouping: param_typing.GroupedColumns,
                                   reference_group: NonNegativeInt = 0, inplace: bool = True,
                                   return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the 'Median of Ratios Normalization' (MRN) method \
        `(Maza et al 2013) <https://doi.org/10.4161%2Fcib.25849>`_. \
        This normalization method uses information about the experimental condition of each sample. \
        To calculate the Median of Ratios scaling factors, you first calculate the weighted mean expression of \
        each gene within the replicates of each experimental condition. You then calculate per gene the ratio between \
        each weighted mean in the experimental condition and those of the reference condition. \
        You then pick the median ratio for each experimental condition, and calculate the scaling factor for each \
        sample by multiplying it with the sample's total number of reads. \
        Finally, the scaling factors are adjusted, for symmetry, so that they multiply to 1.

        :type sample_grouping: nested list of column names
        :param sample_grouping: grouping of the samples into conditions. \
        Each grouping should containg all replicates of the same condition.
        :type reference_group: int (default=0)
        :param reference_group: the index of the sample group to be used as the reference condition. \
        Must be an integer between 0 and the number of sample groups -1.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_median_of_ratios([['cond1','cond2'],['cond3','cond4']])
           Normalized 22 features. Normalized inplace.
        """
        self._validate_is_normalized(expect_normalized=False)

        flat_grouping = parsing.flatten(sample_grouping)
        assert len(flat_grouping) >= len(self._numeric_columns), f"'sample_grouping' must include all columns. " \
                                                                 f"Only {len(flat_grouping)} out of " \
                                                                 f"{len(self._numeric_columns)} " \
                                                                 f"numeric columns were included. "
        assert isinstance(reference_group, int) and reference_group >= 0, \
            f"Invalid value for 'reference_group': {reference_group}"
        assert reference_group < len(sample_grouping), f"'reference_group' value {reference_group} " \
                                                       f"is larger than the number of sample groups!"

        suffix = '_normMRN'
        data = self.df.select(pl.col(self._numeric_columns))
        weighed_expression = data / data.sum()
        weighted_means = []
        for grp in sample_grouping:
            weighted_means.append(weighed_expression.select(pl.col(grp)).mean_horizontal())
        scaling_factors = {}
        for i, grp in enumerate(sample_grouping):
            ratios = weighted_means[i] / weighted_means[reference_group]
            ratios[ratios == np.inf] = np.nan
            median_of_ratios = ratios.median()
            for cond in grp:
                scaling_factors[cond] = median_of_ratios * self.df[cond].sum()
        scaling_factors = pl.DataFrame(scaling_factors)

        # adjust scaling factors to multiply, for symmetry, to 1
        scaling_factors = scaling_factors / gmean(scaling_factors, axis=1)[0]

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    @readable_name('Normalize with pre-calculated scaling factors')
    def normalize_with_scaling_factors(self, scaling_factor_fname: Union[str, Path], inplace: bool = True):

        """
        Normalizes the reads in the CountFilter using pre-calculated scaling factors. \
        Receives a table of sample names and their corresponding scaling factors, \
        and divides each column in the CountFilter by the corresponding scaling factor.

        :type scaling_factor_fname: str or pathlib.Path
        :param scaling_factor_fname: the .csv file which contains scaling factors for the different libraries.
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_with_scaling_factors("tests/test_files/scaling_factors.csv")
           Normalized 22 features. Normalized inplace.

        """
        self._validate_is_normalized(expect_normalized=False)

        suffix = '_normwithscalingfactors'
        if isinstance(scaling_factor_fname, (str, Path)):
            scaling_factors = io.load_table(scaling_factor_fname)
        elif isinstance(scaling_factor_fname, pl.DataFrame):
            scaling_factors = scaling_factor_fname
        else:
            raise TypeError("Invalid type for 'scaling_factor_fname'!")
        new_df = self._norm_scaling_factors(scaling_factors)
        return self._inplace(new_df, opposite=False, inplace=inplace, suffix=suffix, printout_operation='normalize',
                             _is_normalized=True)

    def _get_ma_ref_column(self, ref_column):
        if isinstance(ref_column, int):
            ref_column += 1
            assert ref_column < len(self.columns), "the index of 'ref_column' is larger than the number of columns!"
            ref_column = self.columns[ref_column]
        elif isinstance(ref_column, str):
            if ref_column.lower() == 'auto':
                upper_quartiles = self.df.select(pl.col(self._numeric_columns)).quantile(0.75)
                ref_index = np.argmin(abs((upper_quartiles - upper_quartiles.mean_horizontal()[0]).to_numpy()))
                ref_column = self.columns[ref_index]
            else:
                assert ref_column in self._numeric_columns, "'ref_column' must be a name of a numeric column!"
        else:
            raise TypeError(f"Invalid type for 'ref_column': {type(ref_column)}")
        return ref_column

    def _calculate_ma(self, ref_column: param_typing.ColumnName, columns: List[str]):
        indices = self.df.lazy().filter(pl.col(ref_column) > 0).select(pl.first()).collect()
        data = self.df.lazy().select(pl.col(self._numeric_columns)).filter(pl.col(ref_column) > 0).collect()
        norm_data = data.select(column / data.sum()[column.name] for column in data.iter_columns()).lazy()

        m_data = norm_data.with_columns(pl.all().truediv(pl.col(ref_column)).log(2).replace([np.inf, -np.inf], np.nan))
        a_data = norm_data.with_columns(
            pl.all().mul(pl.col(ref_column)).log(2).truediv(2).replace([np.inf, -np.inf], np.nan))

        weights = norm_data.with_columns(pl.all().mul(-1).add(1)).collect() / data
        weights = (weights + (weights[ref_column]))
        weights = weights / (weights * weights)
        weights = weights.with_columns(pl.all().replace([np.inf, -np.inf], np.nan))

        m_data = indices.with_columns(m_data.collect())
        a_data = indices.with_columns(a_data.collect())
        weights = indices.with_columns(weights)

        return m_data, a_data, weights

    @readable_name('MA plot')
    def ma_plot(self, ref_column: Union[Literal['auto'], param_typing.ColumnName] = 'auto',
                columns: Union[param_typing.ColumnNames, Literal['all']] = 'all', split_plots: bool = False,
                title: Union[str, Literal['auto']] = 'auto', title_fontsize: float = 20,
                label_fontsize: Union[float, Literal['auto']] = 'auto', tick_fontsize: float = 12) -> List[plt.Figure]:
        """
        Generates M-A (log-ratio vs. log-intensity) plots for selected columns in the dataset. \
        This plot is particularly useful for indicating whether a dataset is properly normalized.

        :param ref_column: the column to be used as reference for MA plot. If 'auto', \
        then the reference column will be chosen automatically to be \
        the column whose upper quartile is closest to the mean upper quartile.
        :type ref_column: name of a column or 'auto' (default='auto')
        :param columns: A list of the column names to generate an MA plot for.
        :type columns: str or list of str
        :param split_plots: if True, each individual MA plot will be plotted in its own Figure. \
        Otherwise, all MA plots will be plotted on the same Figure.
        :type split_plots: bool (default=False)
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
         :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :rtype: a list of matplotlib Figures.
        """
        ref_column = self._get_ma_ref_column(ref_column)
        if isinstance(columns, str) and columns.lower() == 'all':
            columns = self._numeric_columns
        else:
            columns = parsing.data_to_list(columns)
        if ref_column in columns:
            columns.remove(ref_column)

        m_data, a_data, _ = self._calculate_ma(ref_column, columns)

        if split_plots:
            subplots = [plt.subplots(constrained_layout=True, figsize=(6.4, 5.6)) for _ in range(len(columns))]
            if label_fontsize == 'auto':
                label_fontsize = 14
        else:
            grid = strategies.SquareStrategy()
            subplots = grid.get_grid(len(columns))
            plt.close()
            fig = plt.figure(figsize=(14, 9))
            if title == 'auto':
                title = 'MA plots'
            fig.suptitle(title, fontsize=title_fontsize)
            if label_fontsize == 'auto':
                label_fontsize = 10 if len(columns) > 15 else 8

        for i, col in enumerate(columns):
            if split_plots:
                ax = subplots[i][1]
                fig = subplots[i][0]
                title = f'MA plot of {col} vs {ref_column}'
            else:
                # noinspection PyUnboundLocalVariable
                ax = fig.add_subplot(subplots[i])
                title = f'{col} vs {ref_column}'
            ax.scatter(a_data[:, i], m_data[:, i], c='black', s=3, alpha=0.3)
            ax.axhline(0, c='green')
            ax.set_title(title)
            ax.set_xlabel('A: \n' + r'$\log_2$ average expression)', fontsize=label_fontsize)
            ax.set_ylabel('M: \n' + r'$\log_2$' + f'("{col}" / "{ref_column}")', fontsize=label_fontsize)
            ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)

        if not split_plots:
            fig.tight_layout(rect=[0, 0.03, 1, 0.92])
        plt.show()

        if split_plots:
            return [subplots[i][0] for i in range(len(columns))]
        return [fig]

    @readable_name('Filter genes with low expression in all columns')
    def filter_low_reads(self, threshold: float = 5, n_samples: PositiveInt = 1, opposite: bool = False,
                         inplace: bool = True):

        """
        Filter out features which are lowly-expressed in all columns, keeping only features with at least 'threshold' \
        reads in at least 'n_samples' columns.

        :type threshold: float
        :param threshold: The minimal number of reads (counts, rpm, rpkm, tpm, etc) a feature should have \
        in at least `n_samples` samples in order not to be filtered out.
        :param n_samples: the minimal number of samples a feature should have at least 'threshold' reads in \
        in order not to be filtered out.
        :type n_samples: positive integer (default=1)
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of CountFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> c.filter_low_reads(5) # remove all rows whose values in all columns are all <5
            Filtered 6 features, leaving 16 of the original 22 features. Filtered inplace.

        """
        validation.validate_threshold(threshold)
        self._validate_is_normalized()
        mask = pl.when(pl.col(self._numeric_columns) >= threshold).then(True).otherwise(False).name.map(lambda x: x)
        mask = self.df.select(pl.col(self._numeric_columns)).with_columns(mask).sum_horizontal() >= n_samples
        new_df = self.df.filter(mask)
        suffix = f"_filt{threshold}reads{n_samples}samples"
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Split into Higly and Lowly expressed genes')
    def split_by_reads(self, threshold: float = 5) -> tuple:

        """
        Splits the features in the CountFilter object into two non-overlapping CountFilter \
        objects, based on their maximum expression level. The first object will contain only highly-expressed \
         features (which have reads over the specified threshold in at least one sample). The second object will \
         contain only lowly-expressed features (which have reads below the specified threshold in all samples).

        :param threshold: The minimal number of reads (counts, RPM, RPKM, TPM etc) a feature needs to have \
        in at least one sample in order to be \
        included in the "highly expressed" object and no the "lowly expressed" object.
        :type threshold: float (default=5)
        :rtype: Tuple[filtering.CountFilter, filtering.CountFilter]
        :return: A tuple containing two CountFilter objects: the first has only highly-expressed features, \
        and the second has only lowly-expressed features.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> low_expression, high_expression = c.split_by_reads(5)
            Filtered 6 features, leaving 16 of the original 22 features. Filtering result saved to new object.
            Filtered 16 features, leaving 6 of the original 22 features. Filtering result saved to new object.

        """
        validation.validate_threshold(threshold)
        self._validate_is_normalized()
        high_expr = self.df.filter(self.df.select(pl.col(self._numeric_columns)).max_horizontal() >= threshold)
        low_expr = self.df.filter(self.df.select(pl.col(self._numeric_columns)).max_horizontal() < threshold)
        return self._inplace(high_expr, opposite=False, inplace=False, suffix=f'_above{threshold}reads'), self._inplace(
            low_expr, opposite=False, inplace=False, suffix=f'_below{threshold}reads')

    @readable_name('Filter genes with low summized expression across all conditions')
    def filter_by_row_sum(self, threshold: float = 5, opposite: bool = False, inplace: bool = True):

        """
        Removes features/rows whose sum is belove 'threshold'.

        :type threshold: float
        :param threshold: The minimal sum a row should have in order not to be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of CountFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> c.filter_by_row_sum(5) # remove all rows whose sum is <5
            Filtered 4 features, leaving 18 of the original 22 features. Filtered inplace.

        """
        validation.validate_threshold(threshold)
        self._validate_is_normalized()

        new_df = self.df.filter(self.df.select(pl.col(self._numeric_columns)).sum_horizontal() >= threshold)
        suffix = f"_filt{threshold}sum"
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('K-Means clustering')
    def split_kmeans(self, n_clusters: Union[PositiveInt, List[PositiveInt], Literal[K_CRITERIA]],
                     n_init: PositiveInt = 3, max_iter: PositiveInt = 300,
                     random_seed: Union[NonNegativeInt, None] = None, power_transform: bool = True,
                     plot_style: Literal['all', 'std_area', 'std_bar'] = 'all',
                     split_plots: bool = False, max_n_clusters_estimate: Union[PositiveInt, Literal['auto']] = 'auto',
                     parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                     gui_mode: bool = False) -> Union[Tuple['CountFilter', ...], Tuple[Tuple['CountFilter', ...], ...]]:
        """
        Clusters the features in the CountFilter object using the K-means clustering algorithm, \
        and then splits those features into multiple non-overlapping CountFilter objects, \
        based on the clustering result.

        :param n_clusters: The number of clusters the algorithm will seek.
        :type n_clusters: int, list of ints, 'gap', 'silhouette', 'calinski_harabasz', 'davies_bouldin', or 'bic'
        :param random_seed: determines random number generation for centroid initialization. \
        Use an int to make the randomness deterministic.
        :type random_seed: Union[int, None] or None (default=None)
        :param n_init: number of time the k-medoids algorithm will be run with different medoid seeds. \
        The final results will be the best output of n_init consecutive runs in terms of inertia.
        :type n_init: int (default=3)
        :param max_iter: maximum number of iterations of the k-medoids algorithm for a single run.
        :type max_iter: int (default=300)
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to clustering.
        :type power_transform: bool (default=True)
        :param plot_style: determines the visual style of the cluster expression plot.
        :type plot_style: 'all', 'std_area', or 'std_bar' (default='all')
        :param split_plots: if True, each discovered cluster will be plotted on its own. \
        Otherwise, all clusters will be plotted in the same Figure.
        :type split_plots: bool (default=False)
        :param max_n_clusters_estimate: the maximum number of clusters to test if trying to automatically estimate \
        the optimal number of clusters. If `max_n_clusters_estimate`='default', \
        an appropriate value will be picked automatically.
        :type max_n_clusters_estimate: int or 'auto' (default='auto')
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :return: if `n_clusters` is an int, returns a tuple of `n_clusters` CountFilter objects, \
        each corresponding to a discovered cluster. \
        If `n_clusters` is a list, returns one tuple of CountFilter objects per value in `n_clusters`.

        :Examples:
            >>> from rnalysis import filtering
            >>> dev_stages = filtering.CountFilter('tests/test_files/elegans_developmental_stages.tsv')
            >>> dev_stages.filter_low_reads(100)
            Filtered 44072 features, leaving 2326 of the original 46398 features. Filtered inplace.
            >>> clusters = dev_stages.split_kmeans(14,power_transform=True)
            Filtered 44072 features, leaving 2326 of the original 46398 features. Filtered inplace.
            Filtered 1801 features, leaving 525 of the original 2326 features. Filtering result saved to new object.
            Filtered 2010 features, leaving 316 of the original 2326 features. Filtering result saved to new object.
            Filtered 2059 features, leaving 267 of the original 2326 features. Filtering result saved to new object.
            Filtered 2102 features, leaving 224 of the original 2326 features. Filtering result saved to new object.
            Filtered 2185 features, leaving 141 of the original 2326 features. Filtering result saved to new object.
            Filtered 2186 features, leaving 140 of the original 2326 features. Filtering result saved to new object.
            Filtered 2200 features, leaving 126 of the original 2326 features. Filtering result saved to new object.
            Filtered 2219 features, leaving 107 of the original 2326 features. Filtering result saved to new object.
            Filtered 2225 features, leaving 101 of the original 2326 features. Filtering result saved to new object.
            Filtered 2225 features, leaving 101 of the original 2326 features. Filtering result saved to new object.
            Filtered 2241 features, leaving 85 of the original 2326 features. Filtering result saved to new object.
            Filtered 2250 features, leaving 76 of the original 2326 features. Filtering result saved to new object.
            Filtered 2259 features, leaving 67 of the original 2326 features. Filtering result saved to new object.
            Filtered 2276 features, leaving 50 of the original 2326 features. Filtering result saved to new object.


        .. figure:: /figures/kmeans_all.png
           :align:   center

           Example plot of split_kmeans(plot_style='all')

        .. figure:: /figures/clustering_PCA_kmeans.png
           :align:   center

           Example plot of split_kmeans()
        """
        runner = clustering.KMeansRunner(self.df.select(pl.col(self._numeric_columns)), power_transform, n_clusters,
                                         max_n_clusters_estimate, random_seed, n_init, max_iter, plot_style,
                                         split_plots, parallel_backend)
        clusterers = runner.run(plot=not gui_mode)
        filt_obj_tuples = []
        for clusterer in clusterers:
            # split the CountFilter object
            filt_obj_tuples.append(
                tuple([self._inplace(self.df.filter(clusterer.labels_ == i), opposite=False, inplace=False,
                                     suffix=f'_kmeanscluster{i + 1}') for i in range(clusterer.n_clusters_)]))
        # if only a single K was calculated, don't return it as a tuple of length 1
        return_val = filt_obj_tuples[0] if len(filt_obj_tuples) == 1 else parsing.data_to_tuple(filt_obj_tuples)
        return (return_val, runner) if gui_mode else return_val

    @readable_name('Hierarchical clustering')
    def split_hierarchical(self, n_clusters: Union[
        PositiveInt, List[PositiveInt], Literal[K_CRITERIA + ('distance',)]],
                           metric: Literal['Euclidean', 'Cosine', 'Pearson', 'Spearman', 'Manhattan',
                           'L1', 'L2', 'Jackknife', 'YS1', 'YR1', 'Sharpened_Cosine'] = 'Euclidean',
                           linkage: Literal['Single', 'Average', 'Complete', 'Ward'] = 'Average',
                           power_transform: bool = True, distance_threshold: Union[float, None] = None,
                           plot_style: Literal['all', 'std_area', 'std_bar'] = 'all', split_plots: bool = False,
                           max_n_clusters_estimate: Union[PositiveInt, Literal['auto']] = 'auto',
                           parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                           gui_mode: bool = False
                           ) -> Union[Tuple['CountFilter', ...], Tuple[Tuple['CountFilter', ...], ...]]:
        """
        Clusters the features in the CountFilter object using the Hierarchical clustering algorithm, \
        and then splits those features into multiple non-overlapping CountFilter objects, \
        based on the clustering result.

        :param n_clusters: The number of clusters the algorithm will seek. If set to 'distance', \
        the algorithm will derive the number of clusters from the distance threshold (see 'distance_threshold').
        :type n_clusters: int, list of ints, 'distance', 'gap', 'silhouette', 'calinski_harabasz', 'davies_bouldin', \
        or 'bic'
        :param metric: the distance metric used to determine similarity between data points. \
        If linkage is 'ward', only the 'Euclidean' metric is accepted. \
        For a full list of supported distance metrics see the user guide.
        :type metric: 'Euclidean', 'l1', 'l2', 'manhattan', or 'cosine',  (default='Euclidean')
        :param linkage: Which linkage criterion to use. The linkage criterion determines which distance to use between \
        sets of observation. The algorithm will merge the pairs of cluster that minimize this criterion.
        :type linkage: 'single', 'Average', 'complete', or 'ward' (default='Average')
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to clustering.
        :type power_transform: bool (default=True)
        :param distance_threshold: a distance threshold above which clusters will not be merged. \
        If a number is specified, `n_clusters` must be None.
        :type distance_threshold: float or None (default=None)
        :param plot_style: determines the visual style of the cluster expression plot.
        :type plot_style: 'all', 'std_area', or 'std_bar' (default='all')
        :param split_plots: if True, each discovered cluster will be plotted on its own. \
        Otherwise, all clusters will be plotted in the same Figure.
        :type split_plots: bool (default=False)
        :param max_n_clusters_estimate: the maximum number of clusters to test if trying to automatically estimate \
        the optimal number of clusters. If `max_n_clusters_estimate`='default', \
        an appropriate value will be picked automatically.
        :type max_n_clusters_estimate: int or 'auto' (default='auto')
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :return: if `n_clusters` is an int, returns a tuple of `n_clusters` CountFilter objects, \
        each corresponding to a discovered cluster. \
        If `n_clusters` is a list, returns one tuple of CountFilter objects per value in `n_clusters`.


        :Examples:
            >>> from rnalysis import filtering
            >>> dev_stages = filtering.CountFilter('tests/test_files/elegans_developmental_stages.tsv')
            >>> dev_stages.filter_low_reads(100)
            Filtered 44072 features, leaving 2326 of the original 46398 features. Filtered inplace.
            >>> clusters = dev_stages.split_hierarchical(n_clusters=13, metric='Euclidean',linkage='ward'
            ...                                         ,power_transform=True)
            Filtered 1718 features, leaving 608 of the original 2326 features. Filtering result saved to new object.
            Filtered 1979 features, leaving 347 of the original 2326 features. Filtering result saved to new object.
            Filtered 2094 features, leaving 232 of the original 2326 features. Filtering result saved to new object.
            Filtered 2110 features, leaving 216 of the original 2326 features. Filtering result saved to new object.
            Filtered 2156 features, leaving 170 of the original 2326 features. Filtering result saved to new object.
            Filtered 2191 features, leaving 135 of the original 2326 features. Filtering result saved to new object.
            Filtered 2195 features, leaving 131 of the original 2326 features. Filtering result saved to new object.
            Filtered 2223 features, leaving 103 of the original 2326 features. Filtering result saved to new object.
            Filtered 2224 features, leaving 102 of the original 2326 features. Filtering result saved to new object.
            Filtered 2238 features, leaving 88 of the original 2326 features. Filtering result saved to new object.
            Filtered 2246 features, leaving 80 of the original 2326 features. Filtering result saved to new object.
            Filtered 2252 features, leaving 74 of the original 2326 features. Filtering result saved to new object.
            Filtered 2286 features, leaving 40 of the original 2326 features. Filtering result saved to new object.

        .. figure:: /figures/hierarchical_all.png
           :align:   center

           Example plot of split_hierarchical(plot_style='all')

                   .. figure:: /figures/clustering_PCA_hierarchical.png
           :align:   center

           Example plot of split_hierarchical()
        """
        runner = clustering.HierarchicalRunner(self.df.select(pl.col(self._numeric_columns)), power_transform,
                                               n_clusters,
                                               max_n_clusters_estimate, metric, linkage, distance_threshold, plot_style,
                                               split_plots, parallel_backend)
        clusterers = runner.run(plot=not gui_mode)
        filt_obj_tuples = []
        for clusterer in clusterers:
            # split the CountFilter object
            this_n_clusters = np.max(np.unique(clusterer.labels_)) + 1
            filt_obj_tuples.append(
                tuple([self._inplace(self.df.filter(clusterer.labels_ == i), opposite=False, inplace=False,
                                     suffix=f'_kmedoidscluster{i + 1}') for i in range(this_n_clusters)]))
        # if only a single K was calculated, don't return it as a tuple of length 1
        return_val = filt_obj_tuples[0] if len(filt_obj_tuples) == 1 else parsing.data_to_tuple(filt_obj_tuples)
        return (return_val, runner) if gui_mode else return_val

    @readable_name('K-Medoids clustering')
    def split_kmedoids(self, n_clusters: Union[PositiveInt, List[PositiveInt], Literal[K_CRITERIA]],
                       n_init: PositiveInt = 3, max_iter: PositiveInt = 300,
                       random_seed: Union[NonNegativeInt, None] = None,
                       metric: Union[str, Literal['Euclidean', 'Cosine', 'Pearson', 'Spearman', 'Manhattan',
                       'L1', 'L2', 'Jackknife', 'YS1', 'YR1', 'Sharpened_Cosine', 'Hamming']] = 'Euclidean',
                       power_transform: bool = True,
                       plot_style: Literal['all', 'std_area', 'std_bar'] = 'all', split_plots: bool = False,
                       max_n_clusters_estimate: Union[PositiveInt, Literal['auto']] = 'auto',
                       parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky', gui_mode: bool = False
                       ) -> Union[Tuple['CountFilter', ...], Tuple[Tuple['CountFilter', ...], ...]]:
        """
        Clusters the features in the CountFilter object using the K-medoids clustering algorithm, \
        and then splits those features into multiple non-overlapping CountFilter objects, \
        based on the clustering result.

        :param n_clusters: The number of clusters the algorithm will seek.
        :type n_clusters: int, list of ints, 'gap', 'silhouette', 'calinski_harabasz', 'davies_bouldin', or 'bic'
        :param random_seed: determines random number generation for centroid initialization. \
        Use an int to make the randomness deterministic.
        :type random_seed: Union[int, None] or None (default=None)
        :param n_init: number of time the k-medoids algorithm will be run with different medoid seeds. \
        The final results will be the best output of n_init consecutive runs in terms of inertia.
        :type n_init: int (default=3)
        :param max_iter: maximum number of iterations of the k-medoids algorithm for a single run.
        :type max_iter: int (default=300)
        :param metric: the distance metric used to determine similarity between data points. \
        For a full list of supported distance metrics see the user guide.
        :type metric: str (default='Euclidean')
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to clustering.
        :type power_transform: bool (default=True)
        :param plot_style: determines the visual style of the cluster expression plot.
        :type plot_style: 'all', 'std_area', or 'std_bar' (default='all')
        :param split_plots: if True, each discovered cluster will be plotted on its own. \
        Otherwise, all clusters will be plotted in the same Figure.
        :type split_plots: bool (default=False)
        :param max_n_clusters_estimate: the maximum number of clusters to test if trying to automatically estimate \
        the optimal number of clusters. If `max_n_clusters_estimate`='default', \
        an appropriate value will be picked automatically.
        :type max_n_clusters_estimate: int or 'auto' (default='auto')
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :return: if `n_clusters` is an int, returns a tuple of `n_clusters` CountFilter objects, \
        each corresponding to a discovered cluster. \
        If `n_clusters` is a list, returns one tuple of CountFilter objects per value in `n_clusters`.

        :Examples:
            >>> from rnalysis import filtering
            >>> dev_stages = filtering.CountFilter('tests/test_files/elegans_developmental_stages.tsv')
            >>> dev_stages.filter_low_reads(100)
            Filtered 44072 features, leaving 2326 of the original 46398 features. Filtered inplace.
            >>> clusters = dev_stages.split_kmedoids(n_clusters=14, metric='spearman', power_transform=True)
            Filtered 1967 features, leaving 359 of the original 2326 features. Filtering result saved to new object.
            Filtered 2020 features, leaving 306 of the original 2326 features. Filtering result saved to new object.
            Filtered 2071 features, leaving 255 of the original 2326 features. Filtering result saved to new object.
            Filtered 2131 features, leaving 195 of the original 2326 features. Filtering result saved to new object.
            Filtered 2145 features, leaving 181 of the original 2326 features. Filtering result saved to new object.
            Filtered 2157 features, leaving 169 of the original 2326 features. Filtering result saved to new object.
            Filtered 2159 features, leaving 167 of the original 2326 features. Filtering result saved to new object.
            Filtered 2182 features, leaving 144 of the original 2326 features. Filtering result saved to new object.
            Filtered 2190 features, leaving 136 of the original 2326 features. Filtering result saved to new object.
            Filtered 2192 features, leaving 134 of the original 2326 features. Filtering result saved to new object.
            Filtered 2229 features, leaving 97 of the original 2326 features. Filtering result saved to new object.
            Filtered 2252 features, leaving 74 of the original 2326 features. Filtering result saved to new object.
            Filtered 2268 features, leaving 58 of the original 2326 features. Filtering result saved to new object.
            Filtered 2275 features, leaving 51 of the original 2326 features. Filtering result saved to new object.


        .. figure:: /figures/kmedoids_all.png
           :align:   center

           Example plot of split_kmedoids(plot_style='all')

        .. figure:: /figures/clustering_PCA_kmedoids.png
           :align:   center

           Example plot of split_kmedoids()
        """
        runner = clustering.KMedoidsRunner(self.df.select(pl.col(self._numeric_columns)), power_transform, n_clusters,
                                           max_n_clusters_estimate, metric, random_seed, n_init, max_iter, plot_style,
                                           split_plots, parallel_backend)
        clusterers = runner.run(plot=not gui_mode)
        filt_obj_tuples = []
        for clusterer in clusterers:
            # split the CountFilter object
            filt_obj_tuples.append(
                tuple([self._inplace(self.df.filter(clusterer.labels_ == i), opposite=False, inplace=False,
                                     suffix=f'_kmedoidscluster{i + 1}') for i in range(clusterer.n_clusters_)]))
        # if only a single K was calculated, don't return it as a tuple of length 1
        return_val = filt_obj_tuples[0] if len(filt_obj_tuples) == 1 else parsing.data_to_tuple(filt_obj_tuples)
        return (return_val, runner) if gui_mode else return_val

    @readable_name('CLICOM (ensemble) clustering')
    def split_clicom(self, *parameter_dicts: dict,
                     replicate_grouping: Union[param_typing.GroupedColumns, Literal['ungrouped']] = 'ungrouped',
                     power_transform: Union[bool, Tuple[bool, bool]] = True,
                     evidence_threshold: param_typing.Fraction = 2 / 3, cluster_unclustered_features: bool = False,
                     min_cluster_size: PositiveInt = 15, plot_style: Literal['all', 'std_area', 'std_bar'] = 'all',
                     split_plots: bool = False, parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                     gui_mode: bool = False) -> Tuple['CountFilter', ...]:
        """
        Clusters the features in the CountFilter object using the modified CLICOM ensemble clustering algorithm \
        `(Mimaroglu and Yagci 2012) <https://doi.org/10.1016/j.eswa.2011.08.059/>`_, \
        and then splits those features into multiple non-overlapping CountFilter objects, \
        based on the clustering result. \
        The CLICOM algorithm incorporates the results of multiple clustering solutions, \
        which can come from different clustering algorithms with differing clustering parameters, \
        and uses these clustering solutions to create a combined clustering solution. \
        Due to the nature of CLICOM, the number of clusters the data will be divided into is determined automatically. \
        This modified version of the CLICOM algorithm can also classify features as noise, \
        which does not belong in any discovered cluster.

        :param replicate_grouping: Allows you to group your data into replicates. Each replicate will be clustered \
        separately, and used as its own clustering setup. This can minimize the influence of batch effects \
        on the clustering results, and take advantage of repeated measures data \
        to improve the accuracy of your clustering. If `replicate_grouping`='ungrouped', the data will be \
        clustered normally as if no replicate data is available. \
        To read more about the theory behind this, see the following publication: \
        https://doi.org/10.1093/bib/bbs057
        :type replicate_grouping: nested list of strings or 'ungrouped' (default='ungrouped')
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to clustering. If both True and False are supplied, \
        RNAlysis will run the initial clustering setups twice: once with a power transform, and once without.
        :type power_transform: True, False, or (True, False) (default=True)
        :param evidence_threshold: determines whether each pair of features can be reliably clustered together. \
        For example, if evidence_threshold=0.5, a pair of features is considered reliably clustered together if \
        they were clustered together in at least 50% of the tested clustering solutions.
        :type evidence_threshold: float between 0 and 1 (default=2/3)
        :param cluster_unclustered_features: if True, RNAlysis will force every feature to be part of a cluster, \
        even if they were not initially determined to reliably belong to any of the discovered clusters. \
        Larger values will lead to fewer clusters, with more features classified as noise.
        :type cluster_unclustered_features: bool (default=False)
        :param min_cluster_size: the minimum size of clusters the algorithm will seek. \
        Larger values will lead to fewer clusters, with more features classified as noise.
        :type min_cluster_size: int (default=15)
        :param parameter_dicts: multiple dictionaries, each corresponding to a clustering setup to be run. \
        Each dictionary must contain a 'method' field with a clustering method supported by RNAlysis \
        ('k-means', 'k-medoids', 'hierarchical', or 'hdbscan'). \
        The other fields of the dictionary should contain your preferred values \
        for each of the clustering algorithm's parameters. \
        Yoy can specify a list of values for each of those parameters, \
        and then RNAlysis will run the clustering algorithm with all legal combinations of parameters you specified. \
        For example, {'method':'k-medoids', 'n_clusters':[3,5], 'metric':['Euclidean', 'cosine']} \
        will run the K-Medoids algorithm four times with the following parameter combinations: \
        (n_clusters=3,metric='Euclidean'), (n_clusters=5, metric='Euclidean'), \
        (n_clusters=3, metric='cosine'), (n_clusters=5, metric='cosine').
        :param plot_style: determines the visual style of the cluster expression plot.
        :type plot_style: 'all', 'std_area', or 'std_bar' (default='all')
        :param split_plots: if True, each discovered cluster will be plotted on its own. \
        Otherwise, all clusters will be plotted in the same Figure.
        :type split_plots: bool (default=False)
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :return: returns a tuple of CountFilter objects, each corresponding to a discovered cluster.


        :Examples:
            >>> from rnalysis import filtering
            >>> dev_stages = filtering.CountFilter('tests/test_files/elegans_developmental_stages.tsv')
            >>> dev_stages.filter_low_reads(100)
            Filtered 44072 features, leaving 2326 of the original 46398 features. Filtered inplace.
            >>> clusters = dev_stages.split_clicom(
            ... {'method': 'hdbscan', 'min_cluster_size': [50, 75, 140], 'metric': ['ys1', 'yr1', 'spearman']},
            ... {'method': 'hierarchical', 'n_clusters': [7, 12], 'metric': ['Euclidean', 'jackknife', 'yr1'],
            ... 'linkage': ['average', 'ward']}, {'method': 'kmedoids', 'n_clusters': [7, 16], 'metric': 'spearman'},
            ... power_transform=True, evidence_threshold=0.5, min_cluster_size=40)
            Found 19 legal clustering setups.
            Running clustering setups: 100%|██████████| 19/19 [00:12<00:00,  1.49 setup/s]
            Generating cluster similarity matrix: 100%|██████████| [00:32<00:00, 651.06it/s]
            Finding cliques: 100%|██████████| 42436/42436 [00:00<00:00, 61385.87it/s]
            Done
            Found 15 clusters of average size 153.60. Number of unclustered genes is 22, which are 0.95% of the genes.
            Filtered 1864 features, leaving 462 of the original 2326 features. Filtering result saved to new object.
            Filtered 2115 features, leaving 211 of the original 2326 features. Filtering result saved to new object.
            Filtered 2122 features, leaving 204 of the original 2326 features. Filtering result saved to new object.
            Filtered 2123 features, leaving 203 of the original 2326 features. Filtering result saved to new object.
            Filtered 2128 features, leaving 198 of the original 2326 features. Filtering result saved to new object.
            Filtered 2167 features, leaving 159 of the original 2326 features. Filtering result saved to new object.
            Filtered 2179 features, leaving 147 of the original 2326 features. Filtering result saved to new object.
            Filtered 2200 features, leaving 126 of the original 2326 features. Filtering result saved to new object.
            Filtered 2204 features, leaving 122 of the original 2326 features. Filtering result saved to new object.
            Filtered 2229 features, leaving 97 of the original 2326 features. Filtering result saved to new object.
            Filtered 2234 features, leaving 92 of the original 2326 features. Filtering result saved to new object.
            Filtered 2238 features, leaving 88 of the original 2326 features. Filtering result saved to new object.
            Filtered 2241 features, leaving 85 of the original 2326 features. Filtering result saved to new object.
            Filtered 2263 features, leaving 63 of the original 2326 features. Filtering result saved to new object.
            Filtered 2279 features, leaving 47 of the original 2326 features. Filtering result saved to new object.


        .. figure:: /figures/clicom_all.png
           :align:   center

           Example plot of split_clicom(plot_style='all')

        .. figure:: /figures/clustering_PCA_clicom.png
           :align:   center

           Example plot of split_clicom()
        """
        if replicate_grouping == 'ungrouped':
            replicate_grouping = [self.columns]
        else:
            assert validation.isinstanceiter(replicate_grouping,
                                             list), "'replicate_grouping' must contain only lists of strings!"
            for grp in replicate_grouping:
                for cond in grp:
                    assert cond in self.columns, f"column '{cond}' does not exist!"

        runner = clustering.CLICOMRunner(self.df.select(pl.col(self._numeric_columns)), replicate_grouping,
                                         power_transform,
                                         evidence_threshold, cluster_unclustered_features, min_cluster_size,
                                         *parameter_dicts, plot_style=plot_style, split_plots=split_plots,
                                         parallel_backend=parallel_backend)
        [clusterer] = runner.run(plot=not gui_mode)
        n_clusters = clusterer.n_clusters_
        if n_clusters == 0:
            print("Found 0 clusters with the given parameters. Please try again with different parameters. ")
        else:
            unclustered = np.count_nonzero(clusterer.labels_ == -1)
            print(f"Found {n_clusters} clusters of average size "
                  f"{(len(clusterer.labels_) - unclustered) / n_clusters  :.2f}. "
                  f"Number of unclustered genes is {unclustered}, "
                  f"which are {100 * (unclustered / len(clusterer.labels_)) :.2f}% of the genes.")

        filt_objs = [self._inplace(self.df.filter(clusterer.labels_ == i), opposite=False, inplace=False,
                                   suffix=f'_clicomcluster{i + 1}') for i in range(n_clusters)]

        return_val = parsing.data_to_tuple(filt_objs)
        return (return_val, runner) if gui_mode else return_val

    @readable_name('HDBSCAN (density) clustering')
    def split_hdbscan(self, min_cluster_size: PositiveInt, min_samples: Union[PositiveInt, None] = 1,
                      metric: Union[str, Literal['Euclidean', 'Cosine', 'Pearson', 'Spearman', 'Manhattan',
                      'L1', 'L2', 'Jackknife', 'YS1', 'YR1', 'Sharpened_Cosine', 'Hamming']] = 'Euclidean',
                      cluster_selection_epsilon: float = 0, cluster_selection_method: Literal['eom', 'leaf'] = 'eom',
                      power_transform: bool = True, plot_style: Literal['all', 'std_area', 'std_bar'] = 'all',
                      split_plots: bool = False, return_probabilities: bool = False,
                      parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky', gui_mode: bool = False
                      ) -> Union[Tuple['CountFilter', ...], List[Union[Tuple['CountFilter', ...], np.ndarray]], None]:
        """
        Clusters the features in the CountFilter object using the HDBSCAN clustering algorithm, \
        and then splits those features into multiple non-overlapping CountFilter objects, \
        based on the clustering result.

        :param min_cluster_size: the minimum size of clusters the algorithm will seek. \
        Larger values will lead to fewer, larger clusters.
        :type min_cluster_size: int
        :param min_samples: the number of samples in a neighbourhood for a point to be considered a core point. \
        Higher values will lead to a more conservative clustering result, \
        with more points being classified as noise. \
        If `min_samples` is None, the algorithm will pick a value automatically
        :type min_samples: int or None (default=1)
        :param metric: the distance metric used to determine similarity between data points. \
        For a full list of supported distance metrics see the user guide.
        :type metric: str (default='Euclidean')
        :param cluster_selection_epsilon: a distance threshold below which clusters will be merged.
        :type cluster_selection_epsilon: float (default=0.0)
        :param cluster_selection_method: The method used to select clusters from the condensed tree. \
        'eom' will use an Excess of Mass algorithm to find the most persistent clusters. \
        'leaf' will select the leaves of the tree, providing the most fine-grained and homogenous clusters.
        :type cluster_selection_method: 'eom' or 'leaf' (default='eom')
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to clustering.
        :type power_transform: bool (default=True)
        :param plot_style: determines the visual style of the cluster expression plot.
        :type plot_style: 'all', 'std_area', or 'std_bar' (default='all')
        :param split_plots: if True, each discovered cluster will be plotted on its own. \
        Otherwise, all clusters will be plotted in the same Figure.
        :type split_plots: bool (default=False)
        :param return_probabilities: if True, the algorithm will return an array containing \
        the probability with which each sample is a member of its assigned cluster, \
        in addition to returning the clustering results. Points which were categorized as noise have probability 0.
        :type return_probabilities: bool (default False)
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :return: if `return_probabilities` is False, returns a tuple of CountFilter objects, \
        each corresponding to a discovered cluster. \
        Otherswise, returns a tuple of CountFilter objects, and a numpy array containing the probability values.

        :Examples:
            >>> from rnalysis import filtering
            >>> dev_stages = filtering.CountFilter('tests/test_files/elegans_developmental_stages.tsv')
            >>> dev_stages.filter_low_reads(100)
            Filtered 44072 features, leaving 2326 of the original 46398 features. Filtered inplace.
            >>> clusters = dev_stages.split_hdbscan(min_cluster_size=75,metric='yr1',power_transform=True)
            Found 14 clusters of average size 141.57. Number of unclustered genes is 344, which are 14.79% of the genes.
            Filtered 2019 features, leaving 307 of the original 2326 features. Filtering result saved to new object.
            Filtered 2122 features, leaving 204 of the original 2326 features. Filtering result saved to new object.
            Filtered 2146 features, leaving 180 of the original 2326 features. Filtering result saved to new object.
            Filtered 2168 features, leaving 158 of the original 2326 features. Filtering result saved to new object.
            Filtered 2173 features, leaving 153 of the original 2326 features. Filtering result saved to new object.
            Filtered 2176 features, leaving 150 of the original 2326 features. Filtering result saved to new object.
            Filtered 2183 features, leaving 143 of the original 2326 features. Filtering result saved to new object.
            Filtered 2192 features, leaving 134 of the original 2326 features. Filtering result saved to new object.
            Filtered 2200 features, leaving 126 of the original 2326 features. Filtering result saved to new object.
            Filtered 2234 features, leaving 92 of the original 2326 features. Filtering result saved to new object.
            Filtered 2238 features, leaving 88 of the original 2326 features. Filtering result saved to new object.
            Filtered 2241 features, leaving 85 of the original 2326 features. Filtering result saved to new object.
            Filtered 2244 features, leaving 82 of the original 2326 features. Filtering result saved to new object.
            Filtered 2246 features, leaving 80 of the original 2326 features. Filtering result saved to new object.

        .. figure:: /figures/hdbscan_all.png
           :align:   center

           Example plot of split_hdbscan(plot_style='all')

        .. figure:: /figures/clustering_PCA_hdbscan.png
           :align:   center

           Example plot of split_hdbscan()
           """
        validation.validate_hdbscan_parameters(min_cluster_size, metric, cluster_selection_method, self.shape[0])

        runner = clustering.HDBSCANRunner(self.df.select(pl.col(self._numeric_columns)), power_transform,
                                          min_cluster_size,
                                          min_samples, metric, cluster_selection_epsilon, cluster_selection_method,
                                          return_probabilities, plot_style, split_plots, parallel_backend)
        if return_probabilities:
            [clusterer], probabilities = runner.run(plot=not gui_mode)
        else:
            [clusterer] = runner.run(plot=not gui_mode)

        if clusterer is None: return  # if hdbscan is not installed, the runner will return None

        n_clusters = clusterer.labels_.max() + 1
        if n_clusters == 0:
            print("Found 0 clusters with the given parameters. Please try again with different parameters. ")
        else:
            unclustered = np.count_nonzero(clusterer.labels_ == -1)
            print(f"Found {n_clusters} clusters of average size "
                  f"{(len(clusterer.labels_) - unclustered) / n_clusters  :.2f}. "
                  f"Number of unclustered genes is {unclustered}, "
                  f"which are {100 * (unclustered / len(clusterer.labels_)) :.2f}% of the genes.")

        filt_objs = parsing.data_to_tuple([self._inplace(self.df.filter(clusterer.labels_ == i), opposite=False,
                                                         inplace=False, suffix=f'_hdbscancluster{i + 1}') for i in
                                           range(n_clusters)])

        # noinspection PyUnboundLocalVariable
        return_val = [filt_objs, probabilities] if return_probabilities else filt_objs
        return (return_val, runner) if gui_mode else return_val

    @readable_name('Hierarchical clustergram plot')
    def clustergram(self, sample_names: Union[param_typing.ColumnNames, Literal['all']] = 'all',
                    metric: Union[Literal['Correlation', 'Cosine', 'Euclidean', 'Jaccard'], str] = 'Euclidean',
                    linkage: Literal['Single', 'Average', 'Complete', 'Ward', 'Weighted', 'Centroid', 'Median'
                    ] = 'Average', title: Union[str, Literal['auto']] = 'auto', title_fontsize: float = 20,
                    tick_fontsize: float = 12, colormap: ColorMap = 'inferno',
                    colormap_label: Union[Literal['auto'], str] = 'auto', cluster_columns: bool = True,
                    log_transform: bool = True, z_score_rows: bool = False
                    ) -> plt.Figure:
        """
        Performs hierarchical clustering and plots a clustergram on the base-2 log of a given set of samples.

        :param z_score_rows: if True, the rows will be z-scored before clustering. \
        This will normalize the rows to have a mean of 0 and a standard deviation of 1, such that \
        genes will be clustered based on the similarity of their expression pattern instead of \
        absolute expression levels.
        :type z_score_rows: bool (default=False)
        :param colormap_label: label for the colorbar
        :type colormap_label: str or 'auto' (default='auto')
        :param cluster_columns: if True, both rows and columns will be clustered. Otherwise, \
        only the rows will be clustered, and columns will maintain their original order.
        :type cluster_columns: bool (default=True)
        :param colormap: the colormap to use in the clustergram.
        :type colormap: str
        :param log_transform: if True, will apply a log transform (log2) to the data before clustering.
        :type log_transform: bool (default=True)
        :type sample_names: 'all' or list.
        :param sample_names: the names of the relevant samples in a list. \
        Example input: ["condition1_rep1", "condition1_rep2", "condition1_rep3", \
        "condition2_rep1", "condition3_rep1", "condition3_rep2"]
        :type metric: 'Euclidean', 'hamming', 'correlation', or any other \
        distance metric available in scipy.spatial.distance.pdist
        :param metric: the distance metric to use in the clustergram. \
        For all possible inputs and their meaning see scipy.spatial.distance.pdist documentation online.
        :type linkage: 'single', 'average', 'complete', 'weighted', 'centroid', 'median' or 'ward'.
        :param linkage: the linkage method to use in the clustergram. \
        For all possible inputs and their meaning see scipy.cluster.hierarchy.linkage documentation online.
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :rtype: A matplotlib Figure.


        .. figure:: /figures/clustergram.png
           :align:   center
           :scale: 40 %

           Example plot of clustergram()

        """
        assert isinstance(metric, str) and isinstance(linkage, str), "Linkage and Metric must be strings!"
        metric = metric.lower()
        linkage = linkage.lower()
        metrics = ['correlation', 'cosine', 'euclidean',
                   'hamming', 'jaccard', 'jensenshannon', 'kulsinski', 'mahalanobis', 'matching', 'minkowski',
                   'rogerstanimoto', 'russellrao', 'sEuclidean', 'sokalmichener', 'sokalsneath', 'sqEuclidean', 'yule']
        linkages = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']
        assert metric in metrics, f"Invalid metric {metric}."
        assert linkage in linkages, f"Invalid linkage {linkage}."

        if sample_names == 'all':
            sample_names = list(self.columns)
        if colormap_label == 'auto':
            colormap_label = r"$\log_2$(Normalized reads + 1)" if log_transform else "Normalized reads"
            if z_score_rows:
                colormap_label += "\nZ-score"

        data = np.log2(self.df[sample_names] + 1) if log_transform else self.df[sample_names]

        print('Calculating clustergram...')
        clustergram = sns.clustermap(data, method=linkage, metric=metric,
                                     cmap=sns.color_palette(colormap, as_cmap=True), yticklabels=False,
                                     cbar_kws=dict(label=colormap_label), col_cluster=cluster_columns,
                                     z_score=0 if z_score_rows else None)

        # set colored borders for colorbar and heatmap
        cbar = clustergram.ax_cbar.get_children()[-1]
        cbar.set_edgecolor('black')
        cbar.set_linewidth(2)
        clustergram.ax_cbar.yaxis.set_label_position('left')
        grid = clustergram.ax_heatmap.get_children()[-1]
        grid.set_edgecolor('black')
        grid.set_linewidth(2)

        if title == 'auto':
            title = f"Clustegram of {self.fname.stem}"
        clustergram.figure.suptitle(title, fontsize=title_fontsize)
        plt.tick_params(axis='both', which='both', labelsize=tick_fontsize)
        try:
            clustergram.figure.set_tight_layout(True)
        except AttributeError:
            pass

        plt.show()
        return clustergram.figure

    @readable_name('Plot expression of specific genes')
    def plot_expression(self, features: Union[List[str], str],
                        samples: Union[param_typing.GroupedColumns, Literal['all']] = 'all',
                        avg_function: Literal['mean', 'median', 'geometric_mean'] = 'mean',
                        spread_function: Literal['sem', 'std', 'gstd', 'gsem', 'iqr', 'range'] = 'sem',
                        bar_colors: ColorList = 'deepskyblue', edge_color: Color = 'black',
                        scatter_color: Color = 'grey',
                        count_unit: str = 'Normalized reads', split_plots: bool = False,
                        jitter: Fraction = 0, group_names: Union[List[str], None] = None,
                        log_scale: bool = False) -> plt.Figure:
        """
        Plot the average expression and spread of the specified features under the specified conditions.
        :type features: str or list of strings
        :param features: the feature/features to plot expression for.
        :param samples: A list of the sample names and/or grouped sample names to be plotted. \
        All specified samples must be present in the CountFilter object. \
        To average multiple replicates of the same condition, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :param avg_function: The function used to calculate the average expression value for each condition.
        :type avg_function: 'mean', 'median', or 'geometric_mean' (default='mean')
        :param spread_function: The function used to calculate the error bars of expression values for each condition.
        :type spread_function: 'sem', 'std', 'gstd', 'gsem', 'iqr', or 'range' (default='sem')
        :param bar_colors: The color or list of colors to use for the bars in the plot.
        :type bar_colors: str or list of color strings (default='deepskyblue')
        :param edge_color: The color of the edges around the bars.
        :type edge_color: str (default='black')
        :param scatter_color: The color of the scatter points representing individual samples.
        :type scatter_color: str (default='grey')
        :type count_unit: str (default='Reads per million')
        :param count_unit: The unit of the count data. Will be displayed in the y axis.
        :type split_plots: bool (default=False)
        :param split_plots: if True, each gene will be plotted in its own Figure. \
        Otherwise, all genes will be plotted in the same Figure.
        :type jitter: float (default=0)
        :param jitter: The amount of jitter to apply to the scatter points. \
        This can help visualize overlapping points.
        :type group_names: list of strings or None (default=None)
        :param group_names: Optionally, specify the names of the groups in the plot. \
        If None, the names of the samples will be used.
        :type log_scale: bool (default=False)
        :param log_scale: If True, the y-axis will be displayed in logarithmic scale.

        .. figure:: /figures/plot_expression.png
           :align:   center
           :scale: 40 %

           Example plot of plot_expression()

        """
        features = parsing.data_to_list(features)
        assert validation.isinstanceiter(features, str), "'features' must be a string or list of strings!"
        for feature in features:
            assert feature in self, f"Supplied feature '{feature}' does not appear in this table. "
        if isinstance(samples, str) and samples.lower() == 'all':
            samples = [[item] for item in self.columns]
        else:
            samples = samples.copy()

        for i in range(len(samples)):
            samples[i] = parsing.data_to_list(samples[i])
            if not validation.isinstanceiter(samples[i], str):
                for j in range(len(samples[i])):
                    samples[i][j] = self.columns[samples[i][j]]
        sample_names = parsing.make_group_names(samples) if group_names is None else group_names
        figs = []
        axes = []
        ylims = []

        if not split_plots:
            g = strategies.SquareStrategy()
            subplots = list(g.get_grid(len(features)))
            plt.close()
            fig = plt.figure()
            fig.tight_layout()
            figs.append(fig)
        for i, feature in enumerate(features):
            if split_plots:
                fig, ax = plt.subplots()
                fig.tight_layout()
                figs.append(fig)
            else:
                ax = fig.add_subplot(subplots[i])
            ax.set_yscale('log' if log_scale else 'linear')
            axes.append(ax)

            feature_df = self.df.filter(pl.first() == feature).drop(cs.first()).transpose(include_header=True)
            if avg_function == 'mean':
                mean = [feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).mean().item() for grp in samples]
            elif avg_function == 'median':
                mean = [feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).median().item() for grp in samples]

            elif avg_function == 'geometric_mean':
                mean = [np.exp(np.mean(np.log(feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).to_numpy())))
                        for grp in samples]
            else:
                raise ValueError(f"Invalid average function {avg_function}.")

            if spread_function == 'sem':
                spread = [
                    0 if len(grp) == 1 else sem(feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).to_numpy())[0]
                    for grp in samples]
            elif spread_function == 'std':
                spread = [0 if len(grp) == 1 else feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).std().item()
                          for grp in samples]
            elif spread_function == 'gstd':
                factors = [
                    1 if len(grp) == 1 else gstd(feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).to_numpy())[
                        0] for grp in samples]
                spread = [[0 if len(grp) == 1 else m - (m / f) for m, f, grp in zip(mean, factors, samples)],
                          [0 if len(grp) == 1 else (m * f) - m for m, f, grp in zip(mean, factors, samples)]]
            elif spread_function == 'gsem':
                factors = [1 if len(grp) == 1 else
                           (np.exp(sem(np.log(feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).to_numpy()))))[
                               0] for grp in samples]
                spread = [[0 if len(grp) == 1 else m - (m / f) for m, f, grp in zip(mean, factors, samples)],
                          [0 if len(grp) == 1 else (m * f) - m for m, f, grp in zip(mean, factors, samples)]]
            elif spread_function == 'iqr':
                spread = [feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).quantile(
                    0.75).item() - feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).quantile(0.25).item() for
                          grp in samples]
            elif spread_function == 'range':
                spread = [
                    [0 if len(grp) == 1 else m - feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).min().item()
                     for m, grp in zip(mean, samples)],
                    [0 if len(grp) == 1 else feature_df.filter(pl.first().is_in(grp)).drop(cs.first()).max().item() - m
                     for m, grp in zip(mean, samples)]]
            else:
                raise ValueError(f"Invalid spread function {spread_function}.")
            points_y = parsing.flatten(
                [parsing.data_to_list(feature_df.filter(pl.first().is_in(grp)).drop(cs.first())) for grp in samples])
            points_x = []
            for j, grouping in enumerate(samples):
                this_group = generic.jitter(len(grouping), jitter) + j
                points_x.extend(parsing.data_to_list(this_group))

            ax.bar(np.arange(len(samples)), mean, yerr=spread, edgecolor=edge_color, width=0.5,
                   color=bar_colors, capsize=6.5, error_kw=dict(capthick=2, lw=2))
            ax.scatter(points_x, points_y, edgecolor=edge_color, facecolor=scatter_color, linewidths=0.9)
            ax.set_xticks(np.arange(len(samples)))
            ax.set_xticklabels(list(sample_names), fontsize=12)
            ax.set_title(feature, fontsize=16)
            plt.ylabel(count_unit, fontsize=14)
            generic.despine(ax)
            ylims.append(ax.get_ylim()[1])
        if not split_plots:
            for ax in axes:
                ax.set_ylim((0.0, max(ylims)))
        plt.show()
        return figs

    @readable_name('Sort table by contribution to a Principal Component (PCA)')
    def sort_by_principal_component(self, component: PositiveInt, ascending: bool = True, power_transform: bool = True,
                                    inplace: bool = True):
        """
         Performs Principal Component Analysis (PCA), and sort the table based on the contribution (loadings) \
         of genes to a specific Principal Component. This type of analysis can help you understand which genes \
         contribute the most to each principal component, particularly using single-list enrichment analysis. .

        :param component: the Principal Component the table should be sorted by.
        :type component: positive int
        :type ascending: bool (default=Trle)
        :param ascending: Sort order: ascending (negative loadings at the top of the list) \
         versus descending (positive loadings at the top of the list).
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to standartization and principal component analysis.
        :type power_transform: bool (default=True)
        :type inplace: bool (default=True)
        :param inplace: If True, perform the operation in-place. \
        Otherwise, returns a sorted copy of the Filter object without modifying the original.
        :return: None if inplace=True, a sorted Filter object otherwise.
        """
        assert 0 < component <= self.shape[0], \
            "'component' must be larger than 0 and equal or lower than the number of genes in the table!"
        self._validate_is_normalized()
        suffix = f'_sortbyPC{component}' + 'powertransform' * bool(power_transform)
        data = self.df[self._numeric_columns].to_numpy().transpose()
        data_standardized = generic.standard_box_cox(data) if power_transform else generic.standardize(data)
        pca_obj = PCA(component)
        pca_obj.fit(data_standardized)
        loading = self.df.select(pl.first()).with_columns(
            pl.DataFrame(pca_obj.components_.T[:, component - 1], schema=['sortOrder']))
        new_df = self.df.join(loading, on=self.df.columns[0], how='left')
        new_df = new_df.sort(pl.col('sortOrder'), descending=not ascending).drop(cs.by_name('sortOrder'))
        return self._inplace(new_df, False, inplace, suffix, 'sort')

    @readable_name('Split table by contribution to Principal Components (PCA)')
    def split_by_principal_components(self, components: Union[PositiveInt, List[PositiveInt]],
                                      gene_fraction: param_typing.Fraction = 0.1, power_transform: bool = True
                                      ) -> Union[
        Tuple['CountFilter', 'CountFilter'], Tuple[Tuple['CountFilter', 'CountFilter'], ...]]:
        """
         Performs Principal Component Analysis (PCA), and split the table based on the contribution (loadings) \
         of genes to specific Principal Components. For each Principal Component specified, *RNAlysis* will find the \
         X% most influential genes on the Principal Component based on their loadings (where X is gene_fraction), \
         (X/2)% from the top and (X/2)% from the bottom. This type of analysis can help you understand which genes \
         contribute the most to each principal component.

        :param components: the Principal Components the table should be filtered by. \
        Each Principal Component will be analyzed separately.
        :type components: int or list of integers
        :param gene_fraction: the total fraction of top influential genes that will be returned. \
        For example, if gene_fraction=0.1, *RNAlysis* will return the top and bottom 5% of genes \
        based on their loadings for any principal component.
        :type gene_fraction: float between 0 and 1 (default=0.1)
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to standartization and principal component analysis.
        :type power_transform: bool (default=True)
        """
        components = parsing.data_to_list(components)
        assert validation.isinstanceiter(components, int), "'components' must be a list of integers!"
        assert all([0 < component <= self.shape[0] for component in components]), \
            "'components' must be larger than 0 and equal or lower than the number of genes in the table!"
        assert isinstance(gene_fraction, float) and 0 <= gene_fraction <= 1, \
            "'gene_fraction' must be a number between 0 and 1!"
        self._validate_is_normalized()

        n_components = max(components)
        n_genes = int(np.floor(gene_fraction * self.shape[0] * 0.5))
        print(f'Extracting the top {n_genes} and bottom {n_genes} genes '
              f'contributing to each specified Principal Component...')

        data = self.df[self._numeric_columns].transpose()
        data_standardized = generic.standard_box_cox(data) if power_transform else generic.standardize(data)

        pca_obj = PCA(n_components)
        pca_obj.fit(data_standardized)
        loadings = self.df.select(pl.first()).with_columns(
            pl.DataFrame(pca_obj.components_.T, schema=[f"{i + 1}" for i in range(n_components)]))
        filt_obj_tuples = []
        for component in components:
            loading = loadings.select(pl.first(), pl.col(str(component))).sort(pl.col(str(component)))
            top = self._inplace(self.df.filter(pl.first().is_in(loading.tail(n_genes).select(pl.first()))),
                                opposite=False, inplace=False, suffix=f'_top{gene_fraction}PC{component}')
            bottom = self._inplace(self.df.filter(pl.first().is_in(loading.head(n_genes).select(pl.first()))),
                                   opposite=False, inplace=False, suffix=f'_bottom{gene_fraction}PC{component}')
            filt_obj_tuples.append((top, bottom))

        return filt_obj_tuples[0] if len(filt_obj_tuples) == 1 else tuple(filt_obj_tuples)

    @readable_name('Principal Component Analysis (PCA) plot')
    def pca(self, samples: Union[param_typing.GroupedColumns, Literal['all']] = 'all', n_components: PositiveInt = 3,
            power_transform: bool = True, labels: bool = True, title: Union[str, Literal['auto']] = 'auto',
            title_fontsize: float = 20, label_fontsize: float = 16, tick_fontsize: float = 12,
            proportional_axes: bool = False, plot_grid: bool = True, legend: Union[List[str], None] = None) -> Tuple[
        PCA, List[plt.Figure]]:
        """
        Performs Principal Component Analysis (PCA), visualizing the principal components that explain the most\
        variance between the different samples. The function will standardize the data prior to PCA, and then plot \
        the requested number of pairwise PCA projections.

        :type samples: 'all' or list.
        :param samples: A list of the sample names and/or grouped sample names to be plotted. \
        All specified samples must be present in the CountFilter object. \
        To draw multiple replicates of the same condition in the same color, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :param n_components: number of Principal Components to plot (minimum is 2). *RNAlysis* will generate a \
        pair-wise scatter plot between every pair of Principal Components.
        :type n_components: int >=2 (default=3)
        :param labels: if True, *RNAlysis* will display labels with the sample names next to each sample on the graph.
        :type labels: bool (default=True)
        :param power_transform: if True, RNAlysis will apply a power transform (Box-Cox) \
        to the data prior to standartization and principal component analysis.
        :type power_transform: bool (default=True)
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
        :param tick_fontsize: determines the font size of the X and Y tick labels, and the sample name labels. .
        :type tick_fontsize: float (default=10)
        :param proportional_axes: if True, the dimensions of the PCA plots will be proportional \
        to the percentage of variance explained by each principal component.
        :type proportional_axes: bool (default=False)
        :param plot_grid: if True, will draw a grid on the PCA plot.
        :type plot_grid: bool (default=True)
        :param legend: if enabled, display a legend on the PCA plot. Each entry in the 'legend' parameter \
        corresponds to one group of samples (one color on the graph), as defined by the parameter 'samples'
        :type legend: list of str, or None (default=None)
        :return: A tuple whose first element is an sklearn.decomposition.pca object, \
        and second element is a list of matplotlib.axis objects.

        .. figure:: /figures/pca.png
           :align:   center
           :scale: 40 %

           Example plot of pca()

        """
        assert isinstance(n_components, int) and n_components >= 2, \
            f"'n_components' must be an integer >=2. Instead got {n_components}."
        self._validate_is_normalized()

        if samples == 'all':
            samples = [[col] for col in self._numeric_columns]
        data = self.df.select(pl.col(parsing.flatten(samples))).transpose()
        data_standardized = generic.standard_box_cox(data) if power_transform else generic.standardize(data)

        pca_obj = PCA()
        pcomps = pca_obj.fit_transform(data_standardized)
        columns = [f'Principal component {i + 1}' for i in range(n_components)]
        principal_df = pl.DataFrame(data=pcomps[:, 0:n_components], schema=columns)
        final_df = principal_df.with_columns(pl.Series(parsing.flatten(samples)).alias('lib'))

        if title == 'auto':
            title = f'PCA plot of {self.fname.stem}'
        pc_var = pca_obj.explained_variance_ratio_

        figs = [self._scree_plot(pc_var, label_fontsize, title_fontsize, tick_fontsize)]
        for first_pc in range(n_components):
            for second_pc in range(first_pc + 1, n_components):
                figs.append(CountFilter._pca_plot(
                    final_df=final_df[
                        [f'Principal component {1 + first_pc}', f'Principal component {1 + second_pc}', 'lib']],
                    pc1_var=pc_var[first_pc], pc2_var=pc_var[second_pc], sample_grouping=samples, labels=labels,
                    label_fontsize=label_fontsize, title=title, title_fontsize=title_fontsize, legend=legend,
                    tick_fontsize=tick_fontsize, proportional_axes=proportional_axes, plot_grid=plot_grid))

        return pca_obj, figs

    @staticmethod
    def _scree_plot(pc_var: List[float], label_fontsize: float, title_fontsize: float,
                    tick_fontsize: float) -> plt.Figure:
        pc_var_percent = [var * 100 for var in pc_var]
        pc_var_cumsum = np.cumsum(pc_var_percent)
        pc_indices = list(range(len(pc_var)))
        pc_names = [f'PC{i + 1}' for i in pc_indices]

        fig = plt.figure(figsize=(9, 9), constrained_layout=True)
        ax = fig.add_subplot(1, 1, 1)
        ax.bar(pc_names, pc_var_percent, edgecolor='black')
        ax.plot(pc_names, pc_var_cumsum, '-o', color='black')
        ax.set_title('Principal Component Analysis - Variance Explained', fontsize=title_fontsize)
        ax.set_ylabel('% Variance explained', fontsize=label_fontsize)
        ax.set_xlabel('Principal Components', fontsize=label_fontsize)
        ax.set_ylim([0, 105])
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)
        ax.yaxis.set_major_formatter('{x}%')
        generic.despine(ax)
        return fig

    @staticmethod
    def _pca_plot(final_df: pl.DataFrame, pc1_var: float, pc2_var: float, sample_grouping: param_typing.GroupedColumns,
                  labels: bool, title: str, title_fontsize: float, label_fontsize: float, tick_fontsize: float,
                  proportional_axes: bool, plot_grid: bool, legend: Union[List[str], None]) -> plt.Figure:
        """
        Internal method, used to plot the results from CountFilter.pca().

        :param final_df: The DataFrame output from pca
        :param pc1_var: Variance explained by the first PC.
        :param pc2_var: Variance explained by the second PC.
        :param sample_grouping: A list of the sample names and/or grouped sample names to be plotted. \
        All specified samples must be present in the DataFrame object. \
        To draw multiple replicates of the same condition in the same color, they can be grouped in an inner list. \
        Example input: \
        [['SAMPLE1A', 'SAMPLE1B', 'SAMPLE1C'], ['SAMPLE2A', 'SAMPLE2B', 'SAMPLE2C'],'SAMPLE3' , 'SAMPLE6']
        :return: an axis object containing the PCA plot.

        """
        if proportional_axes:
            dims = (15, 15 * (pc2_var / pc1_var))
        else:
            dims = (8, 8)
        fig = plt.figure(figsize=dims, constrained_layout=True)
        ax = fig.add_subplot(1, 1, 1)
        if proportional_axes:
            ax.set_aspect(pc2_var / pc1_var)

        ax.set_xlabel(f'{final_df.columns[0]} (explained {pc1_var * 100 :.2f}%)', fontsize=label_fontsize)
        ax.set_ylabel(f'{final_df.columns[1]} (explained {pc2_var * 100 :.2f}%)', fontsize=label_fontsize)
        ax.set_title(title, fontsize=title_fontsize)

        color_generator = generic.color_generator()
        colors = [next(color_generator) for _ in range(len(sample_grouping))]
        colors_per_row = parsing.flatten(
            [[colors[i]] * len(parsing.data_to_list(grp)) for i, grp in enumerate(sample_grouping)])
        sample_names = final_df['lib'].to_list()

        for i, grp in enumerate(sample_grouping):
            label = f'Group {i + 1}' if (legend is None or i >= len(legend)) else legend[i]

            grp = parsing.data_to_list(grp)
            item_inds = []
            for item in grp:
                item_inds.append(sample_names.index(item))
            ax.scatter(final_df[item_inds, 0], final_df[item_inds, 1], color=colors[i], label=label, s=75)

        if legend is not None:
            ax.legend(title="Legend", draggable=True)

        if labels:
            for i, (row) in enumerate(final_df.iter_rows()):
                ax.annotate(row[2], (row[0], row[1]), textcoords='offset pixels',
                            xytext=(label_fontsize, label_fontsize),
                            fontsize=label_fontsize, color=colors_per_row[i])

        ax.grid(plot_grid)
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)
        plt.show()
        return fig

    @readable_name('Scatter plot - sample VS sample')
    def scatter_sample_vs_sample(self, sample1: param_typing.ColumnNames, sample2: param_typing.ColumnNames,
                                 xlabel: Union[str, Literal['auto']] = 'auto',
                                 ylabel: Union[str, Literal['auto']] = 'auto',
                                 title: Union[str, Literal['auto']] = 'auto', title_fontsize: float = 20,
                                 label_fontsize: float = 16, tick_fontsize: float = 12, annotation_fontsize: float = 10,
                                 highlight: Union[Sequence[str], None] = None,
                                 point_color: param_typing.Color = '#6d7178',
                                 highlight_color: param_typing.Color = '#00aaff', opacity: param_typing.Fraction = 0.65,
                                 point_size: float = 10, interactive: bool = True, show_cursor: bool = False
                                 ) -> plt.Figure:
        """
        Generate a scatter plot where every dot is a feature, the x value is log10 of reads \
        (counts, RPM, RPKM, TPM, etc) in sample1, the y value is log10 of reads in sample2. \
        If the plot is generated in interactive mode, data points can be labeled by clicking on them.

        :param sample1: Name of the first sample from the CountFilter object. \
        If sample1 is a list, they will be avarged as replicates.
        :type sample1: string or list of strings
        :param sample2: Name of the second sample from the CountFilter object. \
        If sample2 is a list, they will be averaged as replicates.
        :type sample2: string or list of strings
        :param xlabel: optional. If not specified, sample1 will be used as xlabel.
        :type xlabel: str or 'auto'
        :param ylabel: optional. If not specified, sample2 will be used as ylabel.
        :type ylabel: str or 'auto'
        :param title: optional. If not specified, a title will be generated automatically.
        :type title: str or 'auto'

        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
        :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :param annotation_fontsize: determines the font size of the point annotations created in interactive mode.
        :type annotation_fontsize: float (default=10)
        :param highlight: If specified, the points in the scatter corresponding to the names/features in 'highlight' \
        will be highlighted in red.
        :type highlight: Filter object or iterable of strings
        :param point_color: color of the points in the scatter plot.
        :type point_color: str or tuple of (int, int, int) (default='#6d7178')
        :param highlight_color: color of the highlighted points in the scatter plot.
        :type highlight_color: str or tuple of (int, int, int) (default='#00aaff')
        :param opacity: float between 0 and 1 (default=0.65)
        :type opacity: determines the opacity of the points in the scatter plot. 0 indicates completely transparent, \
        while 1 indicates completely opaque.
        :param point_size: determines the size of the points in the scatter plot
        :type point_size: float (default=10)
        :param interactive: if True, turns on interactive mode. While in interactive mode, you can click on a data \
        point to label it with its gene name/ID, or click on a labeled data point to unlabel it.
        :type interactive: bool (default=True)
        :param show_cursor: if True, show the cursor position on the plot during interactive mode
        :type show_cursor: bool (default=False)
        :return: a matplotlib axis object.

        .. figure:: /figures/rpm_vs_rpm.png
           :align:   center
           :scale: 60 %

           Example plot of scatter_sample_vs_sample()


        """
        self._validate_is_normalized()
        sample1, sample2 = parsing.data_to_list(sample1), parsing.data_to_list(sample2)

        xvals = self.df.select(pl.first()).with_columns(np.log10(self.df.select(pl.col(sample1)).mean_horizontal() + 1))
        yvals = self.df.select(pl.first()).with_columns(np.log10(self.df.select(pl.col(sample2)).mean_horizontal() + 1))
        if xlabel.lower() == 'auto':
            xlabel = r'$\log_{10}$(normalized reads + 1) from ' + f'{sample1}'
        if ylabel.lower() == 'auto':
            ylabel = r'$\log_{10}$(normalized reads + 1) from ' + f'{sample2}'
        if title.lower() == 'auto':
            title = f'{sample1} vs {sample2}'

        if interactive:
            fig = plt.figure(figsize=(8, 8), constrained_layout=True, FigureClass=generic.InteractiveScatterFigure,
                             labels=parsing.data_to_tuple(self.df.select(pl.first())),
                             annotation_fontsize=annotation_fontsize,
                             show_cursor=show_cursor)
            ax = fig.axes[0]
            kwargs = {'picker': 5}
        else:
            fig = plt.figure(figsize=(8, 8), constrained_layout=True)
            ax = fig.add_subplot(1, 1, 1)
            kwargs = {}

        generic.despine(ax)
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.set_title(title, fontsize=title_fontsize)
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)

        if highlight is not None:
            highlight_features = highlight.index_set if validation.isinstanceinh(highlight,
                                                                                 Filter) else parsing.data_to_set(
                highlight)
            highlight_valid = parsing.data_to_list(highlight_features.intersection(self.index_set))
            if len(highlight_valid) < len(highlight_features):
                warnings.warn(
                    f'Out of {len(highlight_features)} features to be highlighted, '
                    f'{len(highlight_features) - len(highlight_valid)} features are missing from the '
                    f'CountFilter object and will not be highlighted.')

        else:
            highlight_valid = []

        base_set = parsing.data_to_list(self.index_set.difference(highlight_valid))
        xvals_base = xvals.filter(pl.first().is_in(base_set)).drop(cs.first())
        yvals_base = yvals.filter(pl.first().is_in(base_set)).drop(cs.first())
        xvals_highlight = xvals.filter(pl.first().is_in(highlight_valid)).drop(cs.first())
        yvals_highlight = yvals.filter(pl.first().is_in(highlight_valid)).drop(cs.first())
        ax.scatter(xvals_base, yvals_base, s=point_size, alpha=opacity, c=point_color, **kwargs)
        ax.scatter(xvals_highlight, yvals_highlight, s=point_size, alpha=opacity, c=highlight_color, **kwargs)
        plt.plot([0, 5], [0, 5], c='green')

        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        axis_lims = (0, max(xlim[1], ylim[1]))
        ax.set_xlim(axis_lims)
        ax.set_ylim(axis_lims)

        plt.show()
        return fig

    @readable_name('Box plot')
    def box_plot(self, samples: Union[param_typing.GroupedColumns, Literal['all']] = 'all', notch: bool = True,
                 scatter: bool = False, ylabel: str = 'log10(Normalized reads + 1)',
                 title: Union[str, Literal['auto']] = 'auto',
                 title_fontsize: float = 20, label_fontsize: float = 16, tick_fontsize: float = 12) -> plt.Figure:
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
        :type notch: bool (default=True)
        :param notch: if True, adds a confidence-interval notch to the box-plot.
        :type scatter: bool (default=False)
        :param scatter: if True, adds a scatter-plot on top of the box-plot.
        :type ylabel: str (default='Log10(Normalized reads + 1)')
        :param ylabel: the label of the Y axis.
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
         :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :rtype: a matplotlib Figure.

        .. figure:: /figures/box_plot.png
           :align:   center
           :scale: 60 %

           Example plot of box_plot(add_scatter=True)


        """
        self._validate_is_normalized()
        if samples == 'all':
            samples_df = self.df.select(pl.col(self._numeric_columns))
        else:
            samples_df = self._avg_subsamples(samples).drop(cs.first())

        samples_df = np.log10(samples_df + 1)
        _ = plt.figure(figsize=(8, 8))

        color_gen = generic.color_generator()
        palette = sns.color_palette([next(color_gen) for _ in range(samples_df.shape[1])], n_colors=samples_df.shape[1])
        ax = sns.boxplot(data=samples_df, notch=notch, palette=palette)
        if scatter:
            _ = sns.stripplot(data=samples_df, palette='dark:#404040', size=3)

        if title == 'auto':
            title = f"Box plot of {self.fname.stem}"
        ax.set_title(title, fontsize=title_fontsize)
        ax.set_xlabel("Samples", fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)

        generic.despine(ax)
        plt.show()
        return ax.figure

    @readable_name('Enhanced box plot')
    def enhanced_box_plot(self, samples: Union[param_typing.GroupedColumns, Literal['all']] = 'all',
                          scatter: bool = False, ylabel: str = 'log10(Normalized reads + 1)',
                          title: Union[str, Literal['auto']] = 'auto', title_fontsize: float = 20,
                          label_fontsize: float = 16, tick_fontsize: float = 12) -> plt.Figure:
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
        :type scatter: bool (default=False)
        :param scatter: if True, adds a scatter-plot on top of the box-plot.
        :type ylabel: str (default='Log10(RPM + 1)')
        :param ylabel: the label of the Y axis.
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
        :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :rtype: matplotlib Figure.

        .. figure:: /figures/enhanced_box_plot.png
           :align:   center
           :scale: 60 %

           Example plot of enhanced_box_plot(add_scatter=True)

        """
        self._validate_is_normalized()
        if samples == 'all':
            samples_df = self.df.drop(cs.first())
        else:
            samples_df = self._avg_subsamples(samples).drop(cs.first())

        samples_df = np.log10(samples_df + 1)
        _ = plt.figure(figsize=(8, 8))

        color_gen = generic.color_generator()
        palette = sns.color_palette([next(color_gen) for _ in range(samples_df.shape[1])], n_colors=samples_df.shape[1])

        ax = sns.boxenplot(data=samples_df, palette=palette)
        if scatter:
            _ = sns.stripplot(data=samples_df, palette='dark:#404040', size=3)
        if title == 'auto':
            title = f"Enhanced Box plot of {self.fname.stem}"
        ax.set_title(title, fontsize=title_fontsize)
        ax.set_xlabel("Samples", fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)
        generic.despine(ax)
        plt.show()
        return ax.figure

    @readable_name('Violin plot')
    def violin_plot(self, samples: Union[param_typing.GroupedColumns, Literal['all']] = 'all',
                    ylabel: str = r'$\log_10$(normalized reads + 1)', title: Union[str, Literal['auto']] = 'auto',
                    title_fontsize: float = 20, label_fontsize: float = 16, tick_fontsize: float = 12) -> plt.Figure:
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
        :type ylabel: str (default='log10(normalized reads + 1)')
        :param ylabel: the label of the Y axis.
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
         :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :rtype: a matplotlib Figure.

        .. figure:: /figures/violin.png
           :align:   center
           :scale: 60 %

           Example plot of violin_plot()


        """
        self._validate_is_normalized()
        if samples == 'all':
            samples_df = self.df.select(pl.col(self._numeric_columns))
        else:
            samples_df = self._avg_subsamples(samples).drop(cs.first())

        samples_df = np.log10(samples_df + 1)
        _ = plt.figure(figsize=(8, 8))
        color_gen = generic.color_generator()
        palette = sns.color_palette([next(color_gen) for _ in range(samples_df.shape[1])], n_colors=samples_df.shape[1])

        ax = sns.violinplot(data=samples_df, palette=palette)
        if title == 'auto':
            title = f"Enhanced Box plot of {self.fname.stem}"
        ax.set_title(title, fontsize=title_fontsize)
        ax.set_xlabel("Samples", fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)
        generic.despine(ax)
        plt.show()
        return ax.figure

    @classmethod
    def from_folder(cls, folder_path: str, save_csv: bool = False, fname: str = None, input_format: str = '.txt'
                    ) -> 'CountFilter':
        """
        Iterates over count .txt files in a given folder and combines them into a single CountFilter table. \
        Can also save the count data table and the uncounted data table to .csv files.

        :param folder_path: str or pathlib.Path. Full path of the folder that contains individual htcount .txt files.
        :param save_csv: bool. If True, the joint DataFrame of count data and uncounted data will be saved \
        to two separate .csv files. The files will be saved in 'folder_path', and named according to the parameters \
        'counted_fname' for the count data, and 'uncounted_fname' for the uncounted data (unaligned, \
        alignment not unique, etc).
        :param fname: str. Name under which to save the combined count data table. Does not need to include \
        the '.csv' suffix.
        :param input_format: the file format of the input files. Default is '.txt'.
        :return: an CountFilter object containing the combined count data from all individual htcount .txt files in the \
        specified folder.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder')
            """
        file_suffix = '.csv'
        if save_csv:
            assert isinstance(fname, str)

            if not fname.endswith(file_suffix):
                fname += file_suffix

            counted_fname = os.path.join(folder_path, fname)

        folder = Path(folder_path)
        counts = pl.DataFrame()
        for item in sorted(folder.iterdir()):
            if item.is_file() and item.suffix == input_format:
                counts = counts.with_columns(pl.read_csv(item, separator='\t').select(pl.last().alias(item.stem)))
        assert len(counts) > 0, f"No valid files with the suffix '{input_format}' were found in '{folder_path}'."

        if save_csv:
            io.save_table(df=counts, filename=counted_fname)

        fname = counted_fname if save_csv else os.path.join(folder.absolute(), folder.name + file_suffix)
        count_filter_obj = cls.from_dataframe(counts, Path(fname))
        return count_filter_obj

    @classmethod
    def from_folder_htseqcount(cls, folder_path: str, norm_to_rpm: bool = False, save_csv: bool = False,
                               counted_fname: str = None, uncounted_fname: str = None,
                               input_format: str = '.txt') -> 'CountFilter':

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
                >>> c = filtering.CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder')

                >>> c = filtering.CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', norm_to_rpm=True) # This will also normalize the CountFilter to reads-per-million (RPM).
               Normalized 10 features. Normalized inplace.

                >>> c = filtering.CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', save_csv=True, counted_fname='name_for_reads_csv_file', uncounted_fname='name_for_uncounted_reads_csv_file') # This will also save the counted reads and uncounted reads as separate .csv files

            """
        file_suffix = '.csv'
        misc_names = ['__no_feature', '__ambiguous', '__alignment_not_unique', '__too_low_aQual', '__not_aligned']
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
        df = pl.DataFrame()
        for item in sorted(folder.iterdir()):
            if item.is_file() and item.suffix == input_format:
                this_df = pl.read_csv(item, separator='\t', has_header=False)
                if len(df) == 0:
                    df = this_df.rename({this_df.columns[1]: item.stem, this_df.columns[0]: ''})
                    continue
                df = df.with_columns(this_df.select(pl.nth(-1).alias(item.stem)))
        assert df.shape[1] > 0, f"Error: no valid files with the suffix '{input_format}' were found in '{folder_path}'."

        uncounted = df.filter(pl.first().is_in(misc_names))
        counts = df.filter(~pl.first().is_in(misc_names))

        if save_csv:
            io.save_table(df=counts, filename=counted_fname)
            io.save_table(df=uncounted, filename=uncounted_fname)

        fname = counted_fname if save_csv else os.path.join(folder.absolute(), folder.name + file_suffix)
        count_filter_obj = cls.from_dataframe(counts, Path(fname))
        if norm_to_rpm:
            count_filter_obj.normalize_to_rpm_htseqcount(uncounted)
        return count_filter_obj


class Pipeline(generic.GenericPipeline):
    """
    A collection of functions to be applied sequentially to Filter objects.


    **Attributes**

    functions: list
        A list of the functions in the Pipeline.
    params: list
        A list of the parameters of the functions in the Pipeline.
    filter_type: Filter object
        The type of Filter objects to which the Pipeline can be applied
    """
    __slots__ = {'filter_type': 'type of filter objects to which the Pipeline can be applied'}
    FILTER_TYPES = {'filter': Filter, 'deseqfilter': DESeqFilter, 'foldchangefilter': FoldChangeFilter,
                    'countfilter': CountFilter}
    FILTER_TYPES_REV = {val: key for key, val in FILTER_TYPES.items()}

    def __init__(self, filter_type: Union[str, 'Filter', 'CountFilter', 'DESeqFilter', 'FoldChangeFilter'] = Filter):
        """
        :param filter_type: the type of Filter object the Pipeline can be applied to.
        :type filter_type: str or Filter object (default=filtering.Filter)

        :Examples:
            >>> from rnalysis import filtering
            >>> pipe = filtering.Pipeline()
            >>> deseq_pipe = filtering.Pipeline('deseqfilter')

        """
        assert isinstance(filter_type,
                          (type, str)), f"'filter_type' must be type of a Filter object, is instead {type(filter_type)}"
        if isinstance(filter_type, str):
            assert filter_type.lower() in self.FILTER_TYPES, f"Invalid filter_type {filter_type}. "
            filter_type = self.FILTER_TYPES[filter_type.lower()]
        else:
            assert filter_type in self.FILTER_TYPES.values(), f"Invalid filter_type {filter_type}"
        self.filter_type = filter_type
        super().__init__()

    def __str__(self):
        string = f"Pipeline for {self.filter_type.readable_name}"
        if len(self) > 0:
            string += ":\n\t" + '\n\t'.join(
                self._readable_func_signature(func, params[0], params[1]) for func, params in
                zip(self.functions, self.params))
        return string

    def __repr__(self):
        string = f"Pipeline('{self.filter_type.__name__}')"
        if len(self) > 0:
            string += ": " + "-->".join(
                self._func_signature(func, params[0], params[1]) for func, params in zip(self.functions, self.params))
        return string

    def __eq__(self, other):
        if self.filter_type == other.filter_type and super().__eq__(other):
            return True
        return False

    def _get_pipeline_dict(self):
        d = super()._get_pipeline_dict()
        d['filter_type'] = self.FILTER_TYPES_REV[self.filter_type]
        return d

    def _init_from_dict(self, pipeline_dict: dict):
        self.filter_type = self.FILTER_TYPES[pipeline_dict['filter_type']]
        self.params = [(parsing.data_to_tuple(p[0]), p[1]) for p in pipeline_dict['params']]
        self.functions = [getattr(self.filter_type, func) for func in pipeline_dict['functions']]

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
        return f"{self.filter_type.__name__}.{super()._func_signature(func, args, kwargs)}"

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

        :Examples:
            >>> from rnalysis import filtering
            >>> pipe = filtering.Pipeline()
            >>> pipe.add_function(filtering.Filter.filter_missing_values)
            Added function 'Filter.filter_missing_values()' to the pipeline.
            >>> pipe.add_function('number_filters', 'col1', 'greater than', value=5, opposite=True)
            Added function 'Filter.number_filters('col1', 'greater than', value=5, opposite=True)' to the pipeline.

        """
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
        super().add_function(func, *args, **kwargs)

    def _apply_filter_norm_sort(self, func: types.FunctionType,
                                filter_object: Union['Filter', 'CountFilter', 'DESeqFilter', 'FoldChangeFilter'],
                                args: tuple, kwargs: dict, inplace: bool, other_outputs: dict, other_cnt):
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
                    temp_res = [
                        self._apply_filter_norm_sort(func, this_obj, args, kwargs, inplace, other_outputs, other_cnt)
                        for this_obj in filter_object]
                    temp_object = []
                    for item in temp_res:
                        if validation.isinstanceinh(item, Filter) or isinstance(item, tuple):
                            temp_object.append(item)
                        elif item is not None:
                            other_outputs[f"{func.__name__}_{other_cnt.setdefault(func.__name__, 1)}"] = item
                            other_cnt[func.__name__] += 1
                    filter_object = parsing.data_to_tuple(temp_object)
                else:
                    filter_object = func(filter_object, *args, **kwargs)
            except (ValueError, AssertionError, TypeError) as e:
                raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")
        else:
            kwargs['inplace'] = True
            try:
                if isinstance(filter_object, tuple):
                    _ = [self._apply_filter_norm_sort(func, this_obj, args, kwargs, inplace, other_outputs, other_cnt)
                         for this_obj in filter_object]
                else:
                    res = func(filter_object, *args, **kwargs)
                    if res is not None:
                        other_outputs[f"{func.__name__}_{other_cnt.setdefault(func.__name__, 1)}"] = res
                        other_cnt[func.__name__] += 1
            except (ValueError, AssertionError, TypeError) as e:
                raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")
        return filter_object

    def _apply_split(self, func: types.FunctionType,
                     filter_object: Union['Filter', 'CountFilter', 'DESeqFilter', 'FoldChangeFilter'],
                     args: tuple, kwargs: dict, other_outputs: dict, other_cnt: dict):
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
                temp_object = [item for item in
                               [self._apply_split(func, this_obj, args, kwargs, other_outputs, other_cnt) for this_obj
                                in filter_object] if item is not None]
                filter_object = parsing.data_to_tuple(temp_object)

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
                    filter_object = parsing.data_to_tuple(temp_object)
                else:
                    raise ValueError(f"Unrecognized output type {type(temp_res)} from function {func}:\n{temp_res}")
        except (ValueError, AssertionError, TypeError) as e:
            raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")
        return filter_object

    def _apply_other(self, func: types.FunctionType,
                     filter_object: Union['Filter', 'CountFilter', 'DESeqFilter', 'FoldChangeFilter'],
                     args: tuple, kwargs: dict, other_outputs: dict, other_cnt: dict, recursive_call: bool = False):
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
        key = f"{func.__name__}_{other_cnt.setdefault(func.__name__, 1)}"
        try:
            if isinstance(filter_object, tuple):
                for this_obj in filter_object:
                    self._apply_other(func, this_obj, args, kwargs, other_outputs, other_cnt, True)
                if len(other_outputs[key]) > 0 and not recursive_call:
                    other_outputs[key] = parsing.data_to_tuple(other_outputs[key])
                    other_cnt[func.__name__] += 1
            else:
                this_output = func(filter_object, *args, **kwargs)
                if this_output is not None:
                    if not recursive_call:
                        other_outputs[key] = this_output
                        other_cnt[func.__name__] += 1
                    else:
                        if key not in other_outputs:
                            other_outputs[key] = []
                        other_outputs[key].append(this_output)
        except (ValueError, AssertionError, TypeError) as e:
            raise e.__class__(f"Invalid function signature {self._func_signature(func, args, kwargs)}")

    @readable_name('Apply Pipeline to a table')
    def apply_to(self, filter_object: Union['Filter', 'CountFilter', 'DESeqFilter', 'FoldChangeFilter'],
                 inplace: bool = True
                 ) -> Union['Filter', Tuple['Filter', dict], Tuple[Tuple['Filter'], dict], dict, None]:

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

        :Examples:
            >>> from rnalysis import filtering
            >>> # create the pipeline
            >>> pipe = filtering.Pipeline('DESeqFilter')
            >>> pipe.add_function(filtering.DESeqFilter.filter_missing_values)
            Added function 'DESeqFilter.filter_missing_values()' to the pipeline.
            >>> pipe.add_function(filtering.DESeqFilter.filter_top_n, by='padj', n=3)
            Added function 'DESeqFilter.filter_top_n(by='padj', n=3)' to the pipeline.
            >>> pipe.add_function('sort', by='baseMean')
            Added function 'DESeqFilter.sort(by='baseMean')' to the pipeline.
            >>> # load the Filter object
            >>> d = filtering.DESeqFilter('tests/test_files/test_deseq_with_nan.csv')
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

        """
        self._validate_pipeline()
        # noinspection PyTypeHints
        assert issubclass(filter_object.__class__,
                          self.filter_type), f"Supplied filter object of type {type(filter_object)} " \
                                             f"mismatches the specified filter_type {self.filter_type}. "

        original_filter_obj = copy.copy(filter_object)
        other_outputs = dict()
        other_cnt = dict()
        # iterate over all functions and arguments
        for func, (args, kwargs) in zip(self.functions, self.params):
            keywords = ['filter', 'normalize', 'sort', 'translate']
            if any([kw in func.__name__ for kw in keywords]):
                filter_object = self._apply_filter_norm_sort(func, filter_object, args, kwargs, inplace, other_outputs,
                                                             other_cnt)
            elif func.__name__.startswith('split'):
                assert not inplace, f"Cannot apply the split function {self._func_signature(func, args, kwargs)} " \
                                    f"when inplace={inplace}!"
                filter_object = self._apply_split(func, filter_object, args, kwargs, other_outputs, other_cnt)
            else:
                self._apply_other(func, filter_object, args, kwargs, other_outputs, other_cnt)

        if not inplace or isinstance(filter_object, tuple):
            if filter_object != original_filter_obj:
                if len(other_outputs) > 0:
                    return filter_object, other_outputs
                return filter_object
        if len(other_outputs) > 0:
            return other_outputs
