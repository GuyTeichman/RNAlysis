"""
This module can perform enrichment analyses on a given set of genomic features and visualize their intersections. \
These include gene ontology/tissue/phenotype enrichment, enrichment for user-defined attributes, \
set visualization ,etc. \
Results of enrichment analyses can be saved to .csv files.
"""
import warnings
from itertools import compress, repeat
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Set, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib_venn as vn
from numba import jit
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as multitest
import tissue_enrichment_analysis as tea
import upsetplot as upset
from ipyparallel import Client
from matplotlib.cm import ScalarMappable
from scipy.stats import hypergeom, ttest_1samp
from statsmodels.stats.descriptivestats import sign_test
from xlmhg import get_xlmhg_test_result as xlmhg_test

from rnalysis import utils
from rnalysis.filtering import Filter


class FeatureSet:
    """ Receives a filtered gene set and the set's name (optional) and preforms various enrichment analyses on them. """
    __slots__ = {'gene_set': 'set of feature names/indices', 'set_name': 'name of the FeatureSet'}
    _go_dicts = {}

    def __init__(self, gene_set: Union[List[str], Set[str], Filter] = None, set_name: str = ''):

        """
        :param gene_set: the set of genomic features to be used in downstream analyses
        :type gene_set: filtering.Filter object, set of strings or list of strings
        :param set_name: name of the FeatureSet
        :type set_name: str


        :Examples:
            >>> from rnalysis import enrichment, filtering
            >>> my_set = enrichment.FeatureSet({'gene1','gene2','gene2'}, 'name of my set')

            >>> filter_obj = filtering.CountFilter('tests/test_files/counted.csv')
            >>> my_other_set = enrichment.FeatureSet(filter_obj, 'name of my other set')

        """
        assert isinstance(set_name, str), f"'set_name' must be of type str, instead got {type(set_name)}."
        if gene_set is None:
            self.gene_set = utils.from_string(
                "Please insert genomic features/indices separated by newline \n"
                "(example: \n'WBGene00000001\nWBGene00000002\nWBGene00000003')", delimiter='\n')
        elif isinstance(gene_set, set):
            pass
        elif isinstance(gene_set, (list, tuple)):
            gene_set = set(gene_set)
        elif utils.isinstanceinh(gene_set, Filter):
            gene_set = gene_set.index_set
        else:
            raise TypeError(f"Error: 'gene_set' must be a set, list or tuple! Is a {type(gene_set)} instead. ")
        self.gene_set = gene_set
        self.set_name = set_name

    def __repr__(self):
        return f"{self.__class__.__name__}: '{self.set_name}'"

    def __len__(self):
        return len(self.gene_set)

    def __contains__(self, item):
        return True if item in self.gene_set else False

    def change_set_name(self, new_name: str):
        """
        Change the 'set_name' of a FeatureSet to a new name.

        :param new_name: the new set name
        :type new_name: str

        """
        assert isinstance(new_name, str), f"New set name must be of type str. Instead, got {type(new_name)}"
        self.set_name = new_name

    def save_txt(self, fname: Union[str, Path]):

        """
        Save the list of features in the FeatureSet object under the specified filename and path.

        :type fname: str or pathlib.Path
        :param fname: full filename/path for the output file. Can include the '.txt' suffix but doesn't have to.

        """
        assert isinstance(fname, (str, Path)), "fname must be str or pathlib.Path!"
        if isinstance(fname, str):
            if not fname.endswith('.txt'):
                fname = fname + '.txt'
        elif isinstance(fname, Path):
            if not fname.suffix == '.txt':
                fname = Path(f"{str(fname.parent)}{fname.name}.txt")
        with open(fname, 'w') as f:
            for gene in self.gene_set:
                f.write(gene + '\n')

    def _set_ops(self, others, op: Callable):

        """
        Performs a given set operation on self and on another object (FeatureSet or set).
        :type others: FeatureSet, set or str
        :param others: Other object to perform set operation with.
        :type: op: Callable (set.union, set.intersection, set.difference or set.symmetric difference)
        :param op: The set operation to be performed.
        :return: A set resulting from the set operation.

        """
        others = list(others)
        for i, other in enumerate(others):
            if isinstance(other, set):
                pass
            elif isinstance(other, FeatureSet):
                others[i] = other.gene_set
            else:
                raise TypeError("'other' must be an FeatureSet object or a set!")
        try:
            return op(self.gene_set, *others)
        except TypeError as e:
            if op == set.symmetric_difference:
                raise TypeError(
                    f"Symmetric difference can only be calculated for two objects, {len(others) + 1} were given!")
            else:
                raise e

    def union(self, *others):

        """
         Calculates the set union of the indices from multiple FeatureSet objects \
        (the indices that exist in at least one of the FeatureSet objects).

        :type others: FeatureSet, RankedSet or set
        :param others: The objects against which the current object will be compared.
        :return: a new FeatureSet with elements from this FeatureSet and all other objects.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000004','WBGene00000005','WBGene00000006'}, 'set name')
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> en.union(s, en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000003', 'WBGene00000004', 'WBGene00000001', 'WBGene00000002', 'WBGene00000006', 'WBGene00000005'}

        """
        return FeatureSet(self._set_ops(others, set.union))

    def intersection(self, *others):

        """
        Calculates the set intersection of the indices from multiple FeatureSet objects \
        (the indices that exist in ALL of the FeatureSet objects).

        :type others: FeatureSet, RankedSet or set
        :param others: The objects against which the current object will be compared.
        :return: a new FeatureSet with elements common to this FeatureSet and all other objects.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000001','WBGene00000002','WBGene00000006'}, 'set name')
            >>> s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> en.intersection(s, en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000001'}

        """
        return FeatureSet(self._set_ops(others, set.intersection))

    def difference(self, *others):

        """
        Calculates the set difference of the indices from multiple FeatureSet objects \
        (the indices that appear in the first FeatureSet object but NOT in the other objects).

        :type others: FeatureSet, RankedSet or set
        :param others: The objects against which the current object will be compared.
        :return: a new FeatureSet with elements in this FeatureSet that are not in the other objects.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000001','WBGene00000002','WBGene00000006'}, 'set name')
            >>> s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> en.difference(s, en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000006'}

        """
        return FeatureSet(self._set_ops(others, set.difference))

    def symmetric_difference(self, other):

        """
        Calculates the set symmetric difference of the indices from two FeatureSet objects \
        (the indices that appear in EXACTLY ONE of the FeatureSet objects, and not both/neither). \
        A-symmetric difference-B is equivalent to (A-difference-B)-union-(B-difference-A).

        :type other: FeatureSet, RankedSet or set
        :param other: A second object against which the current object will be compared.
        :return: a new FeatureSet with elements in either this FeatureSet or the other object, but not both.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000001','WBGene00000002','WBGene00000006'}, 'set name')
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> en.symmetric_difference(en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000002', 'WBGene00000006', 'WBGene00000004'}

        """
        return FeatureSet(self._set_ops([other], set.symmetric_difference))

    @staticmethod
    def _enrichment_save_csv(df: pd.DataFrame, fname: str):

        """
        Internal method, used to save enrichment results to .csv files. Static class method.

        :param df: pandas DataFrame to be saved.
        :param fname: Name and full path under which the DataFrame will be saved

        """
        if fname is None:
            fname = input("Please insert the full name and path to save the file to")
        else:
            assert isinstance(fname, (str, Path))
        if isinstance(fname, Path):
            fname = str(Path)
        utils.save_csv(df, filename=fname + '.csv')

    def go_enrichment(self, mode: str = 'all', alpha: float = 0.05, save_csv: bool = False, fname: str = None):

        """
        Analyzes GO, Tissue and/or Phenotype enrichment for the given group of genomic features. \
        Uses the the Anatomy, Phenotype and Gene Ontology annotations for C. elegans. \
        Corrected p-values are calculated using hypergeometric statistics. \
        For more details see GitHub page of the developers: https://github.com/dangeles/TissueEnrichmentAnalysis

        :type mode: 'go', 'tissue', 'phenotype' or 'all' (default 'all')
        :param mode: the enrichment you wish to perform. 'go' for gene ontology enrichment, \
        'tissue' for tissue enrichment, 'phenotype' for phenotype enrichment, or 'all' for all three.
        :type alpha: float between 0 and 1 (default 0.05)
        :param alpha: Significance threshold.
        :type save_csv: bool (default False)
        :param save_csv: If True, save the result to a csv.
        :type fname: str or pathlib.Path
        :param fname: Name and path in which to save the results. Must be specified if save_csv is True.
        :return: a DataFrame which contains the significant enrichmenet terms

        .. figure::  go_en.png
           :align:   center
           :scale: 40 %

           Example plot of GO enrichment

        .. figure::  tissue_en.png
           :align:   center
           :scale: 40 %

           Example plot of Tissue enrichment

        """
        assert isinstance(alpha, float), "alpha must be a float!"
        assert isinstance(mode, str), "'mode' must be a string!"
        plt.style.use('seaborn-white')
        if mode == 'all':
            d = []
            df_comb = pd.DataFrame()
            for k, arg in enumerate(('go', 'tissue', 'phenotype')):
                print(f'Calculating... {100 * k / 3 :.2f}% done')
                if arg in FeatureSet._go_dicts:
                    d.append(FeatureSet._go_dicts[arg])
                else:
                    d.append(tea.fetch_dictionary(arg))
                    FeatureSet._go_dicts[arg] = d[-1]
                df = tea.enrichment_analysis(self.gene_set, d[-1], alpha=alpha)
                if not df.empty:
                    df_comb = df_comb.append(df)
                    plt.figure()
                    tea.plot_enrichment_results(df, title=f'{arg.capitalize()} Enrichment Analysis', analysis=arg)
                    plt.title(f'{arg.capitalize()} Enrichment Analysis for sample {self.set_name}', fontsize=20)

        else:
            assert (mode == 'go' or mode == 'tissue' or mode == 'phenotype'), "Invalid mode!"
            d = tea.fetch_dictionary(mode)
            df_comb = tea.enrichment_analysis(self.gene_set, d, show=True)
            if not df_comb.empty:
                tea.plot_enrichment_results(df_comb, title=f'{mode.capitalize()} Enrichment Analysis', analysis=mode)
                plt.title(f'{mode.capitalize()} Enrichment Analysis', fontsize=20)

        if save_csv:
            self._enrichment_save_csv(df_comb, fname)
        plt.show()
        return df_comb

    @staticmethod
    def _single_enrichment(gene_set, attribute, attr_ref_df: pd.DataFrame, reps: int):
        assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
        srs = attr_ref_df[attribute]
        bg_array = np.isnan(srs.values)
        obs_array = np.isnan(srs.loc[gene_set].values)
        n = len(gene_set)
        expected_fraction = (bg_array.shape[0] - np.sum(bg_array)) / bg_array.shape[0]
        observed_fraction = (n - np.sum(obs_array)) / n
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = FeatureSet._calc_randomization_pval(n, log2_fold_enrichment, bg_array, reps, observed_fraction)

        return [attribute, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval]

    @staticmethod
    @jit(nopython=True)
    def _calc_randomization_pval(n: int, log2fc: float, bg_array: np.ndarray, reps: int, obs_frac: float):
        ind_range = np.arange(bg_array.shape[0])
        success = 0
        if log2fc >= 0:
            for _ in range(reps):
                success += (n - np.sum(bg_array[np.random.choice(ind_range, n, replace=False)])) / n >= obs_frac

        else:
            for _ in range(reps):
                success += (n - np.sum(bg_array[np.random.choice(ind_range, n, replace=False)])) / n <= obs_frac
        pval = (success + 1) / (reps + 1)
        return pval

    @staticmethod
    def _enrichment_get_attrs(attributes, attr_ref_path):
        if attributes is None:
            attributes = utils.from_string(
                "Please insert attributes separated by newline "
                "(for example: \n'epigenetic_related_genes\nnrde-3 targets\nALG-3/4 class small RNAs')")
        elif isinstance(attributes, (str, int)):
            attributes = [attributes]
        else:
            assert isinstance(attributes, (list, tuple, set)), "'attributes' must be a list, tuple or set!"
            for attr in attributes:
                if isinstance(attr, int):
                    assert attr >= 0, f"Error in attribute number {attr}: index must be non-negative!"
                else:
                    assert isinstance(attr, str), f"Invalid type of attribute {attr}: {type(attr)}"

        try:
            with open(attr_ref_path) as f:
                all_attrs = f.readline().split(',')[1::]
        except FileNotFoundError:
            raise FileNotFoundError(f"Invalid or nonexistent Attribute Reference Table path! path:'{attr_ref_path}'")
        if all_attrs[-1].endswith('\n'):
            all_attrs[-1] = all_attrs[-1][:-1]

        if attributes == ['all']:
            attributes = all_attrs
        elif np.all([True if isinstance(i, int) else False for i in attributes]):
            return [all_attrs[ind] for ind in attributes]
        return attributes

    def _enrichment_get_reference(self, biotype, background_genes, attr_ref_path, biotype_ref_path):
        gene_set = self.gene_set

        attr_ref_df = utils.load_csv(attr_ref_path)
        utils.attr_table_assertions(attr_ref_df)
        attr_ref_df.set_index('gene', inplace=True)

        assert (isinstance(biotype, (str, list, set, tuple)))

        if background_genes is None:
            pass
        else:
            assert isinstance(background_genes,
                              (set, FeatureSet)) or utils.isinstanceinh(background_genes, Filter), \
                f"background_genes must be a set, enrichment.FeatureSet or filtering.Filter; " \
                f"instead got {type(background_genes)}"
            if isinstance(background_genes, FeatureSet):
                background_genes = background_genes.gene_set
            elif utils.isinstanceinh(background_genes, Filter):
                background_genes = background_genes.index_set
            if biotype != 'all':
                warnings.warn(
                    "both 'biotype' and 'background_genes' were specified. Therefore 'biotype' is ignored. ")
                biotype = 'all'

            attr_ref_df = attr_ref_df.loc[background_genes.intersection(set(attr_ref_df.index))]
            if len(attr_ref_df.index) < len(background_genes):
                warnings.warn(
                    f"{len(background_genes) - len(attr_ref_df.index)} indices from the requested "
                    f"background genes do not appear in the Attribute Reference Table, and are therefore ignored. \n"
                    f"This leaves a total of {len(attr_ref_df.index)} background genes. ")
        if biotype == 'all':
            pass
        else:
            biotype_ref_df = utils.load_csv(biotype_ref_path)
            utils.biotype_table_assertions(biotype_ref_df)
            biotype_ref_df.set_index('gene', inplace=True)
            biotype_ref_df.columns = biotype_ref_df.columns.str.lower()
            if isinstance(biotype, (list, tuple, set)):
                mask = pd.Series(np.zeros_like(biotype_ref_df['biotype'].values, dtype=bool),
                                 biotype_ref_df['biotype'].index,
                                 name='biotype')
                for bio in biotype:
                    mask = mask | (biotype_ref_df['biotype'] == bio)
            else:
                biotype_ref_df = biotype_ref_df.loc[biotype_ref_df.index.intersection(attr_ref_df.index)]
                mask = biotype_ref_df['biotype'] == biotype
            attr_ref_df = attr_ref_df.loc[biotype_ref_df[mask].index]
        print(f"{len(attr_ref_df.index)} background genes are used. ")

        not_in_bg = gene_set.difference(set(attr_ref_df.index))
        if len(not_in_bg) > 0:
            gene_set = gene_set.difference(not_in_bg)
            warnings.warn(f"{len(not_in_bg)} genes in the enrichment set do not appear in the "
                          f"Attribute Reference Table and/or the background genes. \n"
                          f"Enrichment will be run on the remaining {len(gene_set)}.")
        attr_ref_df.sort_index(inplace=True)
        return attr_ref_df, gene_set

    def _enrichment_setup(self, biotype: str, background_genes: Union[Iterable[str], str, Iterable[int], int, None],
                          attr_ref_path: str, biotype_ref_path: str, attributes: Union[str, List[str], List[int]]):
        """
        Perform setup for enrichment functions. This function receives most of the input variables \
        from enrichment functions, generates a full list of attributes for enrichment, gets the relevant full path to \
        Attribute/Biotype reference tables, loads the Reference tables into a DataFrame, \
        filters the Attribute DataFrame to include only relevant Attributes, \
        generates a background gene set according to the user's specifications, \
        filters the tested gene set based on the background gene set and/or biotype, \
        and standardizes the data scale input.

        """
        attr_ref_path = utils.get_attr_ref_path(attr_ref_path)
        biotype_ref_path = utils.get_biotype_ref_path(biotype_ref_path)
        attr_ref_df, gene_set = self._enrichment_get_reference(biotype=biotype, background_genes=background_genes,
                                                               attr_ref_path=attr_ref_path,
                                                               biotype_ref_path=biotype_ref_path)

        attributes = self._enrichment_get_attrs(attributes=attributes, attr_ref_path=attr_ref_path)
        return attr_ref_df, gene_set, attributes

    def _enrichment_output(self, enriched_list: list, fdr: float, save_csv: bool, fname: str, plot: bool,
                           return_fig: bool, plot_horizontal: bool = True):
        """
        Formats the enrich list into a results Dataframe, saves the DataFrame to csv if requested, \
        plots the enrichment results, and returns either the Dataframe alone or the Dataframe and the Figure object.
        Called at the end of every enrichment function \
        (enrich_randomization, enrich_randomization_parallel, enrich_statistic...).

        """
        res_df = pd.DataFrame(enriched_list,
                              columns=['name', 'samples', 'obs', 'exp', 'log2_fold_enrichment', 'pval'])
        res_df.replace(-np.inf, -np.max(np.abs(res_df['log2_fold_enrichment'].values)))
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.set_index('name', inplace=True)
        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        if plot:
            fig = self.plot_enrichment_results(res_df, title=f"Enrichment for {self.set_name}",
                                               plot_horizontal=plot_horizontal)
            if return_fig:
                return res_df, fig
        return res_df

    def enrich_randomization_parallel(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                                      fdr: float = 0.05, reps: int = 10000,
                                      biotype: str = 'protein_coding',
                                      background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                                      attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                                      save_csv: bool = False, fname=None, return_fig: bool = False,
                                      plot_horizontal: bool = True, random_seed: int = None):

        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set \
        using a randomization test running in parallel. The attributes are drawn from an Attribute Reference Table. \
        Parallel processing makes this function generally faster than FeatureSet.enrich_randomization. \
        Results should otherwise be the same between the two functions. \
        To use it you must first start a parallel session, using the function 'general.start_parallel_session()'. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. P-values are calculated using a \
        randomization test with the formula p = (successes + 1)/(repeats + 1). \
        This formula results in a positively-biased estimator of the real p-value \
        (a conservative estimate of p-value). When the number of reps approaches infinity, \
        the formula results in an unbiased estimator of the real p-value. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
        :type reps: int larger than 0
        :param reps: How many repetitions to run the randomization for. \
        10,000 is the default. Recommended 10,000 or higher.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        Cannot be specified together with 'biotype'.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :type random_seed: non-negative integer (default None)
        :type random_seed: The random seed used to initialize the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization_parallel()


        .. figure::  plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization_parallel(plot_horizontal = False)

       """
        attr_ref_df, gene_set, attributes, = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)

        client = Client()
        dview = client[:]
        dview.execute("""import numpy as np
              import pandas as pd""")
        if random_seed is not None:
            assert isinstance(random_seed, int) and random_seed >= 0, f"random_seed must be a non-negative integer. " \
                                                                      f"Value {random_seed} invalid."
            dview.execute(f"np.random.seed({random_seed})")
        k = len(attributes)
        gene_set_rep = list(repeat(gene_set, k))
        attr_ref_df_rep = list(repeat(attr_ref_df, k))
        reps_rep = list(repeat(reps, k))

        res = dview.map(FeatureSet._single_enrichment, gene_set_rep, attributes, attr_ref_df_rep, reps_rep)
        enriched_list = res.result()
        return self._enrichment_output(enriched_list, fdr, save_csv, fname, return_fig, True, plot_horizontal)

    def enrich_randomization(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                             reps: int = 10000, biotype: str = 'protein_coding',
                             background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                             attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                             save_csv: bool = False, fname=None, return_fig: bool = False, plot_horizontal: bool = True,
                             random_seed: int = None):

        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set \
        using a randomization test. The attributes are drawn from an Attribute Reference Table. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. P-values are calculated using a \
        randomization test with the formula p = (successes + 1)/(repeats + 1). \
        This formula results in a positively-biased estimator of the real p-value \
        (a conservative estimate of p-value). When the number of reps approaches infinity, \
        the formula results in an unbiased estimator of the real p-value. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
        :type reps: int larger than 0
        :param reps: How many repetitions to run the randomization for. \
        10,000 is the default. Recommended 10,000 or higher.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        Cannot be specified together with 'biotype'.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :type random_seed: non-negative integer (default None)
        :type random_seed: The random seed used to initialize the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization()


        .. figure::  plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization(plot_horizontal = False)

        """
        attr_ref_df, gene_set, attributes = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)
        enriched_list = []
        if random_seed is not None:
            assert isinstance(random_seed, int) and random_seed >= 0, f"random_seed must be a non-negative integer. " \
                                                                      f"Value {random_seed} invalid."
            np.random.seed(random_seed)

        for k, attribute in enumerate(attributes):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            enriched_list.append(self._single_enrichment(gene_set, attribute, attr_ref_df, reps))
            print(f"Finished {k + 1} attributes out of {len(attributes)}")

        return self._enrichment_output(enriched_list, fdr, save_csv, fname, return_fig, True, plot_horizontal)

    def enrich_hypergeometric(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                              biotype: str = 'protein_coding',
                              background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                              attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                              save_csv: bool = False, fname=None, return_fig: bool = False,
                              plot_horizontal: bool = True):

        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set \
        using the Hypergeometric Test. The attributes are drawn from an Attribute Reference Table. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. \
        P-values are calculated using a hypergeometric test: \
        Given M genes in the background set, n genes in the test set, \
        with N genes from the background set belonging to a specific attribute (or 'success') \
        and X genes from the test set belonging to that attribute. \
        If we were to randomly draw n genes from the background set (without replacement), \
        what is the probability of drawing X or more (in case of enrichment)/X or less (in case of depletion) \
        genes belonging to the given attribute? \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment score = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_hypergeometric()


        .. figure::  plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_hypergeometric(plot_horizontal = False)

        """
        attr_ref_df, gene_set, attributes = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)
        enriched_list = []
        for k, attribute in enumerate(attributes):
            srs = attr_ref_df[attribute]
            bg_size = srs.shape[0]
            go_size = srs.notna().sum()
            de_size = len(gene_set)
            go_de_size = srs.loc[gene_set].notna().sum()

            expected_fraction = go_size / bg_size
            observed_fraction = go_de_size / de_size
            log2_fold_enrichment = np.log2(
                observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
            pval = self._calc_hypergeometric_pval(bg_size=bg_size, go_size=go_size,
                                                  de_size=de_size, go_de_size=go_de_size)
            obs, exp = int(de_size * observed_fraction), de_size * expected_fraction

            enriched_list.append(
                (attribute, de_size, obs, exp, log2_fold_enrichment, pval))

        return self._enrichment_output(enriched_list, fdr, save_csv, fname, return_fig, True, plot_horizontal)

    @staticmethod
    def _calc_hypergeometric_pval(bg_size: int, go_size: int, de_size: int, go_de_size: int):

        """
        Performs a hypergeometric test on the given enrichment set. \
        Given M genes in the background set, n genes in the test set, \
        with N genes from the background set belonging to a specific attribute (or 'success') \
        and X genes from the test set belonging to that attribute. \
        If we were to randomly draw n genes from the background set (without replacement), \
        what is the probability of drawing X or more (in case of enrichment)/X or less (in case of depletion) \
        genes belonging to the given attribute?

        :param bg_size: size of the background set. Usually denoted as 'M'.
        :type bg_size: positive int
        :param go_size: number of features in the background set corresponding to the attribute, \
        or number of successes in the population. Usually denoted as 'n'.
        :type go_size: positive int
        :param de_size: size of the differentially-expressed set, or size of test set. usually denoted as 'N'.
        :type de_size: positive int
        :param go_de_size: or number of successes in the test set. Usually denoted as 'x' or 'k'. s
        :type go_de_size: non-negative int
        :return: p-value of the hypergeometric test.
        :rtype: float between 0 and 1

        """
        if go_de_size / de_size < go_size / bg_size:
            return hypergeom.cdf(go_de_size, bg_size, go_size, de_size)
        return hypergeom.sf(go_de_size - 1, bg_size, go_size, de_size)

    @staticmethod
    def plot_enrichment_results(df: pd.DataFrame, en_score_col: str = 'log2_fold_enrichment', title: str = '',
                                ylabel: str = r"$\log_2$(Fold Enrichment)", plot_horizontal: bool = True,
                                center_bars: bool = True):

        """
        Receives a DataFrame output from an enrichment function and plots it in a bar plot. \
        For the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :param df: a pandas DataFrame created by FeatureSet.enrich_randomization.
        :type df: pd.DataFrame
        :param en_score_col: name of the DataFrame column that contains the enrichment scores.
        :type en_score_col: str (default 'log2_fold_enrichment')
        :param title: plot title.
        :type title: str
        :param ylabel: plot ylabel.
        :type ylabel: str
        :param plot_horizontal:
        :type plot_horizontal: bool (default True)
        :param center_bars: if True, centers the bars around Y=0. Otherwise, ylim is determined by min/max values.
        :type center_bars: bool (default True)
        :return: Figure object containing the bar plot
        :rtype: matplotlib.figure.Figure instance
        """
        plt.style.use('seaborn-white')
        # choose functions and parameters according to the graph's orientation (horizontal vs vertical)
        if plot_horizontal:
            figsize = [5.6, 0.4 * (6.4 + df.shape[0])]
            bar_func = plt.Axes.barh
            line_func = plt.Axes.axvline
            cbar_kwargs = dict(location='bottom')
            tick_func = plt.Axes.set_yticks
            ticklabels_func = plt.Axes.set_yticklabels
            ticklabels_kwargs = dict(fontsize=13, rotation=0)
        else:
            figsize = [0.5 * (6.4 + df.shape[0]), 5.6]
            bar_func = plt.Axes.bar
            line_func = plt.Axes.axhline
            cbar_kwargs = dict(location='left')
            tick_func = plt.Axes.set_xticks
            ticklabels_func = plt.Axes.set_xticklabels
            ticklabels_kwargs = dict(fontsize=13, rotation=45)

        # pull names/scores/pvals out to avoid accidentally changing the results DataFrame in-place
        enrichment_names = df.index.values.tolist()
        enrichment_scores = df[en_score_col].values.tolist()
        enrichment_pvalue = df['padj'].values.tolist()

        # set enrichment scores which are 'inf' or '-inf' to be the second highest/lowest enrichment score in the list
        scores_no_inf = [i for i in enrichment_scores if i != np.inf and i != -np.inf and i < 0]
        if len(scores_no_inf) == 0:
            scores_no_inf.append(-1)
        for i in range(len(enrichment_scores)):
            if enrichment_scores[i] == -np.inf:
                enrichment_scores[i] = min(scores_no_inf)
        max_score = max(np.max(np.abs(enrichment_scores)), 3)

        # get color values for bars
        data_color = [(i / 3) * 127.5 for i in enrichment_scores]
        data_color_norm = [i + 127.5 for i in data_color]
        data_color_norm_256 = [int(i) if i != np.inf and i != -np.inf else np.sign(i) * max(np.abs(scores_no_inf)) for i
                               in data_color_norm]
        my_cmap = plt.cm.get_cmap('coolwarm')
        colors = my_cmap(data_color_norm_256)

        # generate bar plot
        fig, ax = plt.subplots(constrained_layout=True, figsize=figsize)
        bar = bar_func(ax, range(len(enrichment_names)), enrichment_scores, color=colors, edgecolor='black',
                       linewidth=1, zorder=2)
        bar.tick_labels = enrichment_names
        bounds = np.array([np.ceil(-max_score) - 1, (np.floor(max_score) + 1) * 1.002])
        # add black line at y=0 and grey lines at every round positive/negative integer in range
        for ind in range(int(bounds[0]), int(bounds[1]) + 1):
            color = 'black' if ind == 0 else 'grey'
            linewidth = 1 if ind == 0 else 0.5
            linestyle = '-' if ind == 0 else '-.'
            line_func(ax, ind, color=color, linewidth=linewidth, linestyle=linestyle, zorder=0)
        # add colorbar
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(*bounds))
        sm.set_array(np.array([]))
        cbar_label_kwargs = dict(label=ylabel, fontsize=16, labelpad=15)
        cbar = fig.colorbar(sm, ticks=range(int(bounds[0]), int(bounds[1]) + 1), **cbar_kwargs)
        cbar.set_label(**cbar_label_kwargs)
        cbar.ax.tick_params(labelsize=14, pad=6)
        # apply xticks
        tick_func(ax, range(len(enrichment_names)))
        ticklabels_func(ax, enrichment_names, **ticklabels_kwargs)
        # title
        ax.set_title(title, fontsize=18)
        # add significance asterisks
        for col, sig in zip(bar, enrichment_pvalue):
            asterisks, fontweight = FeatureSet._get_pval_asterisk(sig)
            if plot_horizontal:
                x = col._width
                y = col.xy[1] + 0.5 * col._height
                valign = 'center'
                halign = 'left' if np.sign(col._width) == 1 else 'right'
                rotation = 270 if np.sign(col._width) == 1 else 90
            else:
                x = col.xy[0] + 0.5 * col._width
                y = col._height
                valign = 'bottom' if np.sign(col._height) == 1 else 'top'
                halign = 'center'
                rotation = 0

            ax.text(x=x, y=y, s=asterisks, fontname='DejaVu Sans', fontweight=fontweight, rotation=rotation,
                    fontsize=12, horizontalalignment=halign, verticalalignment=valign, zorder=1)
        # despine
        _ = [ax.spines[side].set_visible(False) for side in ['top', 'right']]
        # center bars
        if center_bars:
            if plot_horizontal:
                ax.set_xbound(bounds)
                plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            else:
                ax.set_ybound(bounds)
                plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        plt.show()
        return fig

    @staticmethod
    def _get_pval_asterisk(pval: float):
        fontweight = 'bold'
        if pval < 0.0001:
            asterisks = u'\u2217' * 4
        elif pval < 0.001:
            asterisks = u'\u2217' * 3
        elif pval < 0.01:
            asterisks = u'\u2217' * 2
        elif pval < 0.05:
            asterisks = u'\u2217'
        else:
            asterisks = 'ns'
            fontweight = 'normal'
        return asterisks, fontweight

    def enrich_non_categorical(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                               fdr: float = 0.05, parametric_test: bool = False, biotype: str = 'protein_coding',
                               background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                               attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                               plot_log_scale: bool = True, plot_style: str = 'overlap', n_bins: int = 50,
                               save_csv: bool = False, fname=None, return_fig: bool = False):
        attr_ref_df, gene_set, attributes = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)
        enriched_list = []
        for k, attribute in enumerate(attributes):
            srs = attr_ref_df[attribute]
            if not parametric_test:
                exp, obs = srs.median(), srs[gene_set].median()
                _, pval = sign_test(srs[gene_set].values, exp)
            else:
                exp, obs = srs.mean(), srs[gene_set].mean()
                _, pval = ttest_1samp(srs[gene_set], popmean=exp, nan_policy='propagate')

            enriched_list.append(
                (attribute, len(gene_set), obs, exp, np.nan, pval))

        results_df = self._enrichment_output(enriched_list, fdr, save_csv, fname, return_fig, False)
        results_df.dropna(axis=1, inplace=True)
        for attribute, pval in zip(attributes, results_df['padj']):
            self._enrichment_plot_histogram(attribute, gene_set, self.set_name, attr_ref_df[attribute], plot_log_scale,
                                            plot_style, n_bins, parametric_test, pval)
        return results_df

    @staticmethod
    def _enrichment_plot_histogram(attribute: str, gene_set: set, set_name: str, attr_srs: pd.Series,
                                   plot_log_scale: bool, plot_style: str, n_bins: int, parametric_test: bool,
                                   pval: float):
        assert isinstance(n_bins,
                          int) and n_bins > 0, f"'n_bins' must be a positive integer. Instead got {type(n_bins)}."
        # generate observed and expected Series, either linear or in log10 scale
        exp = attr_srs
        obs = exp.loc[gene_set]
        if plot_log_scale:
            xlabel = r"$\log_{10}$" + f"({attribute})"
            exp = np.log10(exp)
            obs = np.log10(obs)
        else:
            xlabel = f"{attribute}"

        # determine bins according to value range and 'n_bins'
        bins = np.linspace(np.min(exp), np.max(exp), n_bins).squeeze()

        # generate histogram according to plot style
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        kwargs = dict(bins=bins, density=True, alpha=0.5, edgecolor='black', linewidth=1)
        colors = ['C0', 'C1']
        if plot_style.lower() == 'interleaved':
            y, x, _ = ax.hist([exp.values, obs.values], **kwargs, color=colors, label=['Expected', 'Observed'])
            max_y_val = np.max(y)
        elif plot_style.lower() == 'overlap':
            y, _, _ = ax.hist(exp.values, **kwargs, color=colors[0], label='Expected')
            y2, _, _ = ax.hist(obs.values, **kwargs, color=colors[1], label='Observed')
            max_y_val = np.max([np.max(y), np.max(y2)])
        else:
            raise ValueError(f"Invalid value for 'plot_style': '{plot_style}'")

        # set either mean or median as the measure of centrality
        if parametric_test:
            x_exp, x_obs = exp.mean(), obs.mean()
            label_exp, label_obs = 'Expected mean', 'Observed mean'
        else:
            x_exp, x_obs = exp.median(), obs.median()
            label_exp, label_obs = 'Expected median', 'Observed median'

        # add lines for mean/median of observed and expected distributions
        ax.vlines(x_exp, ymin=0, ymax=max_y_val * 1.1, color='blue', linestyle='dashed', linewidth=2,
                  label=label_exp)
        ax.vlines(x_obs, ymin=0, ymax=max_y_val * 1.1, color='red', linestyle='dashed', linewidth=2,
                  label=label_obs)

        # add significance notation
        asterisks, fontweight = FeatureSet._get_pval_asterisk(pval)
        ax.vlines([x_exp, x_obs], ymin=max_y_val * 1.12, ymax=max_y_val * 1.16, color='k', linewidth=1)
        ax.hlines(max_y_val * 1.16, xmin=min(x_exp, x_obs), xmax=max(x_exp, x_obs), color='k', linewidth=1)
        ax.text(np.mean([x_exp, x_obs]), max_y_val * 1.17, asterisks, horizontalalignment='center',
                fontweight=fontweight)

        # legend and titles
        ax.legend()
        obs_name = 'Observed' if set_name == '' else set_name
        ax.set_title(f"Histogram of {attribute} - {obs_name} vs Expected", fontsize=17)
        ax.set_ylabel("Probability density", fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        plt.show()

        return fig

    def biotypes(self, ref: str = 'predefined'):

        """
        Returns a DataFrame of the biotypes in the gene set and their count.

        :type ref: str or pathlib.Path (default 'predefined')
        :param ref: Path of the reference file used to determine biotype. \
        Default is the path predefined in the settings file.

        :Examples:
            >>> from rnalysis import enrichment, filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> en = enrichment.FeatureSet(d)
            >>> en.biotypes(ref='tests/biotype_ref_table_for_tests.csv')
                            gene
            biotype
            protein_coding    26
            pseudogene         1
            unknown            1

        """

        ref = utils.get_biotype_ref_path(ref)
        ref_df = utils.load_csv(ref)
        utils.biotype_table_assertions(ref_df)
        ref_df.columns = ref_df.columns.str.lower()
        not_in_ref = pd.Index(self.gene_set).difference(set(ref_df['gene']))
        if len(not_in_ref) > 0:
            warnings.warn(
                f'{len(not_in_ref)} of the features in the Filter object do not appear in the Biotype Reference Table. ')
            ref_df = ref_df.append(pd.DataFrame({'gene': not_in_ref, 'biotype': 'not_in_biotype_reference'}))
        return ref_df.set_index('gene', drop=False).loc[self.gene_set].groupby('biotype').count()


class RankedSet(FeatureSet):
    """
    Receives a ranked gene set, sorted by any biologically-meaningful value (expression level, fold-change, etc)\
     and preforms various enrichment analyses on them. \
     ALl functions that can be applied to FeatureSet objects can also be applied to RankedSet objects.
    """
    __slots__ = {'ranked_genes': 'a vector of feature names/indices ordered by rank'}

    def __init__(self, ranked_genes: Union[Filter, List[str], Tuple[str], np.ndarray], set_name: str = ''):

        if utils.isinstanceinh(ranked_genes, Filter):
            self.ranked_genes = ranked_genes.df.index.values.astype('str', copy=True)
        elif isinstance(ranked_genes, (list, tuple)):
            self.ranked_genes = np.array(ranked_genes, dtype='str')
        elif isinstance(ranked_genes, np.ndarray):
            self.ranked_genes = ranked_genes.astype('str', copy=True)
        elif isinstance(ranked_genes, set):
            raise TypeError("'ranked_genes' must be an array, list, tuple or Filter object, sorted by rank. "
                            "Python sets are not a valid type ofr 'ranked_genes'.")
        elif isinstance(ranked_genes, dict):
            raise TypeError("Generating a RankedSet from a dictionary is not implemented yet.")
        else:
            raise TypeError(f"'ranked_genes' must be an array, list, tuple or Filter object, sorted by rank. "
                            f"Instead got {type(ranked_genes)}.")

        super().__init__(ranked_genes, set_name)
        assert len(self.ranked_genes) == len(self.gene_set), f"'ranked_genes' must have no repeating elements!"

    def _set_ops(self, others, op: Callable):
        warnings.warn("Warning: when performing set operations with RankedSet objects, "
                      "the return type will always be FeatureSet and not RankedSet.")
        super()._set_ops(others, op)

    def enrich_single_list(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                           attr_ref_path: str = 'predefined', save_csv: bool = False, fname=None,
                           return_fig: bool = False, plot_horizontal: bool = True):
        """
        Calculates enrichment and depletion of the sorted RankedSet for user-defined attributes \
        WITHOUT a background set, using the generalized Minimum Hypergeometric Test (XL-mHG, developed by  \
        `Prof. Zohar Yakhini and colleagues <https://dx.doi.org/10.1371/journal.pcbi.0030039/>`_ \
         and generalized by \
         `Florian Wagner <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143196/>`_). \
        The attributes are drawn from an Attribute Reference Table. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. P-values are calculated using a \
        randomization test with the formula p = (successes + 1)/(repeats + 1). \
        This formula results in a positively-biased estimator of the real p-value \
        (a conservative estimate of p-value). When the number of reps approaches infinity, \
        the formula results in an unbiased estimator of the real p-value. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results_single_list.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_single_list()


        .. figure::  plot_enrichment_results_single_list_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_single_list(plot_horizontal = False)

        """
        attr_ref_df, gene_set, attributes = self._enrichment_setup(biotype='all', background_genes=None,
                                                                   attr_ref_path=attr_ref_path, biotype_ref_path='',
                                                                   attributes=attributes)
        if len(gene_set) == len(self.ranked_genes):
            ranked_genes = self.ranked_genes
        else:
            ranked_genes = np.empty((len(gene_set),), dtype=self.ranked_genes.dtype)
            i = 0
            for elem in self.ranked_genes:
                if elem in gene_set:
                    ranked_genes[i] = elem
                    i += 1

        enriched_list = []
        for k, attribute in enumerate(attributes):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            pval, en_score = self._calc_xlmhg_stats(self._xlmhg_index_vector(ranked_genes, attribute, attr_ref_df),
                                                    len(ranked_genes))
            log2_en_score = np.log2(en_score) if en_score > 0 else -np.inf
            enriched_list.append([attribute, len(ranked_genes), log2_en_score, pval])
            print(f"Finished {k + 1} attributes out of {len(attributes)}")

        return self._xlmhg_output(enriched_list, fdr, save_csv, fname, return_fig, plot_horizontal)

    def _xlmhg_output(self, enriched_list: list, fdr: float, save_csv: bool, fname: str, return_fig: bool,
                      plot_horizontal: bool):
        """
        Formats the enrich list into a results Dataframe, saves the DataFrame to csv if requested, \
        plots the enrichment results, and returns either the Dataframe alone or the Dataframe and the Figure object.
        Called at the end of enrich_single_list().

        """
        en_score_col = 'log2_enrichment_score'
        res_df = pd.DataFrame(enriched_list,
                              columns=['name', 'samples', en_score_col, 'pval'])
        res_df.replace(-np.inf, -np.max(np.abs(res_df['log2_enrichment_score'].values)))
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.set_index('name', inplace=True)

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        fig = self.plot_enrichment_results(res_df, en_score_col=en_score_col,
                                           title=f"Single-list enrichment for {self.set_name}",
                                           ylabel=r"$\log_2$(XL-mHG enrichment score)", plot_horizontal=plot_horizontal)
        if return_fig:
            return res_df, fig
        return res_df

    @staticmethod
    def _calc_xlmhg_stats(index_vec: np.ndarray, ranked_genes_len: int):
        index_vec = np.uint16(index_vec)
        rev_index_vec = np.uint16([ranked_genes_len - 1 - index_vec[i - 1] for i in range(len(index_vec), 0, -1)])
        # X = the minimal amount of 'positive' elements above the hypergeometric cutoffs out of all of the positive
        # elements in the ranked set. Determined to be the minimum between x_min and ceil(x_frac * k),
        # where 'k' is the number of 'positive' elements in the ranked set.
        x_frac = 0.25
        x_min = 5
        # L = the lowest possible cutoff (n) to be tested out of the entire list.
        # Determined to be floor(l_frac * N), where 'N' is total number of elements in the ranked set (ranked_genes_len).
        l_frac = 0.25
        # pre-allocate empty array to speed up computation
        table = np.empty((len(index_vec) + 1, ranked_genes_len - len(index_vec) + 1), dtype=np.longdouble)
        res_obj_fwd = xlmhg_test(N=ranked_genes_len, indices=index_vec, L=int(np.floor(l_frac * ranked_genes_len)),
                                 X=min(x_min, int(np.ceil(x_frac * len(index_vec)))), table=table)
        res_obj_rev = xlmhg_test(N=ranked_genes_len, indices=rev_index_vec, L=int(np.floor(l_frac * ranked_genes_len)),
                                 X=min(x_min, int(np.ceil(x_frac * len(index_vec)))), table=table)

        if res_obj_fwd.pval <= res_obj_rev.pval:
            pval, en_score = res_obj_fwd.pval, res_obj_fwd.escore
        else:
            pval, en_score = res_obj_rev.pval, 1 / res_obj_rev.escore
        pval = pval if not np.isnan(pval) else 1
        en_score = en_score if not np.isnan(en_score) else 1
        return pval, en_score

    @staticmethod
    def _xlmhg_index_vector(ranked_genes, attribute, attr_ref_df):
        ranked_srs = attr_ref_df.loc[ranked_genes, attribute]
        assert ranked_srs.shape[0] == len(ranked_genes)
        return np.uint16(np.nonzero(ranked_srs.notna().values)[0])


def _fetch_sets(objs: dict, ref: str = 'predefined'):
    """
    Receives the 'objs' input from enrichment.upset_plot() and enrichment.venn_diagram(), and turns the values in it \
    into python sets.

    :param objs: the 'objs' input given to the function enrichment.upset_plot() or enrichment.venn_diagram().
    :type objs: a dictionary, where the keys are names of sets, and the values are either\
     python sets, FeatureSets or names of columns in the Attribute Reference Table.
    :param ref: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type ref: str or pathlib.Path (default 'predefined')
    :return: a dictionary, where the keys are names of sets and the values are python sets of feature indices.
    """
    assert isinstance(objs, dict), f"objs must be a dictionary. Instaed got {type(objs)}"
    for obj in objs:
        if isinstance(objs[obj], set):
            pass
        elif utils.isinstanceinh(objs[obj], Filter):
            objs[obj] = objs[obj].index_set
        elif isinstance(objs[obj], FeatureSet):
            objs[obj] = objs[obj].gene_set
        elif isinstance(objs[obj], str):
            if 'attr_table' not in locals():
                pth = utils.get_attr_ref_path(ref)
                attr_table = utils.load_csv(pth, 0)
            attr = objs[obj]
            myset = set(attr_table[attr].loc[attr_table[attr].notna()].index)
            objs[obj] = myset
        else:
            raise TypeError
    return objs


def upset_plot(objs: Dict[str, Union[str, FeatureSet, Set[str]]], title: str = '', ref: str = 'predefined'):
    """
    Generate an UpSet plot of 2 or more sets, FeatureSets or attributes from the Attribute Reference Table.


    :param objs: the FeatureSets, python sets or user-defined attributes to plot.
    :type objs: a dictionary with 2 or more entries, where the keys are the names of the sets, and the values are either \
    a FeatureSet, a python set of feature indices, or a name of a column in the Attribute Reference Table. \
    For example: \
    {'first set':{'gene1','gene2','gene3'}, 'second set':'name_of_attribute_from_reference_table'}
    :param title: determines the title of the plot.
    :type title: str
    :param ref: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type ref: str or pathlib.Path (default 'predefined')
    :returns: a dictionary of matplotlib axes, where the keys are 'matrix', 'intersections', 'totals', 'shading'.


        .. figure::  upsetplot.png
           :align:   center
           :scale: 70 %

           Example plot of upset_plot()
    """

    upset_df = _generate_upset_srs(_fetch_sets(objs=objs, ref=ref))
    upsetplot = upset.plot(upset_df)
    plt.title(title)
    return upsetplot


def venn_diagram(objs: Dict[str, Union[str, FeatureSet, Set[str]]], title: str = 'default', ref: str = 'predefined',
                 set_colors: tuple = ('r', 'g', 'b'),
                 alpha: float = 0.4, weighted: bool = True, lines: bool = True, linecolor: str = 'black',
                 linestyle='solid', linewidth=2.0,
                 normalize_to: float = 1.0):
    """
    Generate a Venn diagram of 2 to 3 sets, FeatureSets or attributes from the Attribute Reference Table.

    :param objs: the FeatureSets, python sets or user-defined attributes to plot.
    :type objs: a dictionary with 2-3 entries, where the keys are the names of the sets, and the values are either \
    a FeatureSet, a python set of feature indices, or a name of a column in the Attribute Reference Table. \
    For example: \
    {'first set':{'gene1','gene2','gene3'}, 'second set':'name_of_attribute_from_reference_table'}
    :type title: str
    :param set_colors: determines the colors of the circles in the diagram.
    :param ref: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type ref: str or pathlib.Path (default 'predefined')
    :param title: determines the title of the plot.
    :type set_colors: tuple of matplotlib-format colors, the same size as 'objs'
    :param alpha: determines the opacity of the circles.
    :type alpha: a float between 0 and 1
    :param weighted: if True, the plot will be area-weighted.
    :type weighted: bool (default True)
    :param lines: if True, adds an outline to the circles.
    :type lines: bool (default True)
    :param linecolor: Determines the color of the circles' outline.
    :type linecolor: matplotlib-format color (default 'black')
    :param linestyle: the style of the circles' outline.
    :type linestyle: 'solid' or 'dashed' (default 'solid')
    :param linewidth: the widdth of the circles' outlines.
    :type linewidth: float (default 2.0)
    :param normalize_to:
    :type normalize_to: float (default 1.0)
    :return: a tuple of a VennDiagram object; and a list of 2-3 Circle patches.


        .. figure::  venn.png
           :align:   center
           :scale: 70 %

           Example plot of venn_diagram()
    """
    if len(objs) > 3 or len(objs) < 2:
        raise ValueError(f'Venn can only accept between 2 and 3 sets. Instead got {len(objs)}')
    assert isinstance(title, str), f'Title must be a string. Instead got {type(title)}'
    objs = _fetch_sets(objs=objs, ref=ref)
    if len(objs) == 2:
        func = vn.venn2 if weighted else vn.venn2_unweighted
        func_circles = vn.venn2_circles
        set_colors = set_colors[0:2]
    else:
        func = vn.venn3 if weighted else vn.venn3_unweighted
        func_circles = vn.venn3_circles
    _ = plt.figure()
    plot_obj = func(tuple(objs.values()), tuple(objs.keys()), set_colors=set_colors, alpha=alpha,
                    normalize_to=normalize_to)
    if lines and weighted:
        circle_obj = func_circles(tuple(objs.values()), color=linecolor, linestyle=linestyle, linewidth=linewidth,
                                  normalize_to=normalize_to)
    elif lines and not weighted:
        warnings.warn('Cannot draw lines on an unweighted venn diagram. ')
        circle_obj = None
    else:
        circle_obj = None
    if title == 'default':
        title = 'Venn diagram of ' + ''.join([name + ' ' for name in objs.keys()])
    plt.title(title)
    return plot_obj, circle_obj


def _generate_upset_srs(objs: dict):
    """
    Receives a dictionary of sets from enrichment._fetch_sets(), \
    and reformats it as a pandas Series to be used by the python package 'upsetplot'.

    :param objs: the output of the enrichment._fetch_sets() function.
    :type objs: dict of sets
    :return: a pandas Series in the format requested by the 'upsetplot' package.

    """
    names = list(objs.keys())
    multi_ind = pd.MultiIndex.from_product([[True, False] for _ in range(len(names))], names=names)[:-1]
    srs = pd.Series(index=multi_ind)
    for ind in multi_ind:
        sets = list(compress(names, ind))
        group_size = len(set.intersection(*[objs[s] for s in sets]))
        srs.loc[ind] = group_size
    return srs

# TODO: other types of plots
# TODO: heat map plot of multiple DESEQ files
