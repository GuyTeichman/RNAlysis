import collections
import itertools
import logging
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Iterable, List, Tuple, Union, Collection, Set, Dict, Literal

import joblib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as multitest
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import hypergeom, ttest_1samp, fisher_exact
from statsmodels.stats.descriptivestats import sign_test
from tqdm.auto import tqdm

from rnalysis.utils import ontology, io, parsing, settings, validation, generic, param_typing
from rnalysis.utils.param_typing import PARALLEL_BACKENDS, GRAPHVIZ_FORMATS

try:
    import xlmhglite

    logging.getLogger('xlmhg').setLevel(50)  # suppress warnings from xlmhg module
    HAS_XLMHG = True

except ImportError:  # pragma: no cover
    HAS_XLMHG = False


    class xlmhglite:  # pragma: no cover
        class mHGResult:
            pass

        def get_xlmhg_test_result(self):
            pass


class Size:
    def __init__(self, size: int):
        self.size = size


def dot_scaling_func(values: Union[float, np.ndarray, pd.Series]):
    return 8 * np.array(values) ** 0.75


class sizeHandler:
    def __init__(self, size: int):
        self.size = size

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        radius = np.sqrt((dot_scaling_func(self.size)) / np.pi)
        patch = mpatches.Circle((x0, y0), radius, facecolor='black',
                                edgecolor='black', hatch='xx', lw=3,
                                transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        return patch


class EnrichmentRunner:
    ENRICHMENT_SCORE_YLABEL = r"$\log_2$(Fold Enrichment)"
    SINGLE_SET_ENRICHMENT_SCORE_YLABEL = r"$\log_2$(XL-mHG enrichment score)"

    __slots__ = {'results': 'DataFrame containing enrichment analysis results',
                 'annotation_df': 'DataFrame containing all annotation data per gene',
                 'gene_set': 'the set of genes/genomic features whose enrichment to calculate',
                 'attributes': 'the list of attributes/terms to calculate enrichment for',
                 'alpha': 'the statistical signifiacnce threshold',
                 'attr_ref_path': 'path of the Attribute Reference Table to load, if such table exists',
                 'return_nonsignificant': 'indicates whether to return results which were not found to be '
                                          'statistically significant after enrichment analysis',
                 'save_csv': 'indicates whether the results should be saved to a csv file',
                 'fname': 'name of the file to save results to',
                 'return_fig': 'indicates whether to return a matplotlib Figure after plotting the results',
                 'plot_horizontal': 'indicates whether to plot the results horizontally or vertically',
                 'set_name': 'name of the set of genes/genomic features whose enrichment is calculated',
                 'parallel_backend': 'indicates whether to use parallel processing, and with which backend',
                 'exclude_unannotated_genes': 'indicates whether to discard unannotated genes from the background and enrichment sets',
                 'enrichment_func': 'the function to be used to calculate enrichment p-values',
                 'pvalue_kwargs': 'key-worded arguments for the function which calculates enrichment p-values',
                 'single_set': 'indicates whether enrichment is calculated on a single set of genes '
                               '(without background) or on a set of target genes and a set of background genes',
                 'biotypes': 'the requested biotype of the background gene set',
                 'background_set': 'the background set of genes for enrichment analysis',
                 'biotype_ref_path': 'path of the Biotype Reference Table to load, if such table exists',
                 'random_seed': 'random seed to be used when non-deterministic functions are used',
                 'en_score_col': 'name of the enrichment score column in the results DataFrame',
                 'ranked_genes': 'the set of genes/genomic features whose enrichment to calculate, '
                                 'pre-sorted and ranked by the user',
                 'plot_style': 'plot style',
                 'show_expected': 'show observed/expected values on plot'}
    printout_params = "appear in the Attribute Reference Table"

    def __init__(self, genes: Union[set, np.ndarray], attributes: Union[Iterable, str, int],
                 alpha: param_typing.Fraction, attr_ref_path: str, return_nonsignificant: bool, save_csv: bool,
                 fname: str, return_fig: bool, plot_horizontal: bool, set_name: str,
                 parallel_backend: Literal[PARALLEL_BACKENDS], enrichment_func_name: str, biotypes=None,
                 background_set: set = None, biotype_ref_path: str = None, exclude_unannotated_genes: bool = True,
                 single_set: bool = False, random_seed: int = None, plot_style: Literal['bar', 'lollipop'] = 'bar',
                 show_expected: bool = False, **pvalue_kwargs):
        self.results: pd.DataFrame = pd.DataFrame()
        self.annotation_df: pd.DataFrame = pd.DataFrame()
        self.gene_set = parsing.data_to_set(genes)
        self.attributes = attributes
        self.alpha = alpha
        self.attr_ref_path = settings.get_attr_ref_path(attr_ref_path)
        self.return_nonsignificant = return_nonsignificant
        self.save_csv = save_csv

        if self.save_csv:
            if fname is None:
                self.fname = input("Please insert the full name and path to save the file to")
            else:
                assert isinstance(fname, (str, Path))
                self.fname = str(fname)
        self.return_fig = return_fig
        self.plot_horizontal = plot_horizontal
        self.plot_style = plot_style
        self.show_expected = show_expected
        self.set_name = set_name
        self.parallel_backend = parallel_backend
        self.enrichment_func = self._get_enrichment_func(enrichment_func_name)
        if not self.enrichment_func:
            return
        self.pvalue_kwargs = pvalue_kwargs
        self.exclude_unannotated_genes = exclude_unannotated_genes
        self.single_set = single_set
        if self.single_set:
            assert biotypes is None, "Enrichment in single_set mode does not accept a 'biotypes' argument."
            assert background_set is None, "Enrichment in single_set mode does not accept a 'background_set' argument."
            assert biotype_ref_path is None, \
                "Enrichment in single_set mode does not accept a 'biotype_ref_path' argument."
            assert random_seed is None, "Enrichment in single_set mode does not accept a 'random_seed' argument."
            assert isinstance(genes, np.ndarray), f"Invalid type for argument 'genes' in single_set mode: " \
                                                  f"expected np.ndarray, instead got '{type(genes)}'."

            legal_enrichment_funcs = {'xlmhg'}
            assert enrichment_func_name.lower() in legal_enrichment_funcs, \
                f"Invalid enrichment_func_name for single_set mode: '{enrichment_func_name}'."

            self.biotypes, self.background_set, self.biotype_ref_path, self.random_seed = None, None, None, None
            self.en_score_col = 'log2_enrichment_score'
            self.ranked_genes = genes
            self._update_ranked_genes()
        else:
            assert (isinstance(biotypes, (str, list, set, tuple)))
            self.random_seed = random_seed
            self.biotypes = biotypes
            self.background_set = background_set
            self.biotype_ref_path = settings.get_biotype_ref_path(biotype_ref_path) if biotypes != 'all' else None
            self.en_score_col = 'log2_fold_enrichment'
            self.ranked_genes = None

    @classmethod
    def from_results_df(cls, results: pd.DataFrame, alpha: float, plot_horizontal: bool, set_name: str):
        runner = cls.__new__(cls)
        runner.results = results
        runner.alpha = alpha
        runner.plot_horizontal = plot_horizontal
        runner.set_name = set_name
        return runner

    def run(self, plot: bool = True) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        if not self.enrichment_func:
            return pd.DataFrame()
        self.fetch_annotations()
        self.fetch_attributes()
        if not self.single_set:
            self.get_background_set()
        self.update_gene_set()
        if len(self.gene_set) == 0:
            warnings.warn('After removing unannotated genes and/or genes not in the background set, '
                          'the enrichment set is empty. '
                          'Therefore, RNAlysis will not proceed with enrichment analysis. ')
            return pd.DataFrame()
        self.filter_annotations()
        unformatted_results = self.calculate_enrichment()
        self.format_results(unformatted_results)
        if self.save_csv:
            self.results_to_csv()
        if plot:
            fig = self.plot_results()

            if self.return_fig:
                return self.results, fig
        return self.results

    def _update_ranked_genes(self):
        if len(self.ranked_genes) == len(self.gene_set):
            return
        else:
            ranked_genes = np.empty((len(self.gene_set),), dtype=self.ranked_genes.dtype)
            processed_genes = set()
            i = 0
            for elem in self.ranked_genes:
                if elem in self.gene_set and elem not in processed_genes:
                    ranked_genes[i] = elem
                    processed_genes.add(elem)
                    i += 1

                    if len(processed_genes) == len(self.gene_set):
                        break
            self.ranked_genes = ranked_genes

    def results_to_csv(self):
        io.save_table(self.results, filename=self.fname if self.fname.endswith('.csv') else self.fname + '.csv')

    def get_background_set(self):
        if self.background_set is None:
            self._get_background_set_from_biotype()
        else:
            self._get_background_set_from_set()

        orig_bg_size = len(self.background_set)
        if self.exclude_unannotated_genes:
            self.background_set = self.background_set.intersection(parsing.data_to_set(self.annotation_df.index))
        else:
            ann_missing = self.background_set.difference(parsing.data_to_set(self.annotation_df.index))
            ann_missing_df = pd.DataFrame(index=parsing.data_to_list(ann_missing), columns=self.annotation_df.columns)
            self.annotation_df = pd.concat([self.annotation_df, ann_missing_df])

        if orig_bg_size - len(self.background_set) > 0:
            warning = f"{orig_bg_size - len(self.background_set)} genes out of the requested {orig_bg_size} " \
                      f"background genes do not {self.printout_params}"
            if self.exclude_unannotated_genes:
                warning += f", and are therefore ignored. " \
                           f"\nThis leaves a total of {len(self.background_set)} background genes."
            else:
                warning += ". \nThese genes are not discarded from the background set because you set the paremeter " \
                           "'exclude_unannotated_genes' to False. Note that this is not a recommended practice. "
            warnings.warn(warning)

        print(f"{len(self.background_set)} background genes are used.")

    def _get_enrichment_func(self, pval_func_name: str):
        assert isinstance(pval_func_name, str), f"Invalid type for 'pval_func_name': {type(pval_func_name)}."
        pval_func_name = pval_func_name.lower()
        if pval_func_name == 'fisher':
            return self._fisher_enrichment
        if pval_func_name == 'randomization':
            return self._randomization_enrichment
        elif pval_func_name == 'hypergeometric':
            return self._hypergeometric_enrichment
        elif pval_func_name == 'xlmhg':
            if HAS_XLMHG:
                return self._xlmhg_enrichment
            else:
                if not HAS_XLMHG:
                    warnings.warn("Package 'xlmhg' is not installed. \n"
                                  "If you want to run single-set enrichment analysis, "
                                  "please install package 'xlmhg' and try again. ")

                return False
        else:
            raise ValueError(f"Unknown enrichment function '{pval_func_name}'.")

    def _get_hypergeometric_parameters(self, attribute: str) -> Tuple[int, int, int, int]:
        bg_size = self.annotation_df.shape[0]
        de_size = len(self.gene_set)
        go_size = self.annotation_df[attribute].notna().sum()
        go_de_size = self.annotation_df.loc[list(self.gene_set), attribute].notna().sum()
        return bg_size, de_size, go_size, go_de_size

    def _get_xlmhg_parameters(self, index_vec):
        n = len(self.ranked_genes)
        # X = the minimal amount of 'positive' elements above the hypergeometric cutoffs out of all the positive
        # elements in the ranked set.
        X = 10 if 'x' not in self.pvalue_kwargs else self.pvalue_kwargs['x']
        X = min(X, n)
        # L = the lowest possible cutoff (n) to be tested out of the entire list.
        # Determined to be floor(l_fraction * N), where 'N' is total number of elements in the ranked set (n).
        if 'l_fraction' in self.pvalue_kwargs and self.pvalue_kwargs['l_fraction'] is not None:
            l_fraction = self.pvalue_kwargs['l_fraction']
        else:
            l_fraction = 0.1
        L = int(np.floor(l_fraction * n))
        L = min(L, n)
        # pre-allocate empty array to speed up computation
        table = np.empty((len(index_vec) + 1, n - len(index_vec) + 1), dtype=np.longdouble)
        return n, X, L, table

    def _xlmhg_enrichment(self, attribute: str, **_) -> list:
        index_vec, rev_index_vec = self._generate_xlmhg_index_vectors(attribute)
        n, X, L, table = self._get_xlmhg_parameters(index_vec)

        res_obj_fwd = xlmhglite.get_xlmhg_test_result(N=n, indices=index_vec, X=X, L=L, table=table)
        res_obj_rev = xlmhglite.get_xlmhg_test_result(N=n, indices=rev_index_vec, X=X, L=L, table=table)

        obs, exp, en_score, pval = self._extract_xlmhg_results(res_obj_fwd, res_obj_rev)
        return [attribute, n, obs, exp, en_score, pval]

    @staticmethod
    def _extract_xlmhg_results(result_obj_fwd: xlmhglite.mHGResult, result_obj_rev: xlmhglite.mHGResult
                               ) -> Tuple[int, float, float, float]:
        if result_obj_fwd.pval <= result_obj_rev.pval:
            obs = result_obj_fwd.k
            exp = result_obj_fwd.K * (result_obj_fwd.cutoff / float(result_obj_fwd.N))
            pval = result_obj_fwd.pval
            en_score = result_obj_fwd.escore if not np.isnan(result_obj_fwd.escore) else 1
        else:
            obs = result_obj_rev.k
            exp = result_obj_rev.K * (result_obj_rev.cutoff / float(result_obj_rev.N))
            pval = result_obj_rev.pval
            en_score = 1 / result_obj_rev.escore if not np.isnan(result_obj_rev.escore) else 1
        pval = pval if not np.isnan(pval) else 1
        en_score = en_score if not np.isnan(en_score) else 1
        log2_en_score = np.log2(en_score) if en_score > 0 else -np.inf

        return obs, exp, log2_en_score, pval

    def _generate_xlmhg_index_vectors(self, attribute) -> Tuple[np.ndarray, np.ndarray]:
        n = len(self.ranked_genes)
        ranked_srs = self.annotation_df.loc[self.ranked_genes, attribute]
        assert ranked_srs.shape[0] == n
        index_vec = np.uint16(np.nonzero(ranked_srs.notna().values)[0])
        rev_index_vec = np.uint16([n - 1 - index_vec[i - 1] for i in range(len(index_vec), 0, -1)])
        return index_vec, rev_index_vec

    def _fisher_enrichment(self, attribute: str) -> list:
        bg_size, de_size, go_size, go_de_size = self._get_hypergeometric_parameters(attribute)

        expected_fraction = go_size / bg_size
        observed_fraction = go_de_size / de_size
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_fisher_pval(bg_size=bg_size, de_size=de_size, go_size=go_size, go_de_size=go_de_size)
        obs, exp = int(de_size * observed_fraction), de_size * expected_fraction
        return [attribute, de_size, obs, exp, log2_fold_enrichment, pval]

    def _randomization_enrichment(self, attribute: str, reps: int) -> list:
        bg_array = self.annotation_df[attribute].notna().values
        obs_array = self.annotation_df.loc[list(self.gene_set), attribute].notna().values
        n = len(self.gene_set)
        expected_fraction = np.sum(bg_array) / bg_array.shape[0]
        observed_fraction = np.sum(obs_array) / n
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_randomization_pval(n, log2_fold_enrichment, bg_array, reps, observed_fraction)

        return [attribute, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval]

    def _hypergeometric_enrichment(self, attribute: str):
        bg_size, de_size, go_size, go_de_size = self._get_hypergeometric_parameters(attribute)

        expected_fraction = go_size / bg_size
        observed_fraction = go_de_size / de_size
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_hypergeometric_pval(bg_size=bg_size, de_size=de_size, go_size=go_size, go_de_size=go_de_size)
        obs, exp = int(de_size * observed_fraction), de_size * expected_fraction

        return [attribute, de_size, obs, exp, log2_fold_enrichment, pval]

    @staticmethod
    @generic.numba.jit(nopython=True)
    def _calc_randomization_pval(n: int, log2fc: float, bg_array: np.ndarray, reps: int, obs_frac: float
                                 ) -> float:  # pragma: no cover
        ind_range = np.arange(bg_array.shape[0])
        success = 0
        if log2fc >= 0:
            for _ in range(reps):
                success += np.sum(bg_array[np.random.choice(ind_range, n, replace=False)]) / n >= obs_frac

        else:
            for _ in range(reps):
                success += np.sum(bg_array[np.random.choice(ind_range, n, replace=False)]) / n <= obs_frac
        pval = (success + 1) / (reps + 1)
        return pval

    @staticmethod
    @lru_cache(maxsize=256, typed=False)
    def _calc_fisher_pval(bg_size: int, de_size: int, go_size: int, go_de_size: int) -> float:
        contingency_table = [[go_de_size, go_size - go_de_size],
                             [de_size - go_de_size, bg_size - go_size - de_size + go_de_size]]
        _, pval = fisher_exact(contingency_table)
        return pval

    @staticmethod
    def _calc_hypergeometric_pval(bg_size: int, de_size: int, go_size: int, go_de_size: int) -> float:
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
        :param go_de_size: or number of successes in the differentially-expressed set. Usually denoted as 'x' or 'k'.
        :type go_de_size: non-negative int
        :return: p-value of the hypergeometric test.
        :rtype: float between 0 and 1

        """
        try:
            if go_de_size / de_size < go_size / bg_size:
                return hypergeom.cdf(go_de_size, bg_size, go_size, de_size)
            return hypergeom.sf(go_de_size - 1, bg_size, go_size, de_size)
        except ZeroDivisionError:
            return hypergeom.cdf(go_de_size, bg_size, go_size, de_size)

    def _get_background_set_from_biotype(self):
        if self.biotypes == 'all':
            self.background_set = parsing.data_to_set(self.annotation_df.index)
        else:
            biotype_ref_df = io.load_table(self.biotype_ref_path)
            validation.validate_biotype_table(biotype_ref_df)
            biotype_ref_df.set_index('gene', inplace=True)
            self.biotypes = parsing.data_to_list(self.biotypes)
            mask = pd.Series(np.zeros_like(biotype_ref_df['biotype'].values, dtype=bool),
                             biotype_ref_df['biotype'].index, name='biotype')
            for biotype in self.biotypes:
                mask = mask | (biotype_ref_df['biotype'] == biotype)
            self.background_set = parsing.data_to_set(biotype_ref_df[mask].index)

    def _get_background_set_from_set(self):
        self.background_set = parsing.data_to_set(self.background_set)
        if self.biotypes != 'all':
            warnings.warn("both 'biotype' and 'background_genes' were specified. Therefore 'biotype' is ignored.")

    def update_gene_set(self):
        if self.single_set:
            updated_gene_set = self.gene_set.intersection(parsing.data_to_set(self.annotation_df.index))
            not_annotated = len(self.gene_set) - len(updated_gene_set)
            self.gene_set = updated_gene_set
            self._update_ranked_genes()
            if not_annotated > 0:
                warnings.warn(f"{not_annotated} genes in the enrichment set do are not annotated. \n"
                              f"Enrichment will be computed on the remaining {len(self.gene_set)} genes.")

        else:
            updated_gene_set = self.gene_set.intersection(self.background_set)
            not_in_bg = len(self.gene_set) - len(updated_gene_set)
            self.gene_set = updated_gene_set
            if not_in_bg > 0:
                warning = f"{not_in_bg} genes in the enrichment set do not appear in the background gene set"
                if self.exclude_unannotated_genes:
                    warning += f" and/or do not {self.printout_params}"
                warning += f". \nEnrichment will be computed on the remaining {len(self.gene_set)} genes."
                warnings.warn(warning)

    def fetch_annotations(self):
        self.annotation_df = io.load_table(self.attr_ref_path)
        validation.validate_attr_table(self.annotation_df)
        self.annotation_df.set_index('gene', inplace=True)

    def fetch_attributes(self):
        all_attrs = self.annotation_df.columns
        if self.attributes == 'all':
            self.attributes = parsing.data_to_list(all_attrs)
        else:
            if self.attributes is None:
                attribute_list = parsing.from_string(
                    "Please insert attributes separated by newline "
                    "(for example: \n'epigenetic_related_genes\nnrde-3 targets\nALG-3/4 class small RNAs')")
            else:
                attribute_list = parsing.data_to_list(self.attributes)
            self._validate_attributes(attribute_list, all_attrs)
            self.attributes = [all_attrs[ind] for ind in attribute_list] if \
                validation.isinstanceiter(attribute_list, int) else attribute_list

    @staticmethod
    def _validate_attributes(attribute_list: list, all_attrs: Collection):
        assert isinstance(attribute_list, list), f"Invalid type for 'attribute_list': {type(attribute_list)}."
        all_attrs = parsing.data_to_set(all_attrs)
        for attr in attribute_list:
            if type(attr) in {int}:
                assert attr >= 0, f"Error in attribute number {attr}: index must be non-negative!"
                assert attr < len(all_attrs), f"Attribute index {attr} out of range."
            else:
                assert isinstance(attr, str), f"Invalid type of attribute {attr}: {type(attr)}"
                assert attr in all_attrs, f"Attribute {attr} does not appear in the Attribute Refernce Table."

    def filter_annotations(self):
        if self.single_set:
            self.annotation_df = self.annotation_df.loc[:, self.attributes].sort_index()
        else:
            self.annotation_df = self.annotation_df.loc[list(self.background_set), self.attributes].sort_index()

    def calculate_enrichment(self) -> list:
        self.set_random_seed()
        if self.parallel_backend != 'sequential' and len(self.attributes) > 5:
            return self._calculate_enrichment_parallel()
        return self._calculate_enrichment_serial()

    def set_random_seed(self):
        if self.random_seed is not None:
            assert isinstance(self.random_seed, int) and self.random_seed >= 0, \
                f"random_seed must be a non-negative integer. Value '{self.random_seed}' is invalid."
            np.random.seed(self.random_seed)

    def _calculate_enrichment_serial(self) -> list:
        result = []
        for attribute in tqdm(self.attributes, desc="Calculating enrichment", unit='attributes'):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            result.append(self.enrichment_func(attribute, **self.pvalue_kwargs))
            # print(f"Finished {n_attrs + 1} attributes out of {len(self.attributes)}", end='\r')
        return result

    def _calculate_enrichment_parallel(self) -> list:
        result = generic.ProgressParallel(desc="Calculating enrichment", unit='attribute', n_jobs=-1,
                                          backend=self.parallel_backend)(
            joblib.delayed(self.enrichment_func)(attribute, **self.pvalue_kwargs) for attribute in self.attributes)
        return result

    def format_results(self, unformatted_results_list: list):
        columns = ['name', 'samples', 'obs', 'exp', self.en_score_col, 'pval']
        self.results = pd.DataFrame(unformatted_results_list, columns=columns).set_index('name')
        self._correct_multiple_comparisons()

        # filter non-significant results
        if not self.return_nonsignificant:
            self.results = self.results[self.results['significant']]

    def _correct_multiple_comparisons(self):
        significant, padj = multitest.fdrcorrection(self.results.loc[self.results['pval'].notna(), 'pval'].values,
                                                    alpha=self.alpha)
        self.results.loc[self.results['pval'].notna(), 'padj'] = padj
        self.results['significant'] = False  # set default value as False
        self.results.loc[self.results['padj'].notna(), 'significant'] = significant

    def plot_results(self) -> plt.Figure:
        if self.single_set:
            return self.enrichment_bar_plot(ylabel=self.SINGLE_SET_ENRICHMENT_SCORE_YLABEL,
                                            title=f"Single-list enrichment for gene set '{self.set_name}'")
        return self.enrichment_bar_plot(ylabel=self.ENRICHMENT_SCORE_YLABEL,
                                        title=f"Enrichment for gene set '{self.set_name}'")

    def enrichment_bar_plot(self, n_bars: Union[param_typing.PositiveInt, Literal['all']] = 'all',
                            center_bars: bool = True,
                            ylabel: str = r"$\log_2$(Fold Enrichment)", title: str = 'Enrichment results',
                            ylim: Union[float, Literal['auto']] = 'auto', title_fontsize: float = 18,
                            label_fontsize: float = 13, ylabel_fontsize: float = 16) -> plt.Figure:

        """
        Receives a DataFrame output from an enrichment function and plots it in a bar plot. \
        For the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :param n_bars: number of bars to display in the bar plot. If n_bars='all', \
         all the results will be displayed on the graph. Otherwise, only the top n results will be displayed on the graph.
        :type n_bars: int > 1 or 'all' (default='all')
        :param title: plot title.
        :type title: str
        :param ylabel: plot y-axis label.
        :type ylabel: str (default=r"$\log_2$(Fold Enrichment)")
        :param center_bars: if True, center the bars around Y=0. Otherwise, ylim is determined by min/max values.
        :type center_bars: bool (default=True)
        :param ylim: set the Y-axis limits. If `ylim`='auto', determines the axis limits automatically based on the data. \
        If `ylim` is a number, set the Y-axis limits to [-ylim, ylim].
        :type ylim: float or 'auto' (default='auto')
        :return: Figure object containing the bar plot
        :rtype: matplotlib.figure.Figure instance
        """
        assert self.plot_style in ['bar', 'lollipop'], \
            f"'plot_style' must be 'bar' or 'lollipop', instaed got '{self.plot_style}'."

        # determine number of entries/bars to plot
        if n_bars != 'all':
            assert isinstance(n_bars, int) and n_bars > 0, f"Invalid value for 'n_bars': {n_bars}."
            results = self.results.iloc[:n_bars]
        else:
            results = self.results
        # pull names/scores/pvals out to avoid accidentally changing the results DataFrame in-place
        if 'name' in results.columns:
            enrichment_names = results['name'].values.tolist()
        else:
            enrichment_names = results.index.values.tolist()
        enrichment_scores = results[self.en_score_col].values.tolist()
        enrichment_pvalue = results['padj'].values.tolist()
        enrichment_obs = results['obs'].values.tolist()
        enrichment_exp = results['exp'].values.tolist()

        # choose functions and parameters according to the graph's orientation (horizontal vs vertical)
        if self.plot_horizontal:
            figsize = [10.5, 0.4 * (4.8 + self.results.shape[0])]
            bar_func = plt.Axes.barh if self.plot_style == 'bar' else plt.Axes.hlines
            line_func = plt.Axes.axvline

            cbar_location = 'bottom'
            cbar_orientation = 'horizontal'
            tick_func = plt.Axes.set_yticks
            ticklabels_func = plt.Axes.set_yticklabels
            ticklabels_kwargs = dict(fontsize=label_fontsize, rotation=0)

            for lst in (enrichment_names, enrichment_scores, enrichment_pvalue, enrichment_obs, enrichment_exp):
                lst.reverse()
        else:
            figsize = [0.5 * (4.8 + self.results.shape[0]), 4.2]
            bar_func = plt.Axes.bar if self.plot_style == 'bar' else plt.Axes.vlines
            line_func = plt.Axes.axhline
            cbar_location = 'left'
            cbar_orientation = 'vertical'
            tick_func = plt.Axes.set_xticks
            ticklabels_func = plt.Axes.set_xticklabels
            ticklabels_kwargs = dict(fontsize=label_fontsize, rotation=45)

        # set enrichment scores which are 'inf' or '-inf' to be the second highest/lowest enrichment score in the list
        scores_no_inf = [abs(score) for score in enrichment_scores if score != np.inf and score != -np.inf]
        if len(scores_no_inf) == 0:
            scores_no_inf.append(3)

        for i in range(len(enrichment_scores)):
            if enrichment_scores[i] == -np.inf:
                enrichment_scores[i] = -max(scores_no_inf)
            elif enrichment_scores[i] == np.inf:
                enrichment_scores[i] = max(scores_no_inf)

        if ylim == 'auto':
            if len(scores_no_inf) > 0:
                max_score = max(max(scores_no_inf), 3)
            else:
                max_score = 3
        else:
            max_score = ylim

        # get color values for bars
        data_color_norm = [0.5 * (1 + i / (np.ceil(max_score))) * 255 for i in enrichment_scores]
        data_color_norm_8bit = [int(i) if i != np.inf and i != -np.inf else np.sign(i) * max(np.abs(scores_no_inf)) for
                                i in data_color_norm]
        my_cmap = plt.cm.get_cmap('coolwarm')
        colors = my_cmap(data_color_norm_8bit)

        # generate bar plot
        fig, ax = plt.subplots(tight_layout=True, figsize=figsize)

        args = (ax, range(len(enrichment_names)), enrichment_scores) if self.plot_style == 'bar' else (
            ax, range(len(enrichment_names)), 0, enrichment_scores)
        kwargs = dict(linewidth=1, edgecolor='black') if self.plot_style == 'bar' else dict(linewidth=5)

        graph = bar_func(*args, color=colors, zorder=2, **kwargs)
        graph.tick_labels = enrichment_names

        if self.plot_style == 'lollipop':
            x, y = (enrichment_scores, range(len(enrichment_names))) if self.plot_horizontal else (range(
                len(enrichment_names)), enrichment_scores)
            ax.scatter(x, y, color=colors, zorder=3, s=dot_scaling_func(enrichment_obs))

            legend_sizes = [Size(i) for i in [10, 50, 250, 500]]
            ax.legend(legend_sizes, [f"observed={i.size}" for i in legend_sizes],
                      labelspacing=3.2, borderpad=2.5, draggable=True,
                      handler_map={size: sizeHandler(size.size) for size in legend_sizes}, loc='lower right')

        # determine bounds, and enlarge the bound by a small margin (0.5%) so nothing gets cut out of the figure
        if ylim == 'auto':
            bounds = np.array([np.floor(-max_score), np.ceil(max_score)]) * 1.005
        else:
            bounds = np.array([-max_score, max_score]) * 1.005
        # add black line at y=0 and grey lines at every round positive/negative integer in range
        for ind in range(int(bounds[0]), int(bounds[1]) + 1):
            color = 'black' if ind == 0 else 'grey'
            linewidth = 1 if ind == 0 else 0.5
            linestyle = '-' if ind == 0 else '-.'
            line_func(ax, ind, color=color, linewidth=linewidth, linestyle=linestyle, zorder=0)
        # apply xticks
        tick_func(ax, range(len(enrichment_names)))
        ticklabels_func(ax, enrichment_names, **ticklabels_kwargs)
        ax.tick_params(axis='both', which='major', pad=10)
        # title
        ax.set_title(title, fontsize=title_fontsize)
        # add significance asterisks
        for i, (score, sig, obs, exp) in enumerate(
            zip(enrichment_scores, enrichment_pvalue, enrichment_obs, enrichment_exp)):
            asterisks, fontweight = self._get_pval_asterisk(sig, self.alpha)
            if self.plot_horizontal:
                x = score
                y = i
                valign = 'bottom'
                halign = 'center'
                rotation = 270 if np.sign(score) == 1 else 90
                rotation_mode = 'anchor'
            else:
                x = i
                y = score
                valign = 'bottom' if np.sign(score) == 1 else 'top'
                halign = 'center'
                rotation = 0
                rotation_mode = 'default'

            if self.show_expected:
                extra_text = f"{obs}/{exp:.1f}"
                asterisks = f"{extra_text}\n{asterisks}" if self.plot_horizontal or np.sign(
                    score) == 1 else f"{asterisks}\n{extra_text}"
            ax.text(x=x, y=y, s=asterisks, fontname='DejaVu Sans',
                    fontweight=fontweight, rotation=rotation, fontsize=9,
                    horizontalalignment=halign, verticalalignment=valign, zorder=3, rotation_mode=rotation_mode)

        # despine
        _ = [ax.spines[side].set_visible(False) for side in ['top', 'bottom', 'left', 'right']]
        # center bars
        if center_bars:
            if self.plot_horizontal:
                ax.set_xbound(bounds)
                plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            else:
                ax.set_ybound(bounds)
                plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        # add colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(cbar_location, size=f"{5 * (1 + self.plot_horizontal)}%", pad=0.05)
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(*bounds))
        sm.set_array(np.array([]))
        cbar_label_kwargs = dict(label=ylabel, fontsize=ylabel_fontsize, labelpad=15)
        cbar = fig.colorbar(sm, ticks=range(int(bounds[0]), int(bounds[1]) + 1), cax=cax, orientation=cbar_orientation)
        cbar.set_label(**cbar_label_kwargs)
        cax.tick_params(labelsize=14, pad=6)

        if not self.plot_horizontal:
            cax.yaxis.set_ticks_position('left')
            cax.yaxis.set_label_position('left')

        plt.show()
        return fig

    @staticmethod
    def _get_pval_asterisk(pval: float, alpha: float = 0.05):
        fontweight = 'bold'
        if pval > alpha:
            asterisks = 'ns'
            fontweight = 'normal'
        elif pval < 0.0001:
            asterisks = u'\u2217' * 4
        elif pval < 0.001:
            asterisks = u'\u2217' * 3
        elif pval < 0.01:
            asterisks = u'\u2217' * 2
        else:
            asterisks = u'\u2217'
        return asterisks, fontweight


class NonCategoricalEnrichmentRunner(EnrichmentRunner):
    __slots__ = {
        'parametric_test': 'indicates whether to calculate enrichment using a '
                           'parametric or a-parametric statistical test',
        'plot_log_scale': 'indicates whether to plot the results on a logarithmic or linear scale',
        'plot_style': 'indicates the style of histogram to plot the results in',
        'n_bins': 'number of bins in histogram plot of results'}

    def __init__(self, genes: set, attributes: Union[Iterable, str, int], alpha: param_typing.Fraction, biotypes,
                 background_set: set, attr_ref_path: str, biotype_ref_path: str, save_csv: bool, fname: str,
                 return_fig: bool, plot_log_scale: bool, plot_style: Literal['interleaved', 'overlap'],
                 n_bins: int, set_name: str, parallel_backend: Literal[PARALLEL_BACKENDS], parametric_test: bool):

        assert isinstance(plot_log_scale, bool), f"Invalid type for 'plot_log_scale': '{plot_log_scale}'."
        assert plot_style in {'interleaved', 'overlap'}, f"Invalid value for 'plot_style': '{plot_style}'."
        assert isinstance(n_bins,
                          int) and n_bins > 0, f"'n_bins' must be a positive integer. Instead got {type(n_bins)}."

        enrichment_func_name = 't_test' if parametric_test else 'sign_test'
        super().__init__(genes, attributes, alpha, attr_ref_path, True, save_csv, fname, return_fig, True, set_name,
                         parallel_backend, enrichment_func_name, biotypes, background_set, biotype_ref_path,
                         exclude_unannotated_genes=True, single_set=False)
        self.parametric_test = parametric_test
        self.plot_log_scale = plot_log_scale
        self.plot_style = plot_style
        self.n_bins = n_bins

    def _get_enrichment_func(self, pval_func_name: str):
        assert isinstance(pval_func_name, str), f"Invalid type for 'pval_func_name': {type(pval_func_name)}."
        pval_func_name = pval_func_name.lower()
        if pval_func_name == 't_test':
            return self._one_sample_t_test_enrichment
        elif pval_func_name == 'sign_test':
            return self._sign_test_enrichment
        else:
            raise ValueError(f"Unknown enrichment function '{pval_func_name}'.")

    def _sign_test_enrichment(self, attribute: str) -> list:
        exp = self.annotation_df[attribute].median(skipna=True)
        obs = self.annotation_df.loc[list(self.gene_set), attribute].median(skipna=False)
        if np.isnan(obs):
            pval = np.nan
        else:
            _, pval = sign_test(self.annotation_df.loc[list(self.gene_set), attribute].values, exp)
        return [attribute, len(self.gene_set), obs, exp, pval]

    def _one_sample_t_test_enrichment(self, attribute: str) -> list:
        exp = self.annotation_df[attribute].mean(skipna=True)
        obs = self.annotation_df.loc[list(self.gene_set), attribute].mean(skipna=False)
        if np.isnan(obs):
            pval = np.nan
        else:
            _, pval = ttest_1samp(self.annotation_df.loc[list(self.gene_set), attribute], popmean=exp)
        return [attribute, len(self.gene_set), obs, exp, pval]

    def format_results(self, unformatted_results_list: list):
        columns = ['name', 'samples', 'obs', 'exp', 'pval']
        self.results = pd.DataFrame(unformatted_results_list, columns=columns).set_index('name')
        self._correct_multiple_comparisons()
        if self.results['pval'].isna().any():
            warnings.warn(f"One or more of the genes in the background set contained NaN values in "
                          f"{len(self.results['pval'].isna())} of the attributes. "
                          f"P-values and plots will not be generated for those attributes. ")

    def plot_results(self) -> List[plt.Figure]:
        figs = []
        for attribute, padj in zip(self.attributes, self.results['padj']):
            if not np.isnan(padj):
                fig = self.enrichment_histogram(attribute)
                figs.append(fig)
        return figs

    def enrichment_histogram(self, attribute):
        # generate observed and expected Series, either linear or in log10 scale
        exp = self.annotation_df[attribute]
        obs = exp.loc[list(self.gene_set)]
        if self.plot_log_scale:
            xlabel = r"$\log_{10}$" + f"({attribute})"
            exp = np.log10(exp)
            obs = np.log10(obs)
        else:
            xlabel = f"{attribute}"

        # determine bins according to value range and 'n_bins'
        bins = np.linspace(np.min(exp), np.max(exp), self.n_bins).squeeze()

        # generate histogram according to plot style
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        kwargs = dict(bins=bins, density=True, alpha=0.5, edgecolor='black', linewidth=1)
        colors = ['C0', 'C1']
        if self.plot_style.lower() == 'interleaved':
            y, x, _ = ax.hist([exp.values, obs.values], **kwargs, color=colors, label=['Expected', 'Observed'])
            max_y_val = np.max(y)
        elif self.plot_style.lower() == 'overlap':
            y, _, _ = ax.hist(exp.values, **kwargs, color=colors[0], label='Expected')
            y2, _, _ = ax.hist(obs.values, **kwargs, color=colors[1], label='Observed')
            max_y_val = np.max([np.max(y), np.max(y2)])
        else:
            raise NotImplementedError(f"Plot style '{self.plot_style}' not implemented.")

        # set either mean or median as the measure of centrality
        if self.parametric_test:
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
        asterisks, fontweight = self._get_pval_asterisk(self.results.at[attribute, 'padj'], self.alpha)
        ax.vlines([x_exp, x_obs], ymin=max_y_val * 1.12, ymax=max_y_val * 1.16, color='k', linewidth=1)
        ax.hlines(max_y_val * 1.16, xmin=min(x_exp, x_obs), xmax=max(x_exp, x_obs), color='k', linewidth=1)
        ax.text(np.mean([x_exp, x_obs]), max_y_val * 1.17, asterisks, horizontalalignment='center',
                fontweight=fontweight)

        # legend and titles
        ax.legend()
        obs_name = 'Observed' if self.set_name == '' else self.set_name
        ax.set_title(f"Histogram of {attribute} - {obs_name} vs Expected", fontsize=17)
        ax.set_ylabel("Probability density", fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        plt.show()

        return fig


class KEGGEnrichmentRunner(EnrichmentRunner):
    __slots__ = {'organism': 'the organism name for which to fetch GO Annotations',
                 'taxon_id': 'NCBI Taxon ID for which to fetch GO Annotations',
                 'gene_id_type': 'the type of gene ID index that is used',
                 'return_nonsignificant': 'indicates whether to return results which were not found to be '
                                          'statistically significant after enrichment analysis',
                 'attributes_set': 'set of the attributes/KEGG Pathways for which enrichment should be calculated',
                 'pathway_names_dict': 'a dict with KEGG Pathway IDs as keys and their names as values',
                 'plot_pathway_graphs': 'indicates whether to plot pathway graphs of the statistically significant KEGG pathways',
                 'pathway_graphs_format': 'file format for the generated pathway graphs',
                 'gene_id_translator': 'gene ID translator from KEGG into a user-defined gene ID type'
                 }
    KEGG_DF_QUERIES = {}
    printout_params = "have any KEGG Annotations asocciated with them"

    def __init__(self, genes: Union[set, np.ndarray], organism: Union[str, int], gene_id_type: str, alpha: float,
                 return_nonsignificant: bool, save_csv: bool, fname: str, return_fig: bool, plot_horizontal: bool,
                 plot_pathway_graphs: bool, set_name: str,
                 parallel_backend: Literal[PARALLEL_BACKENDS], enrichment_func_name: str,
                 biotypes=None, background_set: set = None, biotype_ref_path: str = None,
                 exclude_unannotated_genes: bool = True, single_set: bool = False, random_seed: int = None,
                 pathway_graphs_format: Literal[GRAPHVIZ_FORMATS] = 'none',
                 plot_style: Literal['bar', 'lollipop'] = 'bar', show_expected: bool = False, **pvalue_kwargs):
        super().__init__(genes, [], alpha, '', return_nonsignificant, save_csv, fname, return_fig, plot_horizontal,
                         set_name, parallel_backend, enrichment_func_name, biotypes, background_set, biotype_ref_path,
                         exclude_unannotated_genes, single_set, random_seed, plot_style, show_expected, **pvalue_kwargs)
        if not self.enrichment_func:
            return
        (self.taxon_id, self.organism), self.gene_id_type = io.get_taxon_and_id_type(organism, gene_id_type,
                                                                                     self.gene_set, 'KEGG')
        self.pathway_names_dict: dict = {}
        self.plot_pathway_graphs = plot_pathway_graphs
        self.pathway_graphs_format = pathway_graphs_format
        self.gene_id_translator = None

    def _get_annotation_iterator(self):
        return io.KEGGAnnotationIterator(self.taxon_id)

    def format_results(self, unformatted_results_list: list):
        columns = ['KEGG ID', 'name', 'samples', 'obs', 'exp', self.en_score_col, 'pval']
        named_results_list = [[entry[0], self.pathway_names_dict[entry[0]]] + entry[1:] for entry in
                              unformatted_results_list]
        self.results = pd.DataFrame(named_results_list, columns=columns).set_index('KEGG ID')
        self._correct_multiple_comparisons()
        # filter non-significant results
        if not self.return_nonsignificant:
            self.results = self.results[self.results['significant']]

    def _correct_multiple_comparisons(self):
        significant, padj = multitest.fdrcorrection(self.results.loc[self.results['pval'].notna(), 'pval'].values,
                                                    alpha=self.alpha, method='negcorr')
        self.results.loc[self.results['pval'].notna(), 'padj'] = padj
        self.results.loc[self.results['padj'].notna(), 'significant'] = significant

    def fetch_annotations(self):
        # check if annotations for the requested query were previously fetched and cached
        query_key = self._get_query_key()
        if query_key in self.KEGG_DF_QUERIES:
            self.annotation_df, self.pathway_names_dict = self.KEGG_DF_QUERIES[query_key]
            return
        else:
            self.annotation_df, self.pathway_names_dict = self._generate_annotation_df()
            # save query results to KEGG_DF_QUERIES
            self.KEGG_DF_QUERIES[query_key] = self.annotation_df, self.pathway_names_dict

    def _generate_annotation_df(self) -> Tuple[pd.DataFrame, Dict[str, str]]:
        # fetch and process KEGG annotations
        sparse_annotation_dict, pathway_name_dict = self._process_annotations()
        print(f"Found annotations for {len(sparse_annotation_dict)} genes.")

        # translate gene IDs
        translated_sparse_annotation_dict = self._translate_gene_ids(sparse_annotation_dict)

        # get boolean DataFrame for enrichment
        annotation_df = parsing.sparse_dict_to_bool_df(translated_sparse_annotation_dict,
                                                       progress_bar_desc="Generating Gene Ontology Referene Table")
        annotation_df[~annotation_df] = np.nan
        return annotation_df, pathway_name_dict

    def _process_annotations(self) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        desc = f"Fetching KEGG annotations for organism '{self.organism}' (taxon ID:{self.taxon_id})"

        sparse_annotation_dict = {}
        pathway_name_dict = {}
        annotation_iter = self._get_annotation_iterator()
        assert annotation_iter.n_annotations > 0, "No KEGG annotations were found for the given parameters. " \
                                                  "Please try again with a different set of parameters. "
        for pathway_id, pathway_name, annotations in tqdm(annotation_iter, desc=desc,
                                                          total=annotation_iter.n_annotations, unit=' annotations'):

            pathway_name_dict[pathway_id] = pathway_name
            # add annotation to annotation dictionary
            for gene_id in annotations:
                if gene_id not in sparse_annotation_dict:
                    sparse_annotation_dict[gene_id] = (set())
                sparse_annotation_dict[gene_id].add(pathway_id)

        return sparse_annotation_dict, pathway_name_dict

    def _translate_gene_ids(self, sparse_annotation_dict: dict):
        source = 'KEGG'
        translated_sparse_annotation_dict = {}
        sparse_dict_cp = sparse_annotation_dict.copy()
        translator = io.GeneIDTranslator(source, self.gene_id_type).run(
            parsing.data_to_tuple(sparse_annotation_dict.keys()))
        self.gene_id_translator = translator
        for gene_id in sparse_annotation_dict:
            if gene_id in translator:
                translated_sparse_annotation_dict[translator[gene_id]] = sparse_dict_cp.pop(gene_id)
        return translated_sparse_annotation_dict

    def _get_query_key(self):
        return self.taxon_id, self.gene_id_type

    def fetch_attributes(self):
        self.attributes = parsing.data_to_list(self.annotation_df.columns)
        self.attributes_set = parsing.data_to_set(self.attributes)

    def plot_results(self) -> List[plt.Figure]:
        figs = [super().plot_results()]
        if self.plot_pathway_graphs:
            for pathway in self.pathway_names_dict:
                if pathway in self.results.index and self.results.loc[pathway, 'significant']:
                    fig = self.pathway_plot(pathway)
                    if fig is None:
                        break
                    figs.append(fig)
        return figs

    def pathway_plot(self, pathway_id: str) -> bool:
        ylabel = self.SINGLE_SET_ENRICHMENT_SCORE_YLABEL if self.single_set else self.ENRICHMENT_SCORE_YLABEL
        pathway_tree = ontology.fetch_kegg_pathway(pathway_id, self.gene_id_translator)
        if self.single_set:
            highlight_genes = self.gene_set
            # TODO: color by rank/ranking metric?
        else:
            highlight_genes = self.gene_set
        return pathway_tree.plot_pathway(highlight_genes, ylabel=ylabel, graph_format=self.pathway_graphs_format)


class GOEnrichmentRunner(EnrichmentRunner):
    __slots__ = {'dag_tree': 'DAG tree containing the hierarchical structure of all GO Terms',
                 'mod_annotation_dfs': "Additional copies of 'annotation_df' which are "
                                       "actively modified by propagation algorithms",
                 'organism': 'the organism name for which to fetch GO Annotations',
                 'taxon_id': 'NCBI Taxon ID for which to fetch GO Annotations',
                 'gene_id_type': 'the type of gene ID index that is used',
                 'propagate_annotations': 'indicates whether to propagate GO Annotations, and with which algorithm',
                 'aspects': 'the GO Aspects for which GO Annotations should be fetched',
                 'evidence_types': 'the evidence types for which GO Annotations should be fetched',
                 'excluded_evidence_types': 'the evidence types for which GO Annotations should NOT be fetched',
                 'databases': 'the ontology databases from which GO Annotations should be fetched',
                 'excluded_databases': 'the ontology databases from which GO Annotations should NOT be fetched',
                 'qualifiers': 'the evidence types for which GO Annotations should be fetched',
                 'excluded_qualifiers': 'the evidence types for which GO Annotations should NOT be fetched',
                 'return_nonsignificant': 'indicates whether to return results which were not found to be '
                                          'statistically significant after enrichment analysis',
                 'plot_ontology_graph': 'indicates whether to plot ontology graph of the statistically significant GO Terms',
                 'ontology_graph_format': 'file format for the generated ontology graph',
                 'attributes_set': 'set of the attributes/GO Terms for which enrichment should be calculated'}
    printout_params = "have any GO Annotations asocciated with them"
    GOA_DF_QUERIES = {}

    def __init__(self, genes: Union[set, np.ndarray], organism: Union[str, int], gene_id_type: str, alpha: float,
                 propagate_annotations: str, aspects: Union[str, Iterable[str]],
                 evidence_types: Union[str, Iterable[str]], excluded_evidence_types: Union[str, Iterable[str]],
                 databases: Union[str, Iterable[str]], excluded_databases: Union[str, Iterable[str]],
                 qualifiers: Union[str, Iterable[str]], excluded_qualifiers: Union[str, Iterable[str]],
                 return_nonsignificant: bool, save_csv: bool, fname: str, return_fig: bool, plot_horizontal: bool,
                 plot_ontology_graph: bool, set_name: str, parallel_backend: Literal[PARALLEL_BACKENDS],
                 enrichment_func_name: str, biotypes=None, background_set: set = None, biotype_ref_path: str = None,
                 exclude_unannotated_genes: bool = True, single_set: bool = False, random_seed: int = None,
                 ontology_graph_format: Literal[GRAPHVIZ_FORMATS] = 'none',
                 plot_style: Literal['bar', 'lollipop'] = 'bar',
                 show_expected: bool = False, **pvalue_kwargs):

        self.propagate_annotations = propagate_annotations.lower()
        super().__init__(genes, [], alpha, '', return_nonsignificant, save_csv, fname, return_fig, plot_horizontal,
                         set_name, parallel_backend, enrichment_func_name, biotypes, background_set, biotype_ref_path,
                         exclude_unannotated_genes, single_set, random_seed, plot_style, show_expected, **pvalue_kwargs)
        if not self.enrichment_func:
            return
        self.dag_tree: ontology.DAGTree = ontology.fetch_go_basic()
        self.mod_annotation_dfs: Tuple[pd.DataFrame, ...] = tuple()
        (self.taxon_id, self.organism), self.gene_id_type = io.get_taxon_and_id_type(organism, gene_id_type,
                                                                                     self.gene_set, 'UniProtKB')
        self.aspects = aspects
        self.evidence_types = evidence_types
        self.excluded_evidence_types = excluded_evidence_types
        self.databases = databases
        self.excluded_databases = excluded_databases
        self.qualifiers = qualifiers
        self.excluded_qualifiers = excluded_qualifiers
        self.plot_ontology_graph = plot_ontology_graph
        self.ontology_graph_format = ontology_graph_format
        self.attributes_set: set = set()

    def _get_enrichment_func(self, pval_func_name: str):
        enrichment_func = super()._get_enrichment_func(pval_func_name)
        if self.propagate_annotations in {'weight'} and pval_func_name.lower() != 'fisher':
            raise NotImplementedError("The 'weight' propagation algorithm is only compatible with Fisher's Exact test.")
        elif self.propagate_annotations in {'allm'} and pval_func_name.lower() != 'fisher':
            warnings.warn(f"The 'weight' propagation algorithm is only compatible with Fisher's Exact test. "
                          f"Therefore, when calculating 'allm' p-values, the 'weight' method will use "
                          f"Fisher's Exact test, while the rest of the methods will use the '{pval_func_name}' method.")
        return enrichment_func

    def fetch_annotations(self):
        # check if annotations for the requested query were previously fetched and cached
        query_key = self._get_query_key()
        if query_key in self.GOA_DF_QUERIES:
            self.annotation_df = self.GOA_DF_QUERIES[query_key]
            return
        else:
            self.annotation_df = self._generate_annotation_df()
            # save query results to GOA_DF_QUERIES
            self.GOA_DF_QUERIES[query_key] = self.annotation_df

    def _generate_annotation_df(self) -> pd.DataFrame:
        # fetch and process GO annotations
        sparse_annotation_dict, source_to_gene_id_dict = self._process_annotations()
        print(f"Found annotations for {len(sparse_annotation_dict)} genes.")

        # translate gene IDs
        translated_sparse_annotation_dict = self._translate_gene_ids(sparse_annotation_dict, source_to_gene_id_dict)

        # get boolean DataFrame for enrichment
        annotation_df = parsing.sparse_dict_to_bool_df(translated_sparse_annotation_dict,
                                                       progress_bar_desc="Generating Gene Ontology Referene Table")
        return annotation_df

    def _process_annotations(self) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
        if self.propagate_annotations != 'no':
            desc = f"Fetching and propagating GO annotations for organism '{self.organism}' (taxon ID:{self.taxon_id})"
        else:
            desc = f"Fetching GO annotations for organism '{self.organism}' (taxon ID:{self.taxon_id})"

        sparse_annotation_dict = {}
        source_to_gene_id_dict = {}
        annotation_iter = self._get_annotation_iterator()
        assert annotation_iter.n_annotations > 0, "No GO annotations were found for the given parameters. " \
                                                  "Please try again with a different set of parameters. "
        for annotation in tqdm(annotation_iter, desc=desc, total=annotation_iter.n_annotations, unit=' annotations'):
            # extract gene_id, go_id, source from the annotation. skip annotations that don't appear in the GO DAG

            if annotation['annotation_class'] not in self.dag_tree:
                continue
            go_id: str = self.dag_tree[annotation['annotation_class']].id  # replace alt_ids with their main id
            gene_id: str = annotation['bioentity_internal_id']
            source: str = annotation['source']

            # add annotation to annotation dictionary
            if gene_id not in sparse_annotation_dict:
                sparse_annotation_dict[gene_id] = (set())
            sparse_annotation_dict[gene_id].add(go_id)

            # add gene id and source to source dict
            if source not in source_to_gene_id_dict:
                source_to_gene_id_dict[source] = set()
            source_to_gene_id_dict[source].add(gene_id)
            # propagate annotations
            self._propagate_annotation(gene_id, go_id, sparse_annotation_dict)
        return sparse_annotation_dict, source_to_gene_id_dict

    def _get_annotation_iterator(self):
        return io.GOlrAnnotationIterator(self.taxon_id, self.aspects,
                                         self.evidence_types, self.excluded_evidence_types,
                                         self.databases, self.excluded_databases,
                                         self.qualifiers, self.excluded_qualifiers)

    def _propagate_annotation(self, gene_id: str, go_id: str, sparse_annotation_dict: dict):
        if self.propagate_annotations != 'no':
            parents = set(self.dag_tree.upper_induced_graph_iter(go_id))
            sparse_annotation_dict[gene_id].update(parents)

    def _translate_gene_ids(self, sparse_annotation_dict: dict, source_to_gene_id_dict: dict):
        translated_sparse_annotation_dict = {}
        sparse_dict_cp = sparse_annotation_dict.copy()
        for source in source_to_gene_id_dict:
            try:
                translator = io.GeneIDTranslator(source, self.gene_id_type).run(source_to_gene_id_dict[source])
                for gene_id in sparse_dict_cp.copy():
                    if gene_id in translator:
                        translated_sparse_annotation_dict[translator[gene_id]] = sparse_dict_cp.pop(gene_id)

            except AssertionError as e:
                if 'not a valid Uniprot Dataset' in "".join(e.args):
                    warnings.warn(f"Failed to map gene IDs for {len(source_to_gene_id_dict[source])} annotations "
                                  f"from dataset '{source}'.")
                else:
                    raise e
        return translated_sparse_annotation_dict

    def _get_query_key(self):
        return (self.taxon_id, self.gene_id_type, parsing.data_to_tuple(self.aspects, sort=True),
                parsing.data_to_tuple(self.evidence_types, sort=True),
                parsing.data_to_tuple(self.excluded_evidence_types, sort=True),
                parsing.data_to_tuple(self.databases, sort=True),
                parsing.data_to_tuple(self.excluded_databases, sort=True),
                parsing.data_to_tuple(self.qualifiers, sort=True),
                parsing.data_to_tuple(self.excluded_qualifiers, sort=True),
                self.propagate_annotations != 'no')

    def fetch_attributes(self):
        self.attributes = parsing.data_to_list(self.annotation_df.columns)
        self.attributes_set = parsing.data_to_set(self.attributes)

    def _correct_multiple_comparisons(self):
        significant, padj = multitest.fdrcorrection(self.results.loc[self.results['pval'].notna(), 'pval'].values,
                                                    alpha=self.alpha, method='negcorr')
        self.results.loc[self.results['pval'].notna(), 'padj'] = padj
        self.results.loc[self.results['padj'].notna(), 'significant'] = significant

    def plot_results(self) -> List[plt.Figure]:
        n_bars = min(10, len(self.results))
        kwargs = dict()
        if self.single_set:
            kwargs['ylabel'] = self.SINGLE_SET_ENRICHMENT_SCORE_YLABEL
            kwargs['title'] = f"Single-list GO enrichment for gene set '{self.set_name}'\n" \
                              f"top {n_bars} most specific GO terms"
        else:
            kwargs['ylabel'] = self.ENRICHMENT_SCORE_YLABEL
            kwargs['title'] = f"GO enrichment for gene set '{self.set_name}'\n" \
                              f"top {n_bars} most specific GO terms"

        figs = [self.enrichment_bar_plot(n_bars=n_bars, **kwargs)]
        if self.plot_ontology_graph:
            ontology_kwargs = kwargs.copy()
            ontology_kwargs['title'] = f"Single-list GO enrichment for gene set '{self.set_name}'" if self.single_set \
                else f"GO enrichment for gene set '{self.set_name}'"
            figs.extend(self.go_dag_plot(**ontology_kwargs))
        return figs

    def go_dag_plot(self, title, dpi: int = 300, ylabel: str = r"$\log_2$(Fold Enrichment)") -> List[plt.Figure]:
        aspects = ('biological_process', 'cellular_component',
                   'molecular_function') if self.aspects == 'any' else parsing.data_to_tuple(self.aspects)
        figs = []
        for go_aspect in aspects:
            this_title = title + f'\n{go_aspect}'.replace('_', ' ').title()
            fig = self.dag_tree.plot_ontology(go_aspect, self.results, self.en_score_col, this_title, ylabel,
                                              self.ontology_graph_format, dpi)
            if fig is None:
                return figs
            figs.append(fig)

        return figs

    def format_results(self, unformatted_results_dict: dict):
        columns = ['name', 'samples', 'obs', 'exp', self.en_score_col, 'pval']
        self.results = pd.DataFrame.from_dict(unformatted_results_dict, orient='index', columns=columns)
        self._correct_multiple_comparisons()
        # filter non-significant results
        if not self.return_nonsignificant:
            self.results = self.results[self.results['significant']]
        # sort results by specificity (level in DAG tree)
        self.results['dag_level'] = [self.dag_tree[ind].level for ind in self.results.index]
        self.results = self.results.sort_values('dag_level', ascending=False).drop('dag_level', axis=1).rename_axis(
            'go_id')

    def _calculate_enrichment_serial(self) -> dict:
        desc = f"Calculating enrichment for {len(self.attributes)} GO terms " \
               f"using the '{self.propagate_annotations}' method"
        if self.propagate_annotations == 'classic' or self.propagate_annotations == 'no':
            self.mod_annotation_dfs = (self.annotation_df,)
            result = self._go_classic_pvalues_serial(desc)
        elif self.propagate_annotations == 'elim':
            self.mod_annotation_dfs = (self.annotation_df.copy(deep=True),)
            result = self._go_elim_pvalues_serial(desc)
        elif self.propagate_annotations == 'weight':
            self.mod_annotation_dfs = (self.annotation_df.astype("float64", copy=True),)
            result = self._go_weight_pvalues_serial(desc)
        elif self.propagate_annotations == 'all.m':
            result = self._go_allm_pvalues_serial()
        else:
            raise ValueError(f"invalid propagation method '{self.propagate_annotations}'.")
        return result

    def _calculate_enrichment_parallel(self) -> dict:
        desc = f"Calculating enrichment for {len(self.attributes)} GO terms " \
               f"using the '{self.propagate_annotations}' method"
        if self.propagate_annotations == 'classic' or self.propagate_annotations == 'no':
            self.mod_annotation_dfs = (self.annotation_df,)
            result = self._go_classic_pvalues_parallel(desc)
        elif self.propagate_annotations == 'elim':
            self.mod_annotation_dfs = tuple(
                self.annotation_df.loc[:, self._go_level_iterator(namespace)].copy(deep=True) for namespace in
                self.dag_tree.namespaces)
            result = self._go_elim_pvalues_parallel(desc)
        elif self.propagate_annotations == 'weight':
            self.mod_annotation_dfs = tuple(
                self.annotation_df.loc[:, self._go_level_iterator(namespace)].astype("float64", copy=True) for namespace
                in self.dag_tree.namespaces)
            result = self._go_weight_pvalues_parallel(desc)
        elif self.propagate_annotations == 'all.m':
            result = self._go_allm_pvalues_parallel()
        else:
            raise ValueError(f"invalid propagation method '{self.propagate_annotations}'.")
        return result

    def _go_classic_pvalues_serial(self, progress_bar_desc: str = '') -> dict:
        return self._go_classic_on_batch(tqdm(self.attributes, desc=progress_bar_desc, unit=' GO terms'))

    def _go_classic_pvalues_parallel(self, progress_bar_desc: str = '') -> dict:
        # split GO terms into batches of size 1000 to reduce the overhead of parallel jobs.
        # the GO terms in each batch are calculated serially in a single process
        partitioned = parsing.partition_list(self.attributes, 1000)
        return self._parallel_over_grouping(self._go_classic_on_batch, partitioned,
                                            (0 for _ in range(len(partitioned))), progress_bar_desc=progress_bar_desc)

    def _go_classic_on_batch(self, go_term_batch: Iterable[str], mod_df_index: int = 0):
        """
        calculate stats and p-values with the 'classic' propagation algorithm for an iterable batch of GO terms.

        :param go_term_batch: batch of GO terms to calculate enrichment for
        :type go_term_batch: an iterable of str (GO IDs)
        :param mod_df_index:
        :type mod_df_index: int between 0 and 2 (default 0)
        :return:
        :rtype:
        """
        return {go_id: self.enrichment_func(go_id, mod_df_ind=mod_df_index, **self.pvalue_kwargs) for go_id in
                go_term_batch}

    def _go_elim_pvalues_serial(self, progress_bar_desc: str = '') -> dict:
        result = self._go_elim_on_aspect('all', progress_bar_desc=progress_bar_desc)
        return result

    def _parallel_over_grouping(self, func, grouping: Iterable, mod_df_inds: Iterable[int],
                                max_nbytes: Union[str, None] = '1M', progress_bar_desc: str = '') -> dict:
        assert validation.is_method_of_class(func, type(self))
        result_dicts = generic.ProgressParallel(desc=progress_bar_desc, n_jobs=-1, max_nbytes=max_nbytes,
                                                backend=self.parallel_backend)(
            joblib.delayed(func)(group, ind) for group, ind in zip(grouping, mod_df_inds))
        result = {}
        for d in result_dicts:
            result.update(d)
        return result

    def _go_elim_pvalues_parallel(self, progress_bar_desc: str = '') -> dict:
        go_aspects = parsing.data_to_tuple(self.dag_tree.namespaces)
        # max_nbytes is set to None to prevent joblib from turning the mod_annotation_df into a memory map.
        # since this algorithm modifies the dataframe while running, turning the dataframe into a memory map
        # would raise an exception ("...trying to write into a read-only location...").
        return self._parallel_over_grouping(self._go_elim_on_aspect, go_aspects, range(len(go_aspects)),
                                            max_nbytes=None, progress_bar_desc=progress_bar_desc)

    def _go_elim_on_aspect(self, go_aspect: str, mod_df_ind: int = 0, progress_bar_desc: str = None) -> dict:
        result_dict = {}
        marked_nodes = {}
        go_id_iter = self._go_level_iterator(go_aspect) if progress_bar_desc is None \
            else tqdm(self._go_level_iterator(go_aspect), unit=' GO terms', desc=progress_bar_desc,
                      total=len(self.attributes))
        for go_id in go_id_iter:
            if go_id in marked_nodes:  # if this node was marked, remove from it all marked genes
                self.mod_annotation_dfs[mod_df_ind].loc[list(marked_nodes[go_id]), go_id] = False
            result_dict[go_id] = self.enrichment_func(go_id, mod_df_ind=mod_df_ind, **self.pvalue_kwargs)
            # if current GO ID is significantly ENRICHED, mark its ancestors
            if result_dict[go_id][-1] <= self.alpha and result_dict[go_id][-2] > 0:
                new_marked_genes = set(
                    self.mod_annotation_dfs[mod_df_ind][go_id][self.mod_annotation_dfs[mod_df_ind][go_id]].index)
                for ancestor in self.dag_tree.upper_induced_graph_iter(go_id):
                    if ancestor not in marked_nodes:
                        marked_nodes[ancestor] = set()
                    marked_nodes[ancestor] = marked_nodes[ancestor].union(new_marked_genes)
        return result_dict

    def _go_weight_pvalues_serial(self, progress_bar_desc: str = '') -> dict:
        return self._go_weight_on_aspect('all', progress_bar_desc=progress_bar_desc)

    def _go_weight_pvalues_parallel(self, progress_bar_desc: str = '') -> dict:
        go_aspects = parsing.data_to_tuple(self.dag_tree.namespaces)
        # max_nbytes is set to None to prevent joblib from turning the mod_annotation_df into a memory map.
        # since this algorithm modifies the dataframe while running, turning the dataframe into a memory map
        # would raise an exception ("...trying to write into a read-only location...").
        return self._parallel_over_grouping(self._go_weight_on_aspect, go_aspects, range(len(go_aspects)),
                                            max_nbytes=None, progress_bar_desc=progress_bar_desc)

    def _go_weight_on_aspect(self, go_aspect: str, mod_df_ind: int = 0, progress_bar_desc: str = None) -> dict:
        result_dict = {}
        weights = collections.defaultdict(lambda: 1)  # default weight for all nodes is 1
        go_id_iter = self._go_level_iterator(go_aspect) if progress_bar_desc is None \
            else tqdm(self._go_level_iterator(go_aspect), unit=' GO terms', desc=progress_bar_desc,
                      total=len(self.attributes))
        for go_id in go_id_iter:
            children = {child for child in self.dag_tree[go_id].get_children() if child in self.attributes_set}
            self._compute_term_sig(mod_df_ind, go_id, children, weights, result_dict)
        return result_dict

    def _go_allm_pvalues_serial(self):
        outputs = {}
        try:
            for method in ['classic', 'elim', 'weight']:
                self.propagate_annotations = method
                outputs[method] = self._calculate_enrichment_serial()
        finally:
            self.propagate_annotations = 'all.m'  # make sure 'propagate_annotations' is correctly set to its
            # original value if something goes wrong in runtime
        return self._calculate_allm(outputs)

    def _go_allm_pvalues_parallel(self):
        # TODO: progress bar, desc
        outputs = {}
        try:
            for method in ['classic', 'elim', 'weight']:
                self.propagate_annotations = method
                outputs[method] = self._calculate_enrichment_parallel()
        finally:
            self.propagate_annotations = 'all.m'  # make sure 'propagate_annotations' is correctly set to its
            # original value if something goes wrong in runtime
        return self._calculate_allm(outputs)

    def _calculate_allm(self, outputs: Dict[str, Dict[str, list]]) -> Dict[str, tuple]:
        """
        Implementation of 'combining algorithms' (section 2.4) from Alexa et al 2006: \
        https://academic.oup.com/bioinformatics/article/22/13/1600/193669

        :param outputs: a dictionary containing output dictionaries from all 3 propagation algorithms \
        (classic, elim, weight)
        :type outputs: Dict[str, Dict[str, list]]
        :return: a combined output dictionary containing p-values that were averaged on a logarithmic scale.
        :rtype: Dict[str, tuple]
        """
        result = {}
        for go_id in self.attributes:
            pvalue = np.exp(np.mean(
                np.log([outputs['classic'][go_id][-1], outputs['elim'][go_id][-1], outputs['weight'][go_id][-1]])))
            result[go_id] = (*outputs['classic'][go_id][:-1], pvalue)
        return result

    def _go_level_iterator(self, namespace: str = 'all'):
        for go_id in self.dag_tree.level_iter(namespace):
            if go_id in self.attributes_set:  # skip GO IDs that have no annotations whatsoever (direct or inherited)
                yield go_id

    def _compute_term_sig(self, mod_df_ind: int, go_id: str, children: set, weights: dict, result: dict,
                          tolerance: float = 10 ** -50):
        """
        Implementation of the computeTermSig(u, children) function from Alexa et al 2006: \
        https://doi.org/10.1093/bioinformatics/btl140

        :param mod_df_ind: index of the annotation dataframe to be used for the algorithm \
        (in case the algorithm runs in parallel and multiple dataframes are used)
        :type mod_df_ind: int
        :param go_id: the current GO ID (node/'u') on which the function is calculated
        :type go_id: str
        :param children: a set of the children of go_id
        :type children: set
        :param weights: a dictionary of weight multipliers for each GO ID/node to be used in p-value calculations
        :type weights: dict
        :param result: a dictionary containing the results of the enrichment analysis (p-value, fold-change, etc')
        :type result: dict
        :param tolerance: absolute tolerance for float equality in the algorithm
        :type tolerance: float

        """
        # calculate stats OR update p-value for go_id
        if go_id in result:
            result[go_id][-1] = self._calc_fisher_pval(
                *self._get_hypergeometric_parameters(go_id, mod_df_ind=mod_df_ind))
        else:
            result[go_id] = self._fisher_enrichment(go_id, mod_df_ind=mod_df_ind)
        # if go_id has no children left to compare to, stop the calculation here
        if len(children) == 0:
            return
        # calculate ratio of p-values between go_id and its children. if ratio>=1, child is more significant than go_id
        sig_children = set()
        for child in children:
            a = result[child][-1]
            b = result[go_id][-1]
            if np.isclose(a, 0, rtol=0, atol=tolerance):  # if a is close to 0, add tolerance to avoid division by 0
                sig_ratio = b / (a + tolerance)
            elif np.isclose(a, b, rtol=0, atol=tolerance):  # if a and b are nearly identical, set their ratio to 1
                sig_ratio = 1
            else:
                sig_ratio = b / a

            weights[child] = sig_ratio
            if weights[child] >= 1:
                sig_children.add(child)

        # CASE 1: if go_id is more significant than all children, re-weigh the children and recompute their stats
        if len(sig_children) == 0:
            for child in children:
                self.mod_annotation_dfs[mod_df_ind].loc[self.annotation_df[go_id], child] *= weights[child]
                result[child][-1] = self._calc_fisher_pval(
                    *self._get_hypergeometric_parameters(child, mod_df_ind=mod_df_ind))
            return

        # CASE 2: if some children are more significant than parent, re-weigh ancesctors (including 'go_id'),
        # and then recompute stats for go_id
        inclusive_ancestors = {ancestor for ancestor in
                               itertools.chain([go_id], self.dag_tree.upper_induced_graph_iter(go_id)) if
                               ancestor in self.attributes_set}
        for sig_child in sig_children:
            for inclusive_ancestor in inclusive_ancestors:
                self.mod_annotation_dfs[mod_df_ind].loc[self.annotation_df[sig_child], inclusive_ancestor] *= (
                    1 / weights[sig_child])
        # re-run compute_term_sig, only with the children which were not more significant than their parents
        self._compute_term_sig(mod_df_ind, go_id, children.difference(sig_children), weights, result, tolerance)

    def _randomization_enrichment(self, go_id: str, reps: int, mod_df_ind: int = None) -> list:
        mod_df_ind = 0 if mod_df_ind is None else mod_df_ind
        go_name = self.dag_tree[go_id].name
        bg_df = self.mod_annotation_dfs[mod_df_ind][go_id]
        bg_size = self.annotation_df.shape[0]
        n = len(self.gene_set)
        expected_fraction = self.annotation_df[go_id].sum() / bg_size
        observed_fraction = self.annotation_df.loc[list(self.gene_set), go_id].sum() / n
        mod_observed_fraction = bg_df.loc[list(self.gene_set)].sum() / n
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_randomization_pval(n, log2_fold_enrichment, bg_df.values, reps, mod_observed_fraction)
        return [go_name, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval]

    def _xlmhg_enrichment(self, go_id: str, mod_df_ind: int = None, **_) -> list:
        go_name = self.dag_tree[go_id].name
        index_vec, rev_index_vec = self._generate_xlmhg_index_vectors(go_id, mod_df_ind)
        n, X, L, table = self._get_xlmhg_parameters(index_vec)

        res_obj_fwd = xlmhglite.get_xlmhg_test_result(N=n, indices=index_vec, X=X, L=L, table=table)
        res_obj_rev = xlmhglite.get_xlmhg_test_result(N=n, indices=rev_index_vec, X=X, L=L, table=table)

        obs, exp, en_score, pval = self._extract_xlmhg_results(res_obj_fwd, res_obj_rev)
        return [go_name, n, obs, exp, en_score, pval]

    def _generate_xlmhg_index_vectors(self, attribute: str, mod_df_ind: int = None) -> Tuple[np.ndarray, np.ndarray]:
        n = len(self.ranked_genes)
        ranked_srs = self.mod_annotation_dfs[mod_df_ind].loc[self.ranked_genes, attribute]
        assert ranked_srs.shape[0] == len(self.ranked_genes)
        index_vec = np.uint16(np.nonzero(ranked_srs.values)[0])
        rev_index_vec = np.uint16([n - 1 - index_vec[i - 1] for i in range(len(index_vec), 0, -1)])
        return index_vec, rev_index_vec

    def _hypergeometric_enrichment(self, go_id: str, mod_df_ind: int = None) -> list:
        bg_size, de_size, go_size, go_de_size = self._get_hypergeometric_parameters(go_id, mod_df_ind=mod_df_ind)

        go_name = self.dag_tree[go_id].name
        expected_fraction = self.annotation_df[go_id].sum() / bg_size
        observed_fraction = self.annotation_df.loc[list(self.gene_set), go_id].sum() / de_size
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_hypergeometric_pval(bg_size=bg_size, de_size=de_size, go_size=go_size, go_de_size=go_de_size)
        obs, exp = int(de_size * observed_fraction), de_size * expected_fraction
        return [go_name, de_size, obs, exp, log2_fold_enrichment, pval]

    def _fisher_enrichment(self, go_id: str, mod_df_ind: int = None) -> list:
        bg_size, de_size, go_size, go_de_size = self._get_hypergeometric_parameters(go_id, mod_df_ind)

        expected_fraction = self.annotation_df[go_id].sum() / bg_size
        observed_fraction = self.annotation_df.loc[list(self.gene_set), go_id].sum() / de_size
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_fisher_pval(bg_size=bg_size, de_size=de_size, go_size=go_size, go_de_size=go_de_size)
        obs, exp = int(de_size * observed_fraction), de_size * expected_fraction
        return [self.dag_tree[go_id].name, de_size, obs, exp, log2_fold_enrichment, pval]

    def _get_hypergeometric_parameters(self, go_id: str, mod_df_ind: int = None) -> Tuple[int, int, int, int]:
        bg_size = self.mod_annotation_dfs[mod_df_ind].shape[0]
        de_size = len(self.gene_set)
        go_size = int(np.ceil(self.mod_annotation_dfs[mod_df_ind][go_id].sum()))
        go_de_size = int(np.ceil(self.mod_annotation_dfs[mod_df_ind].loc[list(self.gene_set), go_id].sum()))
        return bg_size, de_size, go_size, go_de_size
