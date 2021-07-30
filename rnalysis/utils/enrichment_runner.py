import collections
import itertools
import logging
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Iterable, List, Tuple, Union, Collection, Set, Dict

import joblib
import matplotlib.pyplot as plt
import numba
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as multitest
from matplotlib.cm import ScalarMappable
from scipy.stats import hypergeom, ttest_1samp, fisher_exact
from statsmodels.stats.descriptivestats import sign_test
from tqdm.auto import tqdm
from xlmhg import get_xlmhg_test_result as xlmhg_test

from rnalysis.utils import ontology, io, parsing, ref_tables, validation, generic

logging.getLogger('xlmhg').setLevel(50)  # suppress warnings from xlmhg module


class EnrichmentRunner:
    __slots__ = {'results': 'DataFrame containing enrichment analysis results',
                 'annotation_df': 'DataFrame containing all annotation data per gene',
                 'gene_set': 'the set of genes/genomic features whose enrichment to calculate',
                 'attributes': 'the list of attributes/terms to calculate enrichment for',
                 'alpha': 'the statistical signifiacnce threshold',
                 'attr_ref_path': 'path of the Attribute Reference Table to load, if such table exists',
                 'save_csv': 'indicates whether the results should be saved to a csv file',
                 'fname': 'name of the file to save results to',
                 'return_fig': 'indicates whether to return a matplotlib Figure after plotting the results',
                 'plot_horizontal': 'indicates whether to plot the results horizontally or vertically',
                 'set_name': 'name of the set of genes/genomic features whose enrichment is calculated',
                 'parallel': 'indicates whether to calculate the enrichment using parallel processing',
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
                                 'pre-sorted and ranked by the user'}
    printout_params = "appear in the Attribute Reference Table"

    def __init__(self, genes: Union[set, np.ndarray], attributes: Union[Iterable, str, int], alpha: float,
                 attr_ref_path: str, save_csv: bool, fname: str, return_fig: bool, plot_horizontal: bool, set_name: str,
                 parallel: bool, enrichment_func_name: str, biotypes=None, background_set: set = None,
                 biotype_ref_path: str = None, single_set: bool = False, random_seed: int = None, **pvalue_kwargs):
        self.results: pd.DataFrame = pd.DataFrame()
        self.annotation_df: pd.DataFrame = pd.DataFrame()
        self.gene_set = parsing.data_to_set(genes)
        self.attributes = attributes
        self.alpha = alpha
        self.attr_ref_path = ref_tables.get_attr_ref_path(attr_ref_path)
        self.save_csv = save_csv

        if self.save_csv:
            if fname is None:
                self.fname = input("Please insert the full name and path to save the file to")
            else:
                assert isinstance(fname, (str, Path))
            self.fname = str(fname)
        self.return_fig = return_fig
        self.plot_horizontal = plot_horizontal
        self.set_name = set_name
        self.parallel = parallel
        self.enrichment_func = self._get_enrichment_func(enrichment_func_name)
        self.pvalue_kwargs = pvalue_kwargs
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
            self.biotype_ref_path = ref_tables.get_biotype_ref_path(biotype_ref_path) if biotypes != 'all' else None
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

    def run(self) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        self.fetch_annotations()
        self.fetch_attributes()
        if not self.single_set:
            self.get_background_set()
        self.update_gene_set()
        self.filter_annotations()
        unformatted_results = self.calculate_enrichment()
        self.format_results(unformatted_results)
        fig = self.plot_results()
        if self.save_csv:
            self.results_to_csv()
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
        io.save_csv(self.results, filename=self.fname if self.fname.endswith('.csv') else self.fname + '.csv')

    def get_background_set(self):
        if self.background_set is None:
            self._get_background_set_from_biotype()
        else:
            self._get_background_set_from_set()

        orig_bg_size = len(self.background_set)
        self.background_set = self.background_set.intersection(self.annotation_df.index)
        if orig_bg_size - len(self.background_set) > 0:
            warnings.warn(
                f"{orig_bg_size - len(self.background_set)} genes out of the requested {orig_bg_size} background genes "
                f"do not {self.printout_params}, and are therefore ignored."
                f" \nThis leaves a total of {len(self.background_set)} background genes.")

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
            return self._xlmhg_enrichment
        else:
            raise ValueError(f"Unknown enrichment function '{pval_func_name}'.")

    def _get_hypergeometric_parameters(self, attribute: str) -> Tuple[int, int, int, int]:
        bg_size = self.annotation_df.shape[0]
        de_size = len(self.gene_set)
        go_size = self.annotation_df[attribute].notna().sum()
        go_de_size = self.annotation_df.loc[self.gene_set, attribute].notna().sum()
        return bg_size, de_size, go_size, go_de_size

    def _xlmhg_enrichment(self, attribute: str) -> list:
        n = len(self.ranked_genes)
        index_vec, rev_index_vec = self._xlmhg_index_vectors(attribute)

        # X = the minimal amount of 'positive' elements above the hypergeometric cutoffs out of all of the positive
        # elements in the ranked set. Determined to be the minimum between x_min and ceil(x_frac * k),
        # where 'k' is the number of 'positive' elements in the ranked set.
        x_frac = 0.5
        x_min = 10
        X = min(x_min, int(np.ceil(x_frac * len(index_vec)))) if 'X' not in self.pvalue_kwargs \
            else self.pvalue_kwargs['X']
        # L = the lowest possible cutoff (n) to be tested out of the entire list.
        # Determined to be floor(l_frac * N), where 'N' is total number of elements in the ranked set (n).
        l_frac = 0.1
        L = int(np.floor(l_frac * n)) if 'L' not in self.pvalue_kwargs else self.pvalue_kwargs['L']
        # pre-allocate empty array to speed up computation
        table = np.empty((len(index_vec) + 1, n - len(index_vec) + 1), dtype=np.longdouble)
        res_obj_fwd = xlmhg_test(N=n, indices=index_vec, L=L, X=X, table=table)
        res_obj_rev = xlmhg_test(N=n, indices=rev_index_vec, L=L, X=X, table=table)

        if res_obj_fwd.pval <= res_obj_rev.pval:
            pval, en_score = res_obj_fwd.pval, res_obj_fwd.escore
        else:
            pval, en_score = res_obj_rev.pval, 1 / res_obj_rev.escore
        pval = pval if not np.isnan(pval) else 1
        en_score = en_score if not np.isnan(en_score) else 1
        log2_en_score = np.log2(en_score) if en_score > 0 else -np.inf

        return [attribute, n, log2_en_score, pval]

    def _xlmhg_index_vectors(self, attribute) -> Tuple[np.ndarray, np.ndarray]:
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
        obs_array = self.annotation_df.loc[self.gene_set, attribute].notna().values
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
    @numba.jit(nopython=True)
    def _calc_randomization_pval(n: int, log2fc: float, bg_array: np.ndarray, reps: int, obs_frac: float) -> float:
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
        if go_de_size / de_size < go_size / bg_size:
            return hypergeom.cdf(go_de_size, bg_size, go_size, de_size)
        return hypergeom.sf(go_de_size - 1, bg_size, go_size, de_size)

    def _get_background_set_from_biotype(self):
        if self.biotypes == 'all':
            self.background_set = parsing.data_to_set(self.annotation_df.index)
        else:
            biotype_ref_df = io.load_csv(self.biotype_ref_path)
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
            updated_gene_set = self.gene_set.intersection(self.annotation_df.index)
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
                warnings.warn(f"{not_in_bg} genes in the enrichment set do not {self.printout_params} "
                              f"and/or do not appear in the background gene set. \n"
                              f"Enrichment will be computed on the remaining {len(self.gene_set)} genes.")

    def fetch_annotations(self):
        self.annotation_df = io.load_csv(self.attr_ref_path)
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
            self.annotation_df = self.annotation_df.loc[self.background_set, self.attributes].sort_index()

    def calculate_enrichment(self) -> list:
        self.set_random_seed()
        if self.parallel:
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
        result = generic.ProgressParallel(n_jobs=-1, desc="Calculating enrichment", unit='attribute')(
            joblib.delayed(self.enrichment_func)(attribute, **self.pvalue_kwargs) for attribute in self.attributes)
        return result

    def format_results(self, unformatted_results_list: list):
        if self.single_set:
            columns = ['name', 'samples', self.en_score_col, 'pval']
        else:
            columns = ['name', 'samples', 'obs', 'exp', self.en_score_col, 'pval']
        self.results = pd.DataFrame(unformatted_results_list, columns=columns).set_index('name')
        self._correct_multiple_comparisons()

    def _correct_multiple_comparisons(self):
        significant, padj = multitest.fdrcorrection(self.results.loc[self.results['pval'].notna(), 'pval'].values,
                                                    alpha=self.alpha)
        self.results.loc[self.results['pval'].notna(), 'padj'] = padj
        self.results['significant'] = False  # set default value as False
        self.results.loc[self.results['padj'].notna(), 'significant'] = significant

    def plot_results(self) -> plt.Figure:
        if self.single_set:
            return self.enrichment_bar_plot(ylabel=r"$\log_2$(XL-mHG enrichment score)",
                                            title=f"Single-list enrichment for {self.set_name}")
        return self.enrichment_bar_plot(title=f"Enrichment for {self.set_name}")

    def enrichment_bar_plot(self, n_bars: int = 'all', name_col: str = None, center_bars: bool = True,
                            ylabel: str = r"$\log_2$(Fold Enrichment)",
                            title: str = 'Enrichment results') -> plt.Figure:

        """
        Receives a DataFrame output from an enrichment function and plots it in a bar plot. \
        For the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :param title:
        :type title:
        :param n_bars:
        :type n_bars:
        :param name_col:
        :type name_col:
        :param ylabel: plot ylabel.
        :type ylabel: str
        :param center_bars: if True, centers the bars around Y=0. Otherwise, ylim is determined by min/max values.
        :type center_bars: bool (default True)
        :return: Figure object containing the bar plot
        :rtype: matplotlib.figure.Figure instance
        """
        plt.style.use('seaborn-white')
        # determine number of entries/bars to plot
        if n_bars != 'all':
            assert isinstance(n_bars, int) and n_bars > 0
            results = self.results.iloc[:n_bars]
        else:
            results = self.results
        # pull names/scores/pvals out to avoid accidentally changing the results DataFrame in-place
        enrichment_names = results.index.values.tolist() if name_col is None else results[name_col].values.tolist()
        enrichment_scores = results[self.en_score_col].values.tolist()
        enrichment_pvalue = results['padj'].values.tolist()

        # choose functions and parameters according to the graph's orientation (horizontal vs vertical)
        if self.plot_horizontal:
            figsize = [14, 0.4 * (6.4 + self.results.shape[0])]
            bar_func = plt.Axes.barh
            line_func = plt.Axes.axvline
            cbar_kwargs = dict(location='bottom')
            tick_func = plt.Axes.set_yticks
            ticklabels_func = plt.Axes.set_yticklabels
            ticklabels_kwargs = dict(fontsize=13, rotation=0)
            for lst in (enrichment_names, enrichment_scores, enrichment_pvalue):
                lst.reverse()
        else:
            figsize = [0.5 * (6.4 + self.results.shape[0]), 5.6]
            bar_func = plt.Axes.bar
            line_func = plt.Axes.axhline
            cbar_kwargs = dict(location='left')
            tick_func = plt.Axes.set_xticks
            ticklabels_func = plt.Axes.set_xticklabels
            ticklabels_kwargs = dict(fontsize=13, rotation=45)

        # set enrichment scores which are 'inf' or '-inf' to be the second highest/lowest enrichment score in the list
        scores_no_inf = [i for i in enrichment_scores if i != np.inf and i != -np.inf and i < 0]
        if len(scores_no_inf) == 0:
            scores_no_inf.append(-1)
        for i in range(len(enrichment_scores)):
            if enrichment_scores[i] == -np.inf:
                enrichment_scores[i] = min(scores_no_inf)
        max_score = max(np.max(np.abs(enrichment_scores)), 2)

        # get color values for bars
        data_color_norm = [0.5 * (1 + i / (np.floor(max_score) + 1)) * 255 for i in enrichment_scores]
        data_color_norm_8bit = [int(i) if i != np.inf and i != -np.inf else np.sign(i) * max(np.abs(scores_no_inf)) for
                                i in data_color_norm]
        my_cmap = plt.cm.get_cmap('coolwarm')
        colors = my_cmap(data_color_norm_8bit)

        # generate bar plot
        fig, ax = plt.subplots(constrained_layout=True, figsize=figsize)
        bar = bar_func(ax, range(len(enrichment_names)), enrichment_scores, color=colors, edgecolor='black',
                       linewidth=1, zorder=2)
        bar.tick_labels = enrichment_names
        # determine bound, and enlarge the bound by a small margin (0.2%) so nothing gets cut out of the figure
        bounds = np.array([np.ceil(-max_score) - 1, (np.floor(max_score) + 1)]) * 1.002
        # add black line at y=0 and grey lines at every round positive/negative integer in range
        for ind in range(int(bounds[0]) + 1, int(bounds[1]) + 1):
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
            asterisks, fontweight = self._get_pval_asterisk(sig, self.alpha)
            if self.plot_horizontal:
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
                    fontsize=9, horizontalalignment=halign, verticalalignment=valign, zorder=1)
        # despine
        _ = [ax.spines[side].set_visible(False) for side in ['top', 'right']]
        # center bars
        if center_bars:
            if self.plot_horizontal:
                ax.set_xbound(bounds)
                plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            else:
                ax.set_ybound(bounds)
                plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

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

    def __init__(self, genes: set, attributes: Union[Iterable, str, int], alpha: float, biotypes, background_set: set,
                 attr_ref_path: str, biotype_ref_path: str, save_csv: bool, fname: str,
                 return_fig: bool, plot_log_scale: bool, plot_style: str, n_bins: int, set_name: str,
                 parallel: bool, parametric_test: bool):

        assert isinstance(plot_log_scale, bool), f"Invalid type for 'plot_log_scale': '{plot_log_scale}'."
        assert plot_style in {'interleaved', 'overlap'}, f"Invalid value for 'plot_style': '{plot_style}'."
        assert isinstance(n_bins,
                          int) and n_bins > 0, f"'n_bins' must be a positive integer. Instead got {type(n_bins)}."

        enrichment_func_name = 't_test' if parametric_test else 'sign_test'
        self.parametric_test = parametric_test
        self.plot_log_scale = plot_log_scale
        self.plot_style = plot_style
        self.n_bins = n_bins
        super().__init__(genes, attributes, alpha, attr_ref_path, save_csv, fname, return_fig, True, set_name, parallel,
                         enrichment_func_name, biotypes, background_set, biotype_ref_path, single_set=False)

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
        obs = self.annotation_df.loc[self.gene_set, attribute].median(skipna=False)
        if np.isnan(obs):
            pval = np.nan
        else:
            _, pval = sign_test(self.annotation_df.loc[self.gene_set, attribute].values, exp)
        return [attribute, len(self.gene_set), obs, exp, pval]

    def _one_sample_t_test_enrichment(self, attribute: str) -> list:
        exp = self.annotation_df[attribute].mean(skipna=True)
        obs = self.annotation_df.loc[self.gene_set, attribute].mean(skipna=False)
        if np.isnan(obs):
            pval = np.nan
        else:
            _, pval = ttest_1samp(self.annotation_df.loc[self.gene_set, attribute], popmean=exp)
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
        obs = exp.loc[self.gene_set]
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
                 'plot_go_network': 'indicates whether to plot GO network of the statistically significant GO Terms',
                 'attributes_set': 'set of the attributes/GO Terms for which enrichment should be calculated'}
    printout_params = "have any GO Annotations asocciated with them"
    GOA_DF_QUERIES = {}

    def __init__(self, genes: Union[set, np.ndarray], organism: Union[str, int], gene_id_type: str, alpha: float,
                 propagate_annotations: str, aspects: Union[str, Iterable[str]],
                 evidence_types: Union[str, Iterable[str]], excluded_evidence_types: Union[str, Iterable[str]],
                 databases: Union[str, Iterable[str]], excluded_databases: Union[str, Iterable[str]],
                 qualifiers: Union[str, Iterable[str]], excluded_qualifiers: Union[str, Iterable[str]],
                 return_nonsignificant: bool, save_csv: bool, fname: str, return_fig: bool, plot_horizontal: bool,
                 plot_go_network: bool, set_name: str, parallel: bool, enrichment_func_name: str, biotypes=None,
                 background_set: set = None, biotype_ref_path: str = None, single_set: bool = False,
                 random_seed: int = None, **pvalue_kwargs):
        self.dag_tree: ontology.DAGTree = io.fetch_go_basic()
        self.mod_annotation_dfs: Tuple[pd.DataFrame, ...] = tuple()
        self.organism = organism
        self.taxon_id = None
        self.gene_id_type = gene_id_type
        self.propagate_annotations = propagate_annotations.lower()
        self.aspects = aspects
        self.evidence_types = evidence_types
        self.excluded_evidence_types = excluded_evidence_types
        self.databases = databases
        self.excluded_databases = excluded_databases
        self.qualifiers = qualifiers
        self.excluded_qualifiers = excluded_qualifiers
        self.return_nonsignificant = return_nonsignificant
        self.plot_go_network = plot_go_network
        self.attributes_set: set = set()
        super().__init__(genes, [], alpha, '', save_csv, fname, return_fig, plot_horizontal, set_name, parallel,
                         enrichment_func_name, biotypes, background_set, biotype_ref_path, single_set, random_seed,
                         **pvalue_kwargs)

    def run(self):
        self.get_organism()
        return super().run()

    def _get_enrichment_func(self, pval_func_name: str):
        enrichment_func = super()._get_enrichment_func(pval_func_name)
        if self.propagate_annotations in {'weight'} and pval_func_name.lower() != 'fisher':
            raise NotImplementedError("The 'weight' propagation algorithm is only compatible with Fisher's Exact test.")
        elif self.propagate_annotations in {'allm'} and pval_func_name.lower() != 'fisher':
            warnings.warn(f"The 'weight' propagation algorithm is only compatible with Fisher's Exact test. "
                          f"Therefore, when calculating 'allm' p-values, the 'weight' method will use "
                          f"Fisher's Exact test, while the rest of the methods will use the '{pval_func_name}' method.")
        return enrichment_func

    def get_organism(self):
        if isinstance(self.organism, str) and self.organism.lower() == 'auto':
            self.taxon_id, self.organism = io.infer_taxon_from_gene_ids(self.gene_set)
        else:
            self.taxon_id, self.organism = io.map_taxon_id(self.organism)

    def fetch_annotations(self):
        # check if annotations for the requested query were previously fetched and cached
        query_key = self._get_query_key()
        if query_key in self.GOA_DF_QUERIES:
            self.annotation_df = self.GOA_DF_QUERIES[query_key]
            return
        else:
            self.annotation_df = self._generate_goa_df()
            # save query results to GOA_DF_QUERIES
            self.GOA_DF_QUERIES[query_key] = self.annotation_df

    def _generate_goa_df(self) -> pd.DataFrame:
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
        for annotation in tqdm(annotation_iter, desc=desc, total=annotation_iter.n_annotations, unit=' annotations'):
            # extract gene_id, go_id, source from the annotation
            gene_id: str = annotation['bioentity_internal_id']
            go_id = self.dag_tree[annotation['annotation_class']].id
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
        for source in source_to_gene_id_dict:
            translator = io.map_gene_ids(source_to_gene_id_dict[source], source, self.gene_id_type)
            for gene_id in sparse_annotation_dict.copy():
                if gene_id in translator:
                    translated_sparse_annotation_dict[translator[gene_id]] = sparse_annotation_dict.pop(gene_id)
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

    def plot_results(self) -> Union[plt.Figure, Tuple[plt.Figure, plt.Figure]]:
        n_bars = min(10, len(self.results))
        if self.single_set:
            title = f"Single-list GO enrichment for {self.set_name}" + f"\ntop {n_bars} most specific GO terms"
            bar_plot = self.enrichment_bar_plot(n_bars=n_bars, title=title, ylabel=r"$\log_2$(XL-mHG enrichment score)")
        else:
            title = f"GO enrichment for {self.set_name}" + f"\ntop {n_bars} most specific GO terms"
            bar_plot = self.enrichment_bar_plot(n_bars=n_bars, title=title)
        if self.plot_go_network:
            return bar_plot, self.go_dag_plot()
        return bar_plot

    def go_dag_plot(self):
        # TODO: implement me with pydot/graphviz!
        raise NotImplementedError

    def format_results(self, unformatted_results_dict: dict):
        if self.single_set:
            columns = ['name', 'samples', self.en_score_col, 'pval']
        else:
            columns = ['name', 'samples', 'obs', 'exp', self.en_score_col, 'pval']
        self.results = pd.DataFrame.from_dict(unformatted_results_dict, orient='index', columns=columns)
        self._correct_multiple_comparisons()
        # filter non-significant results
        if not self.return_nonsignificant:
            self.results = self.results[self.results['significant']]
        # sort results by specificity (level in DAG tree)
        self.results['dag_level'] = [self.dag_tree[ind].level for ind in self.results.index]
        self.results = self.results.sort_values('dag_level', ascending=False).drop('dag_level', 1).rename_axis('go_id')

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
        result_dicts = generic.ProgressParallel(n_jobs=-1, max_nbytes=max_nbytes, desc=progress_bar_desc)(
            joblib.delayed(func)(self, group, ind) for group, ind in zip(grouping, mod_df_inds))
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
                self.mod_annotation_dfs[mod_df_ind].loc[marked_nodes[go_id], go_id] = False
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
        observed_fraction = self.annotation_df.loc[self.gene_set, go_id].sum() / n
        mod_observed_fraction = bg_df.loc[self.gene_set].sum() / n
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_randomization_pval(n, log2_fold_enrichment, bg_df.values, reps, mod_observed_fraction)
        return [go_name, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval]

    def _xlmhg_enrichment(self, go_id: str, mod_df_ind: int = None) -> list:
        go_name = self.dag_tree[go_id].name
        res = super()._xlmhg_enrichment(go_id)
        res[0] = go_name
        return res

    def _xlmhg_index_vectors(self, attribute) -> np.ndarray:
        ranked_srs = self.annotation_df.loc[self.ranked_genes, attribute]
        assert ranked_srs.shape[0] == len(self.ranked_genes)
        return np.uint16(np.nonzero(ranked_srs.values)[0])

    def _hypergeometric_enrichment(self, go_id: str, mod_df_ind: int = None) -> list:
        bg_size, de_size, go_size, go_de_size = self._get_hypergeometric_parameters(go_id, mod_df_ind=mod_df_ind)

        go_name = self.dag_tree[go_id].name
        expected_fraction = self.annotation_df[go_id].sum() / bg_size
        observed_fraction = self.annotation_df.loc[self.gene_set, go_id].sum() / de_size
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_hypergeometric_pval(bg_size=bg_size, de_size=de_size, go_size=go_size, go_de_size=go_de_size)
        obs, exp = int(de_size * observed_fraction), de_size * expected_fraction
        return [go_name, de_size, obs, exp, log2_fold_enrichment, pval]

    def _fisher_enrichment(self, go_id: str, mod_df_ind: int = None) -> list:
        bg_size, de_size, go_size, go_de_size = self._get_hypergeometric_parameters(go_id, mod_df_ind)

        expected_fraction = self.annotation_df[go_id].sum() / bg_size
        observed_fraction = self.annotation_df.loc[self.gene_set, go_id].sum() / de_size
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = self._calc_fisher_pval(bg_size=bg_size, de_size=de_size, go_size=go_size, go_de_size=go_de_size)
        obs, exp = int(de_size * observed_fraction), de_size * expected_fraction
        return [self.dag_tree[go_id].name, de_size, obs, exp, log2_fold_enrichment, pval]

    def _get_hypergeometric_parameters(self, go_id: str, mod_df_ind: int = None) -> Tuple[int, int, int, int]:
        bg_size = self.mod_annotation_dfs[mod_df_ind].shape[0]
        de_size = len(self.gene_set)
        go_size = int(np.ceil(self.mod_annotation_dfs[mod_df_ind][go_id].sum()))
        go_de_size = int(np.ceil(self.mod_annotation_dfs[mod_df_ind].loc[self.gene_set, go_id].sum()))
        return bg_size, de_size, go_size, go_de_size
