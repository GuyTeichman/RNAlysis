import abc
import collections
import itertools
import logging
import warnings
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from typing import Iterable, List, Tuple, Union, Collection, Set, Dict, Literal

import joblib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import polars.selectors as cs
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


def dot_scaling_func(values: Union[float, np.ndarray, pl.Series]):
    return 8 * np.array(values) ** 0.75


class SizeHandler:
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


class StatsTest(abc.ABC):
    @abc.abstractmethod
    def run(self, attribute_name: str, annotation_set: Set[str], gene_set: Set[str], background_set: Set[str],
            annotation_set_unfiltered: Union[Set[str], None] = None) -> list:
        raise NotImplementedError

    @staticmethod
    def get_hypergeometric_params(annotation_set: Set[str], gene_set: Set[str], background_set: Set[str]
                                  ) -> Tuple[int, int, int, int, int, float, float]:
        bg_size = len(background_set)
        en_size = len(gene_set)
        attr_size = len(annotation_set)
        en_attr_size = len(gene_set.intersection(annotation_set))
        expected_fraction = attr_size / bg_size
        observed_fraction = en_attr_size / en_size
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        obs, exp = int(en_size * observed_fraction), en_size * expected_fraction

        return bg_size, en_size, attr_size, en_attr_size, obs, exp, log2_fold_enrichment


class TTest(StatsTest):
    def run(self, attribute_name: str, annotation_dict: Dict[str, float], gene_set: Set[str], background_set: Set[str],
            annotation_set_unfiltered=None):
        assert annotation_set_unfiltered is None
        obs_values = [annotation_dict[gene] for gene in gene_set]
        exp = np.nanmean([annotation_dict[gene] for gene in background_set])
        obs = np.mean(obs_values)
        if np.isnan(obs):
            pval = np.nan
        else:
            _, pval = ttest_1samp(obs_values, popmean=exp)
        return [attribute_name, len(gene_set), obs, exp, pval]


class SignTest(StatsTest):
    def run(self, attribute_name: str, annotation_dict: Dict[str, float], gene_set: Set[str], background_set: Set[str],
            annotation_set_unfiltered=None):
        assert annotation_set_unfiltered is None
        obs_values = [annotation_dict[gene] for gene in gene_set]
        exp = np.nanmedian([annotation_dict[gene] for gene in background_set])
        obs = np.median(obs_values)
        if np.isnan(obs):
            pval = np.nan
        else:
            _, pval = sign_test(obs_values, exp)
        return [attribute_name, len(gene_set), obs, exp, pval]


class PermutationTest(StatsTest):
    def __init__(self, n_permutations: int, random_seed: int = None):
        self.n_permutations = n_permutations
        self.random_seed = random_seed

    def run(self, attribute_name: str, annotation_set: Set[str], gene_set: Set[str], background_set: Set[str],
            annotation_set_unfiltered: Union[Set[str], None] = None):
        np.random.seed(self.random_seed)
        bg_size, en_size, attr_size, en_attr_size, obs, exp, log2fc = self.get_hypergeometric_params(
            annotation_set, gene_set, background_set)
        pval = self._calc_permutation_pval(log2fc, self.n_permutations, obs / en_size, bg_size, en_size, attr_size)
        if annotation_set_unfiltered is not None:
            _, _, _, _, obs, exp, log2fc = self.get_hypergeometric_params(annotation_set_unfiltered, gene_set,
                                                                          background_set)
        return [attribute_name, en_size, obs, exp, log2fc, pval]

    @staticmethod
    @generic.numba.jit(nopython=True)
    def _calc_permutation_pval(log2fc: float, reps: int, obs_frac: float, bg_size: int, en_size: int, attr_size: int
                               ) -> float:  # pragma: no cover
        bg_array = np.empty(bg_size, dtype=generic.numba.bool_)
        bg_array[0:attr_size] = True
        bg_array[attr_size:] = False
        ind_range = np.arange(bg_array.shape[0])
        success = 0
        if log2fc >= 0:
            for _ in range(reps):
                success += np.sum(bg_array[np.random.choice(ind_range, en_size, replace=False)]) / en_size >= obs_frac

        else:
            for _ in range(reps):
                success += np.sum(bg_array[np.random.choice(ind_range, en_size, replace=False)]) / en_size <= obs_frac
        pval = (success + 1) / (reps + 1)
        return pval


class FishersExactTest(StatsTest):
    def run(self, attribute_name: str, annotation_set: Set[str], gene_set: Set[str], background_set: Set[str],
            annotation_set_unfiltered: Union[Set[str], None] = None):
        bg_size, en_size, attr_size, en_attr_size, obs, exp, log2fc = self.get_hypergeometric_params(annotation_set,
                                                                                                     gene_set,
                                                                                                     background_set)
        pval = self._calc_fisher_pval(bg_size, en_size, attr_size, en_attr_size)
        if annotation_set_unfiltered is not None:
            _, _, _, _, obs, exp, log2fc = self.get_hypergeometric_params(annotation_set_unfiltered, gene_set,
                                                                          background_set)
        return [attribute_name, en_size, obs, exp, log2fc, pval]

    def run_weight(self, attribute_name: str, annotation_dict: Dict[str, float], gene_set: Set[str],
                   background_set: Set[str]):
        bg_size = len(background_set)
        en_size = len(gene_set)
        attr_size = sum(annotation_dict.values())
        en_attr_size = sum(val for key, val in annotation_dict.items() if key in gene_set)
        pval = self._calc_fisher_pval(bg_size, en_size, attr_size, en_attr_size)
        bg_size, en_size, attr_size, en_attr_size, obs, exp, log2fc = self.get_hypergeometric_params(annotation_dict,
                                                                                                     gene_set,
                                                                                                     background_set)

        return [attribute_name, en_size, obs, exp, log2fc, pval]

    @staticmethod
    @lru_cache(maxsize=256, typed=False)
    def _calc_fisher_pval(bg_size: int, en_size: int, attr_size: int, en_attr_size: int):
        contingency_table = [[en_attr_size, attr_size - en_attr_size],
                             [en_size - en_attr_size, bg_size - attr_size - en_size + en_attr_size]]
        _, pval = fisher_exact(contingency_table)
        return pval


class HypergeometricTest(StatsTest):
    def run(self, attribute_name: str, annotation_set: Set[str], gene_set: Set[str], background_set: Set[str],
            annotation_set_unfiltered: Union[Set[str], None] = None):
        bg_size, en_size, attr_size, en_attr_size, obs, exp, log2fc = self.get_hypergeometric_params(annotation_set,
                                                                                                     gene_set,
                                                                                                     background_set)

        pval = self._calc_hg_pval(bg_size, en_size, attr_size, en_attr_size)
        if annotation_set_unfiltered is not None:
            _, _, _, _, obs, exp, log2fc = self.get_hypergeometric_params(annotation_set_unfiltered, gene_set,
                                                                          background_set)
        return [attribute_name, en_size, obs, exp, log2fc, pval]

    @staticmethod
    def _calc_hg_pval(bg_size: int, en_size: int, attr_size: int, en_attr_size: int):
        try:
            if en_attr_size / en_size < attr_size / bg_size:
                return hypergeom.cdf(en_attr_size, bg_size, attr_size, en_size)
            return hypergeom.sf(en_attr_size - 1, bg_size, attr_size, en_size)
        except ZeroDivisionError:
            return hypergeom.cdf(en_attr_size, bg_size, attr_size, en_size)


class XlmhgTest(StatsTest):
    def __init__(self, min_positive_genes: int = 10, lowest_cutoff: float = 0.25):
        self.min_positive_genes = min_positive_genes
        self.lowest_cutoff = lowest_cutoff

    def run(self, attribute_name: str, annotation_set: Set[str], gene_set: List[str], background_set: Set[str] = None,
            annotation_set_unfiltered: Union[Set[str], None] = None):
        assert background_set is None
        if not HAS_XLMHG:
            warnings.warn("Package 'xlmhglite' is not installed. \n"
                          "If you want to run single-set enrichment analysis, "
                          "please install package 'xlmhglite' and try again. ")
            return [attribute_name, np.nan, np.nan, np.nan, np.nan, np.nan]

        index_vec, rev_index_vec = self._generate_xlmhg_index_vectors(gene_set, annotation_set)
        n, X, L, table = self._get_xlmhg_parameters(index_vec, gene_set)

        res_obj_fwd = xlmhglite.get_xlmhg_test_result(N=n, indices=index_vec, X=X, L=L, table=table)
        res_obj_rev = xlmhglite.get_xlmhg_test_result(N=n, indices=rev_index_vec, X=X, L=L, table=table)

        obs, exp, en_score, pval = self._extract_xlmhg_results(res_obj_fwd, res_obj_rev)
        return [attribute_name, n, obs, exp, en_score, pval]

    def _get_xlmhg_parameters(self, index_vec, ranked_genes):
        n = len(ranked_genes)
        # X = the minimal amount of 'positive' elements above the hypergeometric cutoffs out of all the positive
        # elements in the ranked set.
        X = min(self.min_positive_genes, n)
        # L = the lowest possible cutoff (n) to be tested out of the entire list.
        # Determined to be floor(lowest_cutoff * N), where 'N' is total number of elements in the ranked set (n).
        L = int(np.floor(self.lowest_cutoff * n))
        L = min(L, n)
        # pre-allocate empty array to speed up computation
        table = np.empty((len(index_vec) + 1, n - len(index_vec) + 1), dtype=np.longdouble)
        return n, X, L, table

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

    @staticmethod
    def _generate_xlmhg_index_vectors(ranked_genes, annotation_set) -> Tuple[np.ndarray, np.ndarray]:
        n = len(ranked_genes)
        index_vec = np.uint16([i for i in range(n) if ranked_genes[i] in annotation_set])
        rev_index_vec = np.uint16([n - 1 - index_vec[i - 1] for i in range(len(index_vec), 0, -1)])
        return index_vec, rev_index_vec


class EnrichmentPlotter(abc.ABC):
    __slots__ = ('results', 'en_score_col')

    def __init__(self, results: pl.DataFrame, en_score_col: str):
        assert len(results) > 0, "No enrichment results to plot."
        assert en_score_col in results.columns, f"Column '{en_score_col}' not found in the results DataFrame."
        self.results = results
        self.en_score_col = en_score_col

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

    @abc.abstractmethod
    def run(self):
        pass


class BarPlotter(EnrichmentPlotter):
    __slots__ = (
        'n_bars', 'center_bars', 'ylabel', 'title', 'ylim', 'title_fontsize', 'label_fontsize', 'ylabel_fontsize',
        'plot_style', 'plot_horizontal', 'alpha', 'show_expected')

    def __init__(self, results: pl.DataFrame, en_score_col: str, n_bars: Union[int, Literal['all']],
                 alpha: float, ylabel: str, title: str, ylim: Union[float, Literal['auto']],
                 title_fontsize: float, label_fontsize: float, ylabel_fontsize: float,
                 plot_style: Literal['bar', 'lollipop'], plot_horizontal: bool, show_expected: bool, center_bars: bool):
        super().__init__(results, en_score_col)

        assert plot_style in ['bar', 'lollipop'], \
            f"'plot_style' must be 'bar' or 'lollipop', instead got '{self.plot_style}'."

        if n_bars == 'all':
            n_bars = len(self.results)
        else:
            assert isinstance(n_bars, int) and n_bars > 0, f"Invalid value for 'n_bars': {n_bars}."

        self.n_bars = n_bars
        self.ylabel = ylabel
        self.title = title
        self.ylim = ylim
        self.alpha = alpha
        self.title_fontsize = title_fontsize
        self.label_fontsize = label_fontsize
        self.ylabel_fontsize = ylabel_fontsize
        self.plot_style = plot_style
        self.plot_horizontal = plot_horizontal
        self.center_bars = center_bars
        self.show_expected = show_expected

    def run(self):
        results = self.results.head(self.n_bars)

        if 'name' in results.columns:
            enrichment_names = results['name'].to_list()
        else:
            enrichment_names = results.select(pl.first()).to_series().to_list()
        enrichment_scores = results[self.en_score_col].to_list()
        enrichment_pvalue = results['padj'].to_list()
        enrichment_obs = results['obs'].to_list()
        enrichment_exp = results['exp'].to_list()

        # choose functions and parameters according to the graph's orientation (horizontal vs vertical)
        if self.plot_horizontal:
            figsize = [10.5, 0.4 * (4.8 + self.results.shape[0])]
            bar_func = plt.Axes.barh if self.plot_style == 'bar' else plt.Axes.hlines
            line_func = plt.Axes.axvline

            cbar_location = 'bottom'
            cbar_orientation = 'horizontal'
            tick_func = plt.Axes.set_yticks
            ticklabels_func = plt.Axes.set_yticklabels
            ticklabels_kwargs = dict(fontsize=self.label_fontsize, rotation=0)

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
            ticklabels_kwargs = dict(fontsize=self.label_fontsize, rotation=45)

        # set enrichment scores which are 'inf' or '-inf' to be the second highest/lowest enrichment score in the list
        scores_no_inf = [abs(score) for score in enrichment_scores if score != np.inf and score != -np.inf]
        if len(scores_no_inf) == 0:
            scores_no_inf.append(3)

        for i in range(len(enrichment_scores)):
            if enrichment_scores[i] == -np.inf:
                enrichment_scores[i] = -max(scores_no_inf)
            elif enrichment_scores[i] == np.inf:
                enrichment_scores[i] = max(scores_no_inf)

        if self.ylim == 'auto':
            if len(scores_no_inf) > 0:
                max_score = max(max(scores_no_inf), 3)
            else:
                max_score = 3
        else:
            max_score = self.ylim

        # get color values for bars
        data_color_norm = [0.5 * (1 + i / (np.ceil(max_score))) * 255 for i in enrichment_scores]
        data_color_norm_8bit = [int(i) if i != np.inf and i != -np.inf else np.sign(i) * max(np.abs(scores_no_inf)) for
                                i in data_color_norm]
        my_cmap = plt.colormaps.get_cmap('coolwarm')
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
                      handler_map={size: SizeHandler(size.size) for size in legend_sizes}, loc='lower right')

        # determine bounds, and enlarge the bound by a small margin (0.5%) so nothing gets cut out of the figure
        if self.ylim == 'auto':
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
        ax.set_title(self.title, fontsize=self.title_fontsize)
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
        if self.center_bars:
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
        cbar_label_kwargs = dict(label=self.ylabel, fontsize=self.ylabel_fontsize, labelpad=15)
        cbar = fig.colorbar(sm, ticks=range(int(bounds[0]), int(bounds[1]) + 1), cax=cax, orientation=cbar_orientation)
        cbar.set_label(**cbar_label_kwargs)
        cax.tick_params(labelsize=14, pad=6)

        if not self.plot_horizontal:
            cax.yaxis.set_ticks_position('left')
            cax.yaxis.set_label_position('left')

        plt.show()
        return fig


class HistogramPlotter(EnrichmentPlotter):
    __slots__ = ('annotations', 'attributes', 'alpha', 'set_name', 'gene_set', 'background_set', 'plot_log_scale',
                 'n_bins', 'plot_style', 'parametric_test')

    def __init__(self, results: pl.DataFrame, annotations: pl.DataFrame, attributes: List[str], en_score_col: str,
                 alpha: float, set_name: str, gene_set: set, background_set: str, plot_log_scale: bool, n_bins: int,
                 plot_style: str, parametric_test: bool):
        super().__init__(results, en_score_col)

        assert all([attr in annotations.columns for attr in
                    attributes]), f"Attributes {attributes} not found in annotations DataFrame."

        self.alpha = alpha
        self.annotations = annotations.select(pl.first(), pl.col(attributes))
        self.attributes = attributes
        self.set_name = set_name
        self.gene_set = gene_set
        self.background_set = background_set
        self.plot_log_scale = plot_log_scale
        self.n_bins = n_bins
        self.plot_style = plot_style
        self.parametric_test = parametric_test

    def _plot_histogram(self, attribute: str):
        # generate observed and expected Series, either linear or in log10 scale
        exp = [self.annotations[attribute][gene] for gene in self.background_set]
        obs = [self.annotations[attribute][gene] for gene in self.gene_set]
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
            y, x, _ = ax.hist([exp, obs], **kwargs, color=colors, label=['Expected', 'Observed'])
            max_y_val = np.max(y)
        elif self.plot_style.lower() == 'overlap':
            y, _, _ = ax.hist(exp, **kwargs, color=colors[0], label='Expected')
            y2, _, _ = ax.hist(obs, **kwargs, color=colors[1], label='Observed')
            max_y_val = np.max([np.max(y), np.max(y2)])
        else:
            raise NotImplementedError(f"Plot style '{self.plot_style}' not implemented.")

        # set either mean or median as the measure of centrality
        if self.parametric_test:
            x_exp, x_obs = np.nanmean(exp), np.mean(obs)
            label_exp, label_obs = 'Expected mean', 'Observed mean'
        else:
            x_exp, x_obs = np.nanmedian(exp), np.median(obs)
            label_exp, label_obs = 'Expected median', 'Observed median'

        # add lines for mean/median of observed and expected distributions
        ax.vlines(x_exp, ymin=0, ymax=max_y_val * 1.1, color='blue', linestyle='dashed', linewidth=2,
                  label=label_exp)
        ax.vlines(x_obs, ymin=0, ymax=max_y_val * 1.1, color='red', linestyle='dashed', linewidth=2,
                  label=label_obs)

        # add significance notation
        asterisks, fontweight = self._get_pval_asterisk(
            self.results.filter(pl.first() == attribute).select('padj').item(), self.alpha)
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

    def run(self):
        figs = []
        for attr in self.attributes:
            figs.append(self._plot_histogram(attr))
        return figs


class KEGGPlotter(EnrichmentPlotter):
    def __init__(self):
        pass

    def run(self):
        pass


class GODagPlotter(EnrichmentPlotter):
    def __init__(self):
        pass

    def run(self):
        pass


class MultiPlotter(EnrichmentPlotter):
    def __init__(self, results: pl.DataFrame, en_score_col: str, plotter_objs: List[EnrichmentPlotter]):
        super().__init__(results, en_score_col)
        self.plotter_objs = plotter_objs

    def run(self):
        figs = []
        for plotter in self.plotter_objs:
            figs.extend(parsing.data_to_list(plotter.run()))
        return figs


class EnrichmentRunner:
    ENRICHMENT_SCORE_YLABEL = "$\\log_2$(Fold Enrichment)"
    SINGLE_SET_ENRICHMENT_SCORE_YLABEL = r"$\log_2$(XL-mHG enrichment score)"

    __slots__ = {'results': 'DataFrame containing enrichment analysis results',
                 'annotations': 'dict containing all annotation data',
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
                 'stats_test': 'the function to be used to calculate enrichment p-values',
                 'single_set': 'indicates whether enrichment is calculated on a single set of genes '
                               '(without background) or on a set of target genes and a set of background genes',
                 'background_set': 'the background set of genes for enrichment analysis',
                 'en_score_col': 'name of the enrichment score column in the results DataFrame',
                 'ranked_genes': 'the set of genes/genomic features whose enrichment to calculate, '
                                 'pre-sorted and ranked by the user',
                 'plot_style': 'plot style',
                 'show_expected': 'show observed/expected values on plot',
                 'annotated_genes': 'set of annotated genes'}
    printout_params = "appear in the Attribute Reference Table"

    def __init__(self, genes: Union[set, np.ndarray], attributes: Union[Iterable, str, int],
                 alpha: param_typing.Fraction, attr_ref_path: str, return_nonsignificant: bool, save_csv: bool,
                 fname: str, return_fig: bool, plot_horizontal: bool, set_name: str,
                 parallel_backend: Literal[PARALLEL_BACKENDS], stats_test: StatsTest,
                 background_set: set = None, exclude_unannotated_genes: bool = True,
                 single_set: bool = False, plot_style: Literal['bar', 'lollipop'] = 'bar',
                 show_expected: bool = False):
        self.results: pl.DataFrame = pl.DataFrame()
        self.annotations: dict = dict()
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
        self.stats_test = stats_test
        if self.stats_test is None:
            return
        self.exclude_unannotated_genes = exclude_unannotated_genes
        self.annotated_genes = None
        self.single_set = single_set
        if self.single_set:
            assert background_set is None, "Enrichment in single_set mode does not accept a 'background_set' argument."
            assert isinstance(genes, np.ndarray), f"Invalid type for argument 'genes' in single_set mode: " \
                                                  f"expected np.ndarray, instead got '{type(genes)}'."

            assert isinstance(stats_test, XlmhgTest), \
                f"Invalid enrichment function for single_set mode: '{stats_test}'."

            self.background_set = None
            self.en_score_col = 'log2_enrichment_score'
            self.ranked_genes = genes
            self._update_ranked_genes()
        else:
            self.background_set = background_set
            self.en_score_col = 'log2_fold_enrichment'
            self.ranked_genes = None

    @classmethod
    def from_results_df(cls, results: pl.DataFrame, alpha: float, plot_horizontal: bool, set_name: str):
        runner = cls.__new__(cls)
        runner.results = results
        runner.alpha = alpha
        runner.plot_horizontal = plot_horizontal
        runner.set_name = set_name
        return runner

    def run(self, plot: bool = True) -> Union[pl.DataFrame, Tuple[pl.DataFrame, plt.Figure]]:
        if self.stats_test is None or not self.validate_statistical_test():
            return pl.DataFrame()
        self.fetch_annotations()
        self.get_attribute_names()
        if not self.single_set:
            self.get_background_set()
        self.update_gene_set()
        if len(self.gene_set) == 0:
            warnings.warn('After removing unannotated genes and/or genes not in the background set, '
                          'the enrichment set is empty. '
                          'Therefore, RNAlysis will not proceed with enrichment analysis. ')
            return pl.DataFrame()
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
        self.background_set = parsing.data_to_set(self.background_set)

        orig_bg_size = len(self.background_set)

        if self.exclude_unannotated_genes:
            if self.annotated_genes is None:
                annotated_genes = set()
                for s in self.annotations.values():
                    annotated_genes.update(s)
            else:
                annotated_genes = self.annotated_genes
            self.background_set = self.background_set.intersection(annotated_genes)

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

    def validate_statistical_test(self):
        if not isinstance(self.stats_test, StatsTest):
            warnings.warn(f"Invalid statistical test: {self.stats_test}.")
            return False
        if isinstance(self.stats_test, XlmhgTest) and not HAS_XLMHG:
            warnings.warn("Package 'xlmhglite' is not installed. \n"
                          "If you want to run single-set enrichment analysis, "
                          "please install package 'xlmhglite' and try again. ")
            return False
        return True

    def _get_background_set_from_set(self):
        self.background_set = parsing.data_to_set(self.background_set)

    def update_gene_set(self):
        if self.single_set:
            updated_gene_set = self.gene_set.intersection(set.union(*[i for i in self.annotations.values()]))
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
        annotation_df = io.load_table(self.attr_ref_path)
        validation.validate_attr_table(annotation_df)
        self.annotated_genes = parsing.data_to_set(annotation_df.select(pl.first()))
        self.annotations = {}
        for attr in annotation_df.columns[1:]:
            self.annotations[attr] = parsing.data_to_set(
                annotation_df.filter(pl.col(attr).is_not_null()).select(pl.first()))

    def get_attribute_names(self):
        all_attrs = parsing.data_to_list(self.annotations.keys())
        if self.attributes == 'all':
            self.attributes = all_attrs
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
        self.annotations = {attr: self.annotations[attr] for attr in self.attributes}
        if not self.single_set:
            for attr in self.annotations:
                if isinstance(self.annotations[attr], set):
                    self.annotations[attr].intersection_update(self.background_set)
                else:
                    self.annotations[attr] = {k: v for k, v in self.annotations[attr].items() if
                                              k in self.background_set}

    def calculate_enrichment(self) -> list:
        if self.parallel_backend != 'sequential' and len(self.attributes) > 5:
            return self._calculate_enrichment_parallel()
        return self._calculate_enrichment_serial()

    def _calculate_enrichment_serial(self) -> list:
        gene_set = self.ranked_genes if self.single_set else self.gene_set
        result = []
        for attribute in tqdm(self.attributes, desc="Calculating enrichment", unit='attributes'):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            this_res = self.stats_test.run(attribute, self.annotations[attribute], gene_set, self.background_set)
            result.append(this_res)
        return result

    def _calculate_enrichment_parallel(self) -> list:
        gene_set = self.ranked_genes if self.single_set else self.gene_set
        result = generic.ProgressParallel(desc="Calculating enrichment", unit='attribute', n_jobs=-1,
                                          backend=self.parallel_backend)(
            joblib.delayed(self.stats_test.run)(attribute, self.annotations[attribute], gene_set,
                                                self.background_set) for attribute in self.attributes)
        return result

    def format_results(self, unformatted_results_list: list):
        columns = ['name', 'samples', 'obs', 'exp', self.en_score_col, 'pval']
        self.results = pl.DataFrame(unformatted_results_list, schema=columns)
        self._correct_multiple_comparisons()

        # filter non-significant results
        if not self.return_nonsignificant:
            self.results = self.results.filter(pl.col('significant') == True)

    def _correct_multiple_comparisons(self):
        results_nulls = self.results.filter(pl.col('pval').is_nan())
        self.results = self.results.filter(pl.col('pval').is_not_nan())

        significant, padj = multitest.fdrcorrection(self.results['pval'].to_list(), alpha=self.alpha)
        self.results = pl.concat([self.results.with_columns(padj=padj, significant=significant), results_nulls],
                                 how='align')

    def plot_results(self) -> plt.Figure:
        if len(self.results) == 0:
            return plt.Figure()
        if self.single_set:
            return self.enrichment_bar_plot(ylabel=self.SINGLE_SET_ENRICHMENT_SCORE_YLABEL,
                                            title=f"Single-list enrichment for gene set '{self.set_name}'")
        return self.enrichment_bar_plot(ylabel=self.ENRICHMENT_SCORE_YLABEL,
                                        title=f"Enrichment for gene set '{self.set_name}'")

    def enrichment_bar_plot(self, n_bars: Union[int, Literal['all']] = 'all',
                            center_bars: bool = True,
                            ylabel: str = "$\\log_2$(Fold Enrichment)", title: str = 'Enrichment results',
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
        :type ylabel: str (default="$\\log_2$(Fold Enrichment)")
        :param center_bars: if True, center the bars around Y=0. Otherwise, ylim is determined by min/max values.
        :type center_bars: bool (default=True)
        :param ylim: set the Y-axis limits. If `ylim`='auto', determines the axis limits automatically based on the data. \
        If `ylim` is a number, set the Y-axis limits to [-ylim, ylim].
        :type ylim: float or 'auto' (default='auto')
        :param title_fontsize: font size for the plot title.
        :type title_fontsize: float (default=18)
        :param label_fontsize: font size for the attribute labels on the plot.
        :type label_fontsize: float (default=13)
        :param ylabel_fontsize: font size for the y-axis colorbar label.
        :type ylabel_fontsize: float (default=16)
        :return: Figure object containing the bar plot
        :rtype: matplotlib.figure.Figure instance
        """
        assert self.plot_style in ['bar', 'lollipop'], \
            f"'plot_style' must be 'bar' or 'lollipop', instead got '{self.plot_style}'."
        assert len(self.results) > 0, "No enrichment results to plot."
        # determine number of entries/bars to plot
        if n_bars != 'all':
            assert isinstance(n_bars, int) and n_bars >= 0, f"Invalid value for 'n_bars': {n_bars}."
            results = self.results.head(n_bars)
        else:
            results = self.results
        # pull names/scores/pvals out to avoid accidentally changing the results DataFrame in-place
        if 'name' in results.columns:
            enrichment_names = results['name'].to_list()
        else:
            enrichment_names = results.select(pl.first()).to_series().to_list()
        enrichment_scores = results[self.en_score_col].to_list()
        enrichment_pvalue = results['padj'].to_list()
        enrichment_obs = results['obs'].to_list()
        enrichment_exp = results['exp'].to_list()

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
        my_cmap = plt.colormaps.get_cmap('coolwarm')
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
                      handler_map={size: SizeHandler(size.size) for size in legend_sizes}, loc='lower right')

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

    def __init__(self, genes: set, attributes: Union[Iterable, str, int], alpha: param_typing.Fraction,
                 background_set: set, attr_ref_path: str, save_csv: bool, fname: str,
                 return_fig: bool, plot_log_scale: bool, plot_style: Literal['interleaved', 'overlap'],
                 n_bins: int, set_name: str, parallel_backend: Literal[PARALLEL_BACKENDS], parametric_test: bool):

        assert isinstance(plot_log_scale, bool), f"Invalid type for 'plot_log_scale': '{plot_log_scale}'."
        assert plot_style in {'interleaved', 'overlap'}, f"Invalid value for 'plot_style': '{plot_style}'."
        assert isinstance(n_bins,
                          int) and n_bins > 0, f"'n_bins' must be a positive integer. Instead got {type(n_bins)}."

        stats_test = TTest() if parametric_test else SignTest()
        super().__init__(genes, attributes, alpha, attr_ref_path, True, save_csv, fname, return_fig, True, set_name,
                         parallel_backend, stats_test, background_set, exclude_unannotated_genes=True, single_set=False)
        self.parametric_test = parametric_test
        self.plot_log_scale = plot_log_scale
        self.plot_style = plot_style
        self.n_bins = n_bins

    def fetch_annotations(self):
        annotation_df = io.load_table(self.attr_ref_path)
        validation.validate_attr_table(annotation_df)
        self.annotations = {}
        for attr in annotation_df.columns[1:]:
            self.annotations[attr] = {key: val[0] for key, val in
                                      annotation_df.select(cs.first() | cs.by_name(attr)).rows_by_key(
                                          annotation_df.columns[0]).items()}
            for gene_id, val in self.annotations[attr].items():
                if val is None:
                    self.annotations[attr][gene_id] = np.nan

    def format_results(self, unformatted_results_list: list):
        columns = ['name', 'samples', 'obs', 'exp', 'pval']
        self.results = pl.DataFrame(unformatted_results_list, schema=columns, orient='row')
        self._correct_multiple_comparisons()
        if self.results['pval'].is_nan().any():
            warnings.warn(f"One or more of the genes in the background set contained NaN values in "
                          f"{len(self.results['pval'].is_nan())} of the attributes. "
                          f"P-values and plots will not be generated for those attributes. ")
        if not self.return_nonsignificant:
            self.results = self.results.filter(pl.col('significant') == True)

    def plot_results(self) -> List[plt.Figure]:
        if len(self.results) == 0:
            return []
        figs = []
        for attribute in self.attributes:
            padj = self.results.filter(pl.first() == attribute).select('padj').item()
            if (padj is not None) and (not np.isnan(padj)):
                fig = self.enrichment_histogram(attribute)
                figs.append(fig)
        return figs

    def enrichment_histogram(self, attribute):
        # generate observed and expected Series, either linear or in log10 scale
        exp = [self.annotations[attribute][gene] for gene in self.background_set]
        obs = [self.annotations[attribute][gene] for gene in self.gene_set]
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
            y, x, _ = ax.hist([exp, obs], **kwargs, color=colors, label=['Expected', 'Observed'])
            max_y_val = np.max(y)
        elif self.plot_style.lower() == 'overlap':
            y, _, _ = ax.hist(exp, **kwargs, color=colors[0], label='Expected')
            y2, _, _ = ax.hist(obs, **kwargs, color=colors[1], label='Observed')
            max_y_val = np.max([np.max(y), np.max(y2)])
        else:
            raise NotImplementedError(f"Plot style '{self.plot_style}' not implemented.")

        # set either mean or median as the measure of centrality
        if self.parametric_test:
            x_exp, x_obs = np.nanmean(exp), np.mean(obs)
            label_exp, label_obs = 'Expected mean', 'Observed mean'
        else:
            x_exp, x_obs = np.nanmedian(exp), np.median(obs)
            label_exp, label_obs = 'Expected median', 'Observed median'

        # add lines for mean/median of observed and expected distributions
        ax.vlines(x_exp, ymin=0, ymax=max_y_val * 1.1, color='blue', linestyle='dashed', linewidth=2,
                  label=label_exp)
        ax.vlines(x_obs, ymin=0, ymax=max_y_val * 1.1, color='red', linestyle='dashed', linewidth=2,
                  label=label_obs)

        # add significance notation
        asterisks, fontweight = self._get_pval_asterisk(
            self.results.filter(pl.first() == attribute).select('padj').item(), self.alpha)
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
                 plot_pathway_graphs: bool, set_name: str, parallel_backend: Literal[PARALLEL_BACKENDS],
                 stats_test: StatsTest, background_set: set = None, exclude_unannotated_genes: bool = True,
                 single_set: bool = False, pathway_graphs_format: Literal[GRAPHVIZ_FORMATS] = 'none',
                 plot_style: Literal['bar', 'lollipop'] = 'bar', show_expected: bool = False):
        super().__init__(genes, [], alpha, '', return_nonsignificant, save_csv, fname, return_fig, plot_horizontal,
                         set_name, parallel_backend, stats_test, background_set, exclude_unannotated_genes, single_set,
                         plot_style, show_expected)
        if not self.stats_test:
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
        super().format_results(unformatted_results_list)
        self.results = self.results.rename({'name': 'KEGG ID'})
        self.results.insert_column(1, pl.Series(
            [self.pathway_names_dict[kegg_id] for kegg_id in self.results['KEGG ID']]).alias('name'))

    def _correct_multiple_comparisons(self):
        results_nulls = self.results.filter(pl.col('pval').is_nan())
        self.results = self.results.filter(pl.col('pval').is_not_nan())
        significant, padj = multitest.fdrcorrection(self.results['pval'].to_list(), alpha=self.alpha, method='negcorr')
        self.results = pl.concat([self.results.with_columns(padj=padj, significant=significant), results_nulls],
                                 how='align')

    def fetch_annotations(self):
        # check if annotations for the requested query were previously fetched and cached
        query_key = self._get_query_key()
        if query_key in self.KEGG_DF_QUERIES:
            self.annotations, self.pathway_names_dict = self.KEGG_DF_QUERIES[query_key]
            return
        else:
            self.annotations, self.pathway_names_dict = self._generate_annotation_dict()
            # save query results to KEGG_DF_QUERIES
            self.KEGG_DF_QUERIES[query_key] = self.annotations, self.pathway_names_dict

    def _generate_annotation_dict(self) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        # fetch and process KEGG annotations
        annotation_dict, pathway_name_dict = self._process_annotations()
        print(f"Found annotations for {len(annotation_dict)} genes.")
        # translate gene IDs
        translated_annotation_dict = self._translate_gene_ids(annotation_dict)
        return translated_annotation_dict, pathway_name_dict

    def _process_annotations(self) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        desc = f"Fetching KEGG annotations for organism '{self.organism}' (taxon ID:{self.taxon_id})"

        annotation_dict = {}
        pathway_name_dict = {}
        annotation_iter = self._get_annotation_iterator()
        assert annotation_iter.n_annotations > 0, "No KEGG annotations were found for the given parameters. " \
                                                  "Please try again with a different set of parameters. "
        for pathway_id, pathway_name, annotations in tqdm(annotation_iter, desc=desc,
                                                          total=annotation_iter.n_annotations, unit=' annotations'):
            pathway_name_dict[pathway_id] = pathway_name
            # add annotation to annotation dictionary
            if pathway_id not in annotation_dict:
                annotation_dict[pathway_id] = set()
            annotation_dict[pathway_id].update(annotations)
        return annotation_dict, pathway_name_dict

    def _translate_gene_ids(self, annotation_dict: dict):
        translated_dict = {}
        source = 'KEGG'

        translator = io.GeneIDTranslator(source, self.gene_id_type).run(
            set.union(*annotation_dict.values()))
        self.gene_id_translator = translator
        for kegg_id in annotation_dict:
            translated_dict[kegg_id] = set()
            for gene_id in annotation_dict[kegg_id]:
                if gene_id in translator:
                    translated_dict[kegg_id].add(translator[gene_id])
        return translated_dict

    def _get_query_key(self):
        return self.taxon_id, self.gene_id_type

    def get_attribute_names(self):
        self.attributes = parsing.data_to_list(self.annotations.keys())
        self.attributes_set = parsing.data_to_set(self.attributes)

    def plot_results(self) -> List[plt.Figure]:
        if len(self.results) == 0:
            return []
        figs = [super().plot_results()]
        if self.plot_pathway_graphs:
            for pathway in self.pathway_names_dict:
                if pathway in self.results.select(pl.first()) and self.results[pathway, 'significant']:
                    fig = self.pathway_plot(pathway)
                    if fig is None:
                        break
                    figs.append(fig)
        return figs

    def pathway_plot(self, pathway_id: str) -> plt.Figure:
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
                 'mutable_annotations': "Additional copies of 'annotations' which are "
                                        "actively mutated by p-value propagation algorithms",
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
                 stats_test: StatsTest, background_set: set = None, exclude_unannotated_genes: bool = True,
                 single_set: bool = False, ontology_graph_format: Literal[GRAPHVIZ_FORMATS] = 'none',
                 plot_style: Literal['bar', 'lollipop'] = 'bar', show_expected: bool = False):

        self.propagate_annotations = propagate_annotations.lower()
        super().__init__(genes, [], alpha, '', return_nonsignificant, save_csv, fname, return_fig, plot_horizontal,
                         set_name, parallel_backend, stats_test, background_set,
                         exclude_unannotated_genes, single_set, plot_style, show_expected)
        if not self.stats_test:
            return
        self.dag_tree: ontology.DAGTree = ontology.fetch_go_basic()
        self.mutable_annotations: Tuple[dict, ...] = tuple()
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

    def fetch_annotations(self):
        # check if annotations for the requested query were previously fetched and cached
        query_key = self._get_query_key()
        if query_key in self.GOA_DF_QUERIES:
            self.annotations = self.GOA_DF_QUERIES[query_key]
            return
        else:
            self.annotations = self._generate_annotation_dict()
            # save query results to GOA_DF_QUERIES
            self.GOA_DF_QUERIES[query_key] = self.annotations

    def _generate_annotation_dict(self) -> Dict[str, Set[str]]:
        # fetch and process GO annotations
        annotation_dict, source_to_gene_id_dict = self._process_annotations()
        print(f"Found annotations for {len(annotation_dict)} genes.")
        # translate gene IDs
        translated_annotation_dict = self._translate_gene_ids(annotation_dict, source_to_gene_id_dict)
        return translated_annotation_dict

    def _process_annotations(self) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
        if self.propagate_annotations != 'no':
            desc = f"Fetching and propagating GO annotations for organism '{self.organism}' (taxon ID:{self.taxon_id})"
        else:
            desc = f"Fetching GO annotations for organism '{self.organism}' (taxon ID:{self.taxon_id})"

        annotation_dict = {}
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
            if go_id not in annotation_dict:
                annotation_dict[go_id] = (set())
            annotation_dict[go_id].add(gene_id)

            # add gene id and source to source dict
            if source not in source_to_gene_id_dict:
                source_to_gene_id_dict[source] = set()
            source_to_gene_id_dict[source].add(gene_id)
            # propagate annotations
            self._propagate_annotation(gene_id, go_id, annotation_dict)
        return annotation_dict, source_to_gene_id_dict

    def _get_annotation_iterator(self):
        return io.GOlrAnnotationIterator(self.taxon_id, self.aspects,
                                         self.evidence_types, self.excluded_evidence_types,
                                         self.databases, self.excluded_databases,
                                         self.qualifiers, self.excluded_qualifiers)

    def _propagate_annotation(self, gene_id: str, go_id: str, annotation_dict: dict):
        if self.propagate_annotations != 'no':
            parents = set(self.dag_tree.upper_induced_graph_iter(go_id))
            for parent in parents:
                if parent not in annotation_dict:
                    annotation_dict[parent] = set()
                annotation_dict[parent].add(gene_id)

    def _translate_gene_ids(self, annotation_dict: dict, source_to_gene_id_dict: dict):
        # TODO: profile performance and maybe optimize this
        translated_dict = {}
        translators = []
        for source in source_to_gene_id_dict:
            try:
                translator = io.GeneIDTranslator(source, self.gene_id_type).run(source_to_gene_id_dict[source])
                translators.append(translator)
            except AssertionError as e:
                if 'not a valid Uniprot Dataset' in "".join(e.args):
                    warnings.warn(f"Failed to map gene IDs for {len(source_to_gene_id_dict[source])} annotations "
                                  f"from dataset '{source}'.")
                else:
                    raise e
        for go_id in annotation_dict:
            translated_dict[go_id] = set()
            for gene_id in annotation_dict[go_id]:
                for translator in translators:
                    if gene_id in translator:
                        translated_dict[go_id].add(translator[gene_id])
                        break
        return translated_dict

    def _get_query_key(self):
        return (self.taxon_id, self.gene_id_type, parsing.data_to_tuple(self.aspects, sort=True),
                parsing.data_to_tuple(self.evidence_types, sort=True),
                parsing.data_to_tuple(self.excluded_evidence_types, sort=True),
                parsing.data_to_tuple(self.databases, sort=True),
                parsing.data_to_tuple(self.excluded_databases, sort=True),
                parsing.data_to_tuple(self.qualifiers, sort=True),
                parsing.data_to_tuple(self.excluded_qualifiers, sort=True),
                self.propagate_annotations != 'no')

    def validate_statistical_test(self):
        super().validate_statistical_test()
        if self.propagate_annotations in {'weight'} and not isinstance(self.stats_test, FishersExactTest):
            warnings.warn("The 'weight' propagation algorithm is only compatible with Fisher's Exact test.")
            return False
        if self.propagate_annotations in {'allm'} and not isinstance(self.stats_test, FishersExactTest):
            warnings.warn(f"The 'weight' propagation algorithm is only compatible with Fisher's Exact test. ")
            return False
        return True

    def get_attribute_names(self):
        self.attributes = parsing.data_to_list(self.annotations.keys())
        self.attributes_set = parsing.data_to_set(self.attributes)

    def _correct_multiple_comparisons(self):
        results_nulls = self.results.filter(pl.col('pval').is_nan())
        self.results = self.results.filter(pl.col('pval').is_not_nan())
        significant, padj = multitest.fdrcorrection(self.results['pval'].to_list(), alpha=self.alpha, method='negcorr')
        self.results = pl.concat([self.results.with_columns(padj=padj, significant=significant), results_nulls],
                                 how='align')

    def plot_results(self) -> List[plt.Figure]:
        if len(self.results) == 0:
            return []
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

    def go_dag_plot(self, title, dpi: int = 300, ylabel: str = "$\\log_2$(Fold Enrichment)") -> List[plt.Figure]:
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
        super().format_results(parsing.data_to_list(unformatted_results_dict.values()))
        self.results = self.results.rename({'name': 'GO ID'})
        self.results.insert_column(1, pl.Series([self.dag_tree[go_id].name for go_id in self.results['GO ID']]).alias(
            'name'))
        # sort results by specificity (level in DAG tree)
        self.results = self.results.with_columns(pl.Series(
            [self.dag_tree[ind].level for ind in self.results.select(pl.first()).to_series().to_list()]).alias(
            'dag_level'))
        self.results = self.results.sort(pl.col('dag_level'), descending=True).drop('dag_level')

    def _calculate_enrichment_serial(self) -> dict:
        desc = f"Calculating enrichment for {len(self.attributes)} GO terms " \
               f"using the '{self.propagate_annotations}' method"
        if self.propagate_annotations == 'classic' or self.propagate_annotations == 'no':
            self.mutable_annotations = (self.annotations,)
            result = self._go_classic_pvalues_serial(desc)
        elif self.propagate_annotations == 'elim':
            self.mutable_annotations = (deepcopy(self.annotations),)
            result = self._go_elim_pvalues_serial(desc)
        elif self.propagate_annotations == 'weight':
            self.mutable_annotations = ({attr: {v: 1.0 for v in val} for attr, val in self.annotations.items()},)
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
            self.mutable_annotations = (self.annotations,)
            result = self._go_classic_pvalues_parallel(desc)
        elif self.propagate_annotations == 'elim':
            self.mutable_annotations = tuple(
                {attr: deepcopy(self.annotations[attr]) for attr in self._go_level_iterator(namespace) if
                 attr in self.annotations} for namespace in
                self.dag_tree.namespaces)
            result = self._go_elim_pvalues_parallel(desc)
        elif self.propagate_annotations == 'weight':
            self.mutable_annotations = tuple({attr: {v: 1.0 for v in val} for attr, val in self.annotations.items() if
                                              attr in set(self._go_level_iterator(namespace))} for namespace in
                                             self.dag_tree.namespaces)
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

    def _go_classic_on_batch(self, go_term_batch: Iterable[str], annotation_idx: int = 0):
        """
        calculate stats and p-values with the 'classic' propagation algorithm for an iterable batch of GO terms.

        :param go_term_batch: batch of GO terms to calculate enrichment for
        :type go_term_batch: an iterable of str (GO IDs)
        :param annotation_idx:
        :type annotation_idx: int between 0 and 2 (default 0)
        :return:
        :rtype:
        """
        gene_set = self.ranked_genes if self.single_set else self.gene_set
        annotations = self.mutable_annotations[annotation_idx]
        return {go_id: self.stats_test.run(go_id, annotations[go_id], gene_set, self.background_set) for go_id in
                go_term_batch}

    def _go_elim_pvalues_serial(self, progress_bar_desc: str = '') -> dict:
        result = self._go_elim_on_aspect('all', progress_bar_desc=progress_bar_desc)
        return result

    def _parallel_over_grouping(self, func, grouping: Iterable, annotation_indices: Iterable[int],
                                max_nbytes: Union[str, None] = '1M', progress_bar_desc: str = '') -> dict:
        assert validation.is_method_of_class(func, type(self))
        result_dicts = generic.ProgressParallel(desc=progress_bar_desc, n_jobs=-1, max_nbytes=max_nbytes,
                                                backend=self.parallel_backend)(
            joblib.delayed(func)(group, ind) for group, ind in zip(grouping, annotation_indices))
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

    def _go_elim_on_aspect(self, go_aspect: str, annotation_idx: int = 0, progress_bar_desc: str = None) -> dict:
        result_dict = {}
        marked_nodes = {}
        annotations = self.mutable_annotations[annotation_idx]
        gene_set = self.ranked_genes if self.single_set else self.gene_set
        go_id_iter = self._go_level_iterator(go_aspect) if progress_bar_desc is None \
            else tqdm(self._go_level_iterator(go_aspect), unit=' GO terms', desc=progress_bar_desc,
                      total=len(self.attributes))
        for go_id in go_id_iter:
            if go_id in marked_nodes:  # if this node was marked, remove from it all marked genes
                annotations[go_id].difference_update(marked_nodes[go_id])
            result_dict[go_id] = self.stats_test.run(go_id, annotations[go_id], gene_set, self.background_set,
                                                     self.annotations[go_id])
            # if current GO ID is significantly ENRICHED, mark its ancestors
            if result_dict[go_id][-1] <= self.alpha and result_dict[go_id][-2] > 0:
                new_marked_genes = annotations[go_id].intersection(annotations[go_id])
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

    def _go_weight_on_aspect(self, go_aspect: str, annotation_idx: int = 0, progress_bar_desc: str = None) -> dict:
        result_dict = {}
        mut_annotation = self.mutable_annotations[annotation_idx]
        weights = collections.defaultdict(lambda: 1)  # default weight for all nodes is 1
        go_id_iter = self._go_level_iterator(go_aspect) if progress_bar_desc is None \
            else tqdm(self._go_level_iterator(go_aspect), unit=' GO terms', desc=progress_bar_desc,
                      total=len(self.attributes))
        for go_id in go_id_iter:
            children = {child for child in self.dag_tree[go_id].get_children() if child in self.attributes_set}
            self._compute_term_sig(mut_annotation, self.gene_set, self.background_set, self.dag_tree,
                                   self.attributes_set, go_id, children, weights, result_dict)
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

    @staticmethod
    def _compute_term_sig(mut_annotation: Dict[str, Dict[str, int]], gene_set: Set[str], background_set: Set[str],
                          dag_tree: ontology.DAGTree, attributes_set: Set[str],
                          go_id: str, children: set, weights: dict, result: dict, tolerance: float = 10 ** -50):
        """
        Implementation of the computeTermSig(u, children) function from Alexa et al 2006: \
        https://doi.org/10.1093/bioinformatics/btl140

        :param annotation_idx: index of the annotation dataframe to be used for the algorithm \
        (in case the algorithm runs in parallel and multiple dataframes are used)
        :type annotation_idx: int
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
        test = FishersExactTest()
        this_res = test.run_weight(go_id, mut_annotation[go_id], gene_set, background_set)
        if go_id in result:
            result[go_id][-1] = this_res[-1]
        else:
            result[go_id] = this_res
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
                for gene_id in mut_annotation[child]:
                    mut_annotation[child][gene_id] *= weights[child]
                result[child][-1] = test.run_weight(child, mut_annotation[child], gene_set, background_set)[-1]
            return

        # CASE 2: if some children are more significant than parent, re-weigh ancesctors (including 'go_id'),
        # and then recompute stats for go_id
        inclusive_ancestors = {ancestor for ancestor in
                               itertools.chain([go_id], dag_tree.upper_induced_graph_iter(go_id)) if
                               ancestor in attributes_set}
        for sig_child in sig_children:
            for inclusive_ancestor in inclusive_ancestors:
                for gene_id in mut_annotation[inclusive_ancestor]:
                    mut_annotation[inclusive_ancestor][gene_id] *= (1 / weights[sig_child])
        # re-run compute_term_sig, only with the children which were not more significant than their parents
        GOEnrichmentRunner._compute_term_sig(mut_annotation, gene_set, background_set, dag_tree, attributes_set, go_id,
                                             children.difference(sig_children), weights, result, tolerance)
