from abc import ABC

import itertools
import numpy as np
from sklearn_extra.cluster import KMedoids
from hdbscan import HDBSCAN
import matplotlib.pyplot as plt
import numba
import pairwisedist as pwdist
import pandas as pd
import seaborn as sns
from grid_strategy import strategies
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances, silhouette_score
from typing import List, Set, Dict, Tuple, Union, Iterable, Generic, Callable

from rnalysis.utils import generic, parsing, validation


class KMedoidsIter:
    __slots__ = {'n_clusters': 'number of clusters to find', 'metric': 'distance metric',
                 'n_init': 'number of initializations to run', 'init': 'initialization algorithm',
                 'max_iter': 'max number of algorithm iterations per initialization',
                 'random_state': 'random state of the random number generator',
                 'clusterer': 'the KMedoids clusterer object', 'inertia_': "the clustering solution's inertia",
                 'cluster_centers_': "the clustering solution's cluster centers",
                 'medoid_indices_': "the clustering solution's medoid indices",
                 'labels_': "the clustering solution's point labels"}

    def __init__(self, n_clusters: int, metric: str = 'euclidean', init: str = 'k-medoids++', max_iter: int = 300,
                 n_init: int = 10, random_state: int = None):
        assert isinstance(n_init, int), f"'n_init' must be an integer, is {type(n_init)} instead."
        assert isinstance(metric, str), f"'metric' must be a string, is {type(metric)} instead."
        assert n_init > 0, f"'n_init' must be a positive integer. Input {n_init} is invalid. "
        self.n_clusters = n_clusters
        self.metric = metric
        self.n_init = n_init
        self.init = init
        self.max_iter = max_iter
        self.random_state = random_state
        self.clusterer = KMedoids(n_clusters=self.n_clusters, metric=self.metric, init=self.init,
                                  max_iter=self.max_iter, random_state=random_state)
        self.inertia_ = None
        self.cluster_centers_ = None
        self.medoid_indices_ = None
        self.labels_ = None

    def __repr__(self):
        return repr(self.clusterer)

    def __str__(self):
        return str(self.clusterer)

    def fit(self, x):
        inertias = np.zeros(self.n_init)
        clusterers = []
        for i in range(self.n_init):
            if self.random_state is not None:
                clusterers.append(
                    KMedoids(n_clusters=self.n_clusters, metric=self.metric, init=self.init, max_iter=self.max_iter,
                             random_state=self.random_state + i).fit(x))
            else:
                clusterers.append(KMedoids(n_clusters=self.n_clusters, metric=self.metric, init=self.init,
                                           max_iter=self.max_iter).fit(x))
            inertias[i] = clusterers[i].inertia_
        best_clusterer = clusterers[int(np.argmax(inertias))]
        self.clusterer = best_clusterer
        self.inertia_ = self.clusterer.inertia_
        self.cluster_centers_ = self.clusterer.cluster_centers_
        self.medoid_indices_ = self.clusterer.medoid_indices_
        self.labels_ = self.clusterer.labels_
        return self

    def predict(self, x):
        return self.clusterer.predict(x)

    def fit_predict(self, x):
        self.fit(x)
        return self.predict(x)


class ClusteringRunner:
    precomputed_metrics = {'spearman': pwdist.spearman_distance, 'pearson': pwdist.pearson_distance,
                           'ys1': pwdist.ys1_distance, 'yr1': pwdist.yr1_distance,
                           'jackknife': pwdist.jackknife_distance}

    def __init__(self, data: pd.DataFrame, power_transform: bool, metric: Tuple[str, str] = None,
                 plot_style: str = 'none', split_plots: bool = False):
        assert plot_style.lower() in {'all', 'std_bar', 'std_area', 'none'}, \
            f"Invalid value for 'plot_style': '{plot_style}'. " \
            f"'plot_style' must be 'all', 'std_bar', 'std_area' or 'none'. "
        assert isinstance(split_plots, bool), f"'split_plots' must be of type bool, instead got {type(split_plots)}"

        self.clusterer_class: type
        self.clusterer_kwargs: dict
        self.clusterers: list = []

        self.data: pd.DataFrame = data
        self.power_transform: bool = power_transform
        self.transform: Callable = generic.standard_box_cox if power_transform else generic.standardize

        if metric is not None:
            metric_arg_name, metric_value = metric
            metric_value = metric_value.lower()
            validation.validate_clustering_parameters(metric_value)
            if metric_value in self.precomputed_metrics:
                def transform(x):
                    return self.precomputed_metrics[metric_value](self.transform(x))

                self.transform = transform
                self.metric = 'precomputed'
            else:
                self.metric = metric_value.lower()
            self.clusterer_kwargs[metric_arg_name] = self.metric

        self.plot_style: str = plot_style.lower()
        self.split_plots: bool = split_plots

    def run(self):
        raise NotImplementedError

    def plot_clustering(self, n_clusters: int, data: pd.DataFrame, labels: np.ndarray, centers: np.ndarray, title: str
                        ) -> Union[Tuple[plt.Figure, List[plt.Axes]], None]:
        if self.plot_style == 'none':
            return

        if self.split_plots:
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
        color_generator = generic.color_generator()
        for i, subplot in enumerate(subplots):
            if self.split_plots:
                axes.append(subplot[1])
                fig = subplot[0]
                fig.suptitle(title, fontsize=18)
            else:
                # noinspection PyUnboundLocalVariable
                axes.append(fig.add_subplot(subplot))
            mean = centers[i, :]
            stdev = data.loc[labels == i, :].T.std(axis=1, ddof=1)
            x = np.arange(data.shape[1]) + 0.5
            color = next(color_generator)
            max_y = max(max_y, np.max(stdev + mean))
            min_y = min(min_y, np.min(mean - stdev))
            axes[-1].axhline(color='grey', linewidth=2, linestyle='--')
            if self.plot_style.lower() == 'std_area':
                axes[-1].plot(x, mean, marker='o', linewidth=2, color=color, markeredgecolor='black')
                axes[-1].fill_between(x, mean + stdev, mean - stdev, alpha=0.2, color=color)
            else:
                axes[-1].errorbar(x, mean, marker='o', linewidth=2, color=color, yerr=stdev, markeredgecolor='black')
                if self.plot_style.lower() == 'all':
                    vals = data.loc[labels == i, :].T
                    axes[-1].plot(x, vals, color=color, alpha=0.05, linewidth=0.35)

            axes[-1].set_title(f"Cluster number {i + 1} ({np.count_nonzero(labels == i)} genes)")
            axes[-1].set_ylabel('Standardized expression', fontsize=ylabel_fontsize)
            axes[-1].set_xticks(x)
            axes[-1].set_xticklabels(data.columns, fontsize=ticklabels_fontsize)
        for ax in axes:
            ax.set_ylim((min_y, max_y))
        if not self.split_plots:
            fig.tight_layout(rect=[0, 0.03, 1, 0.92])
        plt.show()
        return fig, axes

    def clustering_labels_to_binary_format(self) -> np.ndarray:
        pass


class ClusteringRunnerWithNClusters(ClusteringRunner, ABC):
    def __init__(self, data: pd.DataFrame, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'default', plot_style: str = 'none',
                 split_plots: bool = False, metric: Tuple[str, str] = None):
        super().__init__(data, power_transform, metric, plot_style, split_plots)
        self.n_clusters: list = self.parse_n_clusters(n_clusters)
        self.max_n_clusters_estimate = min(20, data.shape[
            0] // 4) if max_n_clusters_estimate == 'default' else max_n_clusters_estimate

    def parse_n_clusters(self, n_clusters) -> List[int]:
        # get best n_clusters value using gap statistic/silhouette method, if requested by the user
        max_clusters = min(20, self.data.shape[0] // 4)  # if max_clusters == 'default' else max_clusters
        if isinstance(n_clusters, str) and n_clusters.lower() == 'silhouette':
            n_clusters = self.silhouette()
        elif isinstance(n_clusters, str) and n_clusters.lower() in {'gap', 'gap statistic', 'gap_statistic'}:
            n_clusters = self.gap_statistic()

        n_clusters = parsing.data_to_list(n_clusters)
        # make sure all values of n_clusters are in valid range
        n_clusters, n_clusters_copy = itertools.tee(n_clusters)
        assert np.all([isinstance(item, int) and 2 <= item <= len(self.data.index) for item in n_clusters_copy]), \
            f"Invalid value for n_clusters: '{n_clusters}'. n_clusters must be 'gap', 'silhouette', " \
            f"or an integer/Iterable of integers in range 2 <= n_clusters <= n_features."
        return parsing.data_to_list(n_clusters)

    def gap_statistic(self, random_state: int = None, n_refs: int = 10) -> int:
        raw_data = self.data.loc[:, self._numeric_columns].values
        # determine the range of k values to be tested
        n_clusters_range = np.arange(1, self.max_n_clusters_estimate + 1)
        print(f"Estimating the optimal number of clusters using the Gap Statistic method in range "
              f"{2}:{self.max_n_clusters_estimate}...")
        # set the random state, if one was supplied
        if random_state is not None:
            np.random.seed(random_state)
        # calculate the SVD of the observed data (X), and transform via X' = x_tag = dot(X, V)
        # note: the data is centered by sklearn.decomposition.PCA, no need to pre-center it.
        pca = PCA(random_state=random_state).fit(raw_data)
        x_tag = pca.transform(raw_data)
        # determine the ranges of the columns of X' (x_tag)
        a, b = x_tag.min(axis=0, keepdims=True), x_tag.max(axis=0, keepdims=True)
        # transform the observed data using Box-Cox, and then standardize it
        data = self.transform(raw_data)
        # generate 'n_refs' random reference arrays:
        # draw uniform features Z' over the ranges of the columns of X', and back-transform via Z = dot(Z', V.T), then
        # transform the random reference data using Box-Cox, and then standardize it
        refs = [self.transform(pca.inverse_transform(
            np.random.random_sample(size=raw_data.shape) * (b - a) + a))
            for _ in range(n_refs)]
        # allocate empty arrays for observed/expected log(inertia), gap scores Gap(K) and gap error S(K)
        log_disp_obs = np.zeros((len(n_clusters_range)))
        log_disp_exp = np.zeros((len(n_clusters_range)))
        gap_scores = np.zeros((len(n_clusters_range)))
        gap_error = np.zeros((len(n_clusters_range)))

        # iterate over K values in the range
        for ind, this_n_clusters in enumerate(n_clusters_range):
            # init the clusterer with given arguments
            clusterer = self.clusterer_class(n_clusters=this_n_clusters, **self.clusterer_kwargs)
            # cluster each of the n_refs reference arrays into Ki clusters
            ref_disps = np.zeros(n_refs)
            for ref_ind, ref in enumerate(refs):
                # calculate dispersion/inertia for reference data clustered into Ki clusters
                ref_disps[ref_ind] = self._compute_dispersion(clusterer.fit(ref), ref, this_n_clusters)
            # calculate dispersion/inertia for the observed data clustered into Ki clusters
            disp = self._compute_dispersion(clusterer.fit(data), data, this_n_clusters)
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
        k_candidates = n_clusters_range[:-1][diff >= 0]
        # if no K fulfilled the condition (the gap function is monotonically increasing),
        # return the last value of K as the best K
        best_n_clusters = k_candidates[0] if len(k_candidates) > 0 else n_clusters_range[-1]
        n_clusters_ind = np.argmax(n_clusters_range == best_n_clusters)
        print(f"Using the Gap Statistic method, {best_n_clusters} was chosen as the best number of clusters (K).")
        if len(k_candidates) >= 2:
            print(f"Other potentially good values of K that were found: {list(k_candidates[1:])}. ")

        # plot results
        fig = self._plot_gap_statistic(n_clusters_range, log_disp_obs, log_disp_exp, gap_scores, gap_error,
                                       best_n_clusters, n_clusters_ind)

        return best_n_clusters

    @staticmethod
    def _compute_dispersion(clusterer, data, n_clusters: int):
        # known as Wk in the Gap Statistic paper
        if hasattr(clusterer, 'inertia_'):
            return clusterer.inertia_
        labels = clusterer.labels_
        dispersion = 0
        for r in range(n_clusters):
            this_label = labels == r
            n = np.sum(this_label)
            # don't calculate within-cluster dispersion for empty clusters or clusters with a single member
            if n > 1:
                dispersion += np.sum(pairwise_distances(data[this_label, :]) ** 2) / (2 * n)
        return dispersion

    @staticmethod
    def _plot_gap_statistic(n_clusters_range, log_disp_obs, log_disp_exp, gap_scores, gap_error, n_clusters: int,
                            n_clusters_ind
                            ) -> plt.Figure:
        fig, (ax_inertia, ax) = plt.subplots(1, 2, figsize=(14, 9))
        ax_inertia.plot(n_clusters_range, log_disp_obs, '-o')
        ax_inertia.plot(n_clusters_range, log_disp_exp, '-o')
        ax_inertia.legend(['Observed', 'Expected'])
        ax_inertia.set_ylabel(r"$\ln$(dispersion)", fontsize=15)
        ax_inertia.set_xlabel("Number of clusters (K)", fontsize=15)
        ax_inertia.set_xticks(n_clusters_range)
        ax.errorbar(n_clusters_range, gap_scores, yerr=gap_error, marker='o', color='r')
        ax.set_title("Gap Statistic method for estimating optimal number of clusters")
        ax.set_ylabel('Gap Value', fontsize=15)
        ax.set_xlabel("Number of clusters (K)", fontsize=15)
        ax.annotate(f'Best K={n_clusters}',
                    xy=(n_clusters, (gap_scores[n_clusters_ind] - gap_error[n_clusters_ind]) / 1.05),
                    xytext=(n_clusters, (gap_scores[n_clusters_ind] - gap_error[n_clusters_ind]) / 1.2),
                    arrowprops=dict(facecolor='black'))
        ax.set_xticks(n_clusters_range)
        sns.despine()
        plt.show()
        return fig

    def silhouette(self) -> int:
        print(f"Estimating the optimal number of clusters using the Silhouette method in range "
              f"{2}:{self.max_n_clusters_estimate}...")
        data = self.transform(self.data.loc[:, self._numeric_columns].values)
        sil_scores = []
        k_range = np.arange(2, self.max_n_clusters_estimate + 1)
        for n_clusters in k_range:
            clusterer = self.clusterer_class(n_clusters=n_clusters, **self.clusterer_kwargs)
            sil_scores.append(silhouette_score(data, clusterer.fit_predict(data)))
        fig = plt.figure(figsize=(7, 9))
        ax = fig.add_subplot(111)
        ax.plot(k_range, sil_scores, '--o', color='r')
        ax.set_ylim(-1.0, 1.0)
        ax.set_title("Silhouette method for estimating optimal number of clusters")
        ax.set_ylabel('Average silhouette score')
        ax.set_xlabel("Number of clusters (k)")
        ax.axhline(0, color='grey', linewidth=2, linestyle='--')
        best_n_clusters = k_range[int(np.argmax(sil_scores))]
        ax.annotate(f'Best n_clusters={best_n_clusters}', xy=(best_n_clusters, max(sil_scores)),
                    xytext=(best_n_clusters + 1, 0.75),
                    arrowprops=dict(facecolor='black', shrink=0.15))
        ax.set_xticks(k_range)
        print(f"Using the Silhouette method, {best_n_clusters} was chosen as the best number of clusters (k).")
        return best_n_clusters


class KMeansRunner(ClusteringRunnerWithNClusters):
    clusterer_class = KMeans

    def __init__(self, data, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'default', random_state: int = None, n_init: int = 3,
                 max_iter: int = 300, plot_style: str = 'none', split_plots: bool = False):
        self.random_state = random_state
        self.n_init = n_init
        self.max_iter = max_iter
        self.clusterer_kwargs = dict(n_init=self.n_init, max_iter=self.max_iter, random_state=self.random_state)

        super(KMeansRunner, self).__init__(data, power_transform, n_clusters, max_n_clusters_estimate, plot_style,
                                           split_plots)

    def run(self) -> List[KMeans]:
        self.clusterers = []
        # generate standardized data for plots
        data_for_plot = generic.standardize(self.data)
        # iterate over all K values, generating clustering results and plotting them
        for this_n_clusters in self.n_clusters:
            this_n_clusters = int(this_n_clusters)
            clusterer = self.clusterer_class(n_clusters=this_n_clusters, **self.clusterer_kwargs).fit(
                self.transform(self.data.values))
            # plot results
            centers = clusterer.cluster_centers_
            self.plot_clustering(this_n_clusters, data_for_plot, clusterer.labels_, centers,
                                 title=f"Results of K-Means Clustering for n_clusters={this_n_clusters} and "
                                       f"power_transform={self.power_transform}")
            self.clusterers.append(clusterer)
        return self.clusterers


class KMedoidsRunner(ClusteringRunnerWithNClusters):
    clusterer_class = KMedoidsIter

    def __init__(self, data, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'default', metric: str = 'euclidean',
                 random_state: int = None, n_init: int = 3, max_iter: int = 300, plot_style: str = 'none',
                 split_plots: bool = False):
        self.random_state = random_state
        self.n_init = n_init
        self.max_iter = max_iter
        self.clusterer_kwargs = dict(n_init=self.n_init, max_iter=self.max_iter, random_state=self.random_state)
        super(KMedoidsRunner, self).__init__(data, power_transform, n_clusters, max_n_clusters_estimate, plot_style,
                                             split_plots, ('metric', metric))

    def run(self) -> List[KMedoidsIter]:
        self.clusterers = []
        # generate standardized data for plots
        data_for_plot = generic.standardize(self.data)
        # iterate over all K values, generating clustering results and plotting them
        for this_n_clusters in self.n_clusters:
            this_n_clusters = int(this_n_clusters)
            clusterer = self.clusterer_class(n_clusters=this_n_clusters, **self.clusterer_kwargs).fit(
                self.transform(self.data.values))
            # get cluster centers
            centers = data_for_plot.values[clusterer.medoid_indices_, :]
            # plot results
            self.plot_clustering(this_n_clusters, data_for_plot, clusterer.labels_, centers,
                                 title=f"Results of K-Medoids Clustering for n_clusters={this_n_clusters}, "
                                       f"metric='{self.metric}', power_transform={self.power_transform}")
            self.clusterers.append(clusterer)
        return self.clusterers


class HierarchicalRunner(ClusteringRunnerWithNClusters):
    clusterer_class = AgglomerativeClustering

    def __init__(self, data, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'default', metric: str = 'euclidean',
                 linkage: str = 'average', distance_threshold: float = None, plot_style: str = 'none',
                 split_plots: bool = False):

        assert linkage.lower() in {'ward', 'complete', 'average', 'single'}
        self.linkage = linkage.lower()
        self.distance_threshold = distance_threshold
        if n_clusters is None:
            if distance_threshold is None:
                raise ValueError('Neither n_clusters or distance_threshold were provided.')
            self.clusterer_kwargs = dict(n_clusters=None, linkage=self.linkage,
                                         distance_threshold=self.distance_threshold)
            super(ClusteringRunnerWithNClusters, self).__init__(data, power_transform, ('affinity', metric), plot_style,
                                                                split_plots)
        else:
            if distance_threshold is not None:
                raise ValueError('Both n_clusters and distance_threshold were provided.')
            self.clusterer_kwargs = dict(linkage=self.linkage)
            super().__init__(data, power_transform, n_clusters, max_n_clusters_estimate, plot_style, split_plots,
                             ('affinity', metric))

    def run(self) -> List[AgglomerativeClustering]:
        self.clusterers = []
        # generate standardized data for plots
        data_for_plot = generic.standardize(self.data)

        if self.distance_threshold is None:
            # iterate over all K values, generating clustering results and plotting them
            for this_n_clusters in self.n_clusters:
                this_n_clusters = int(this_n_clusters)
                clusterer = self.clusterer_class(n_clusters=this_n_clusters, **self.clusterer_kwargs).fit(
                    self.transform(self.data.values))
                # get cluster centers
                centers = np.array(
                    [data_for_plot.values[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(this_n_clusters)])
                # plot results
                title = f"Results of Hierarchical Clustering for n_clusters={this_n_clusters}, "
                f"metric='{self.metric}', \nlinkage='{self.linkage}', power_transform={self.power_transform}"
                self.clusterers.append(clusterer)
                self.plot_clustering(this_n_clusters, data_for_plot, clusterer.labels_, centers, title=title)
        else:
            clusterer = self.clusterer_class(**self.clusterer_kwargs).fit(self.transform(self.data.values))
            # get cluster centers
            centers = np.array(
                [data_for_plot.values[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(clusterer.n_clusters_)])
            self.clusterers.append(clusterer)
            # plot results
            title = f"Results of Hierarchical Clustering for distance_threshold={self.distance_threshold}, "
            f"\nmetric='{self.metric}', linkage='{self.linkage}', " f"power_transform={self.power_transform}"
            self.plot_clustering(clusterer.n_clusters_, data_for_plot, clusterer.labels_, centers, title=title)

        return self.clusterers


class HDBSCANRunner(ClusteringRunner):
    clusterer_class = HDBSCAN

    def __init__(self, data, power_transform: bool, min_cluster_size: int = 5, min_samples: int = 1,
                 metric: str = 'euclidean', cluster_selection_epsilon: float = 0, cluster_selection_method: str = 'eom',
                 return_probabilities: bool = False, plot_style: str = 'none', split_plots: bool = False):
        self.min_cluster_size = min_cluster_size
        self.min_samples = min_samples
        self.cluster_selection_epsilon = cluster_selection_epsilon
        self.cluster_selection_method = cluster_selection_method
        self.return_probabilities = return_probabilities
        self.clusterer_kwargs = dict(min_cluster_size=self.min_cluster_size, min_samples=self.min_samples,
                                     cluster_selection_epsilon=self.cluster_selection_epsilon,
                                     cluster_selection_method=self.cluster_selection_method)
        super().__init__(data, power_transform, ('metric', metric), plot_style, split_plots)

    def run(self):
        self.clusterers = []
        # generate standardized data for plots
        data_for_plot = generic.standardize(self.data)
        # calculate clustering result
        clusterer = self.clusterer_class(**self.clusterer_kwargs)
        clusterer.fit(self.transform(self.data.values))
        self.clusterers.append(clusterer)
        # extract clustering result
        n_clusters = clusterer.labels_.max() + 1
        probabilities = clusterer.probabilities_
        unclustered = np.count_nonzero(clusterer.labels_ == -1)

        if n_clusters != 0:
            # get cluster centers
            centers = np.array([self.data.loc[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(n_clusters)])
            # plot results
            title = f"Results of HDBSCAN Clustering for min_cluster_size={self.min_cluster_size}, "
            f"min_samples = {self.min_samples}, metric='{self.metric}', \nepsilon={self.cluster_selection_epsilon}, "
            f"method='{self.cluster_selection_method}'and power_transform={self.power_transform}"
            self.plot_clustering(n_clusters, data_for_plot, clusterer.labels_, centers, title=title)

        if self.return_probabilities:
            return self.clusterers, probabilities
        return self.clusterers


class CLICOMRunner(ClusteringRunner):
    def __init__(self, data, power_transform: bool, plot_style: str, split_plots: bool,
                 ):
        super(CLICOMRunner, self).__init__(data, power_transform, plot_style=plot_style, split_plots=split_plots)
        # TODO: which arguments does ensemble clustering get???
        self.clusterers: CLICOM

    def run(self):
        pass


class CLICOM:
    def __init__(self, clustering_solutions: List[np.ndarray], threshold: float, cluster_wise_cliques: bool = False):
        self.clustering_solutions = clustering_solutions
        self.threshold = threshold
        self.cluster_wise_cliques = cluster_wise_cliques
        if self.cluster_wise_cliques:
            self.adj_mat = self.get_cluster_similarity_matrix()
        else:
            self.adj_mat = self.get_evidence_accumulation_matrix()
        self.binary_mat = self.adj_mat >= self.threshold
        self.clique_set = set()
        self.clustering_solution = None
        self.labels = None

    def run(self):
        self.find_cliques()
        self.cliques_to_clusters()
        if self.cluster_wise_cliques:
            self.majority_voter()
        else:
            raise NotImplementedError
        self.format_output()

    def find_cliques(self):
        # fast_cliquer algorithm
        pass

    def cliques_to_clusters(self, allowed_overlap: float = 0.1):
        pass

    def majority_voter(self):
        pass

    def format_output(self):
        pass

    def get_evidence_accumulation_matrix(self) -> np.ndarray:
        pass

    def get_cluster_similarity_matrix(self) -> np.ndarray:
        pass

    def inter_cluster_similarity(self, a: Tuple[int, int], b: Tuple[int, int]) -> float:
        # ECS method
        pass

    def cumulative_cluster_similarity(self, clique: set) -> float:
        # CECS method
        pass
