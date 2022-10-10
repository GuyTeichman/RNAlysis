import functools
import itertools
import warnings
from abc import ABC
from typing import List, Set, Tuple, Union, Callable

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pairwisedist as pwdist
import pandas as pd
import seaborn as sns
import sklearn.metrics.pairwise as sklearn_pairwise
from grid_strategy import strategies
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances, silhouette_score
from sklearn_extra.cluster import KMedoids
from tqdm.auto import tqdm

from rnalysis.utils import generic, parsing, validation

try:
    import hdbscan

    HAS_HDBSCAN = True
except ImportError:
    HAS_HDBSCAN = False


class BinaryFormatClusters:
    __slots__ = {'clustering_solutions': 'list of clustering solutions dividing a set of items into clusters',
                 'len_index': 'the number of clusters in each clustering solution',
                 'n_clusters': 'the total number of clusters within all clustering solutions',
                 'n_features': 'the total number of items/genomic features which are split into clusters',
                 'n_solutions': "the number of clustering solutions contained within 'clustering_solutions'",
                 'cluster_sets': "a set representation of all clusters"}

    def __init__(self, clustering_solutions: List[np.ndarray] = None):
        if clustering_solutions is None:
            self.clustering_solutions = None
            self.len_index = None
            self.n_clusters = None
            self.n_features = None
            self.n_solutions = None
            self.cluster_sets = None
        else:
            self._validate_clustering_solutions(clustering_solutions)
            self.clustering_solutions = clustering_solutions
            self.len_index = [sol.shape[0] for sol in clustering_solutions]
            self.n_clusters = sum(self.len_index)
            self.n_features = self.clustering_solutions[0].shape[1]
            self.n_solutions = len(self.clustering_solutions)
            self.cluster_sets: List[Set[int]] = [set(np.where(cluster)[0]) for cluster in self]

    @staticmethod
    def _validate_clustering_solutions(clustering_solutions: List[np.ndarray]):
        assert isinstance(clustering_solutions, list), \
            f"'clustering_solutions' must be a list; instead got {type(clustering_solutions)}"
        assert len(clustering_solutions) > 0, f"'clustering_solutions' must contain at least one clustering solution"
        assert validation.isinstanceiter(clustering_solutions, np.ndarray), \
            f"'clustering_solutions' must exclusively contain numpy arrays"
        for solution in clustering_solutions:
            # in each clustering solution, every feature must be included in one cluster at most
            # clustering algorithms such as HDBSCAN can classify some features as 'noise',
            # and therefore we accept features that are included in 0 clusters.
            assert np.all(solution.sum(axis=0) <= 1), \
                'each feature must be included in one cluster at most per clustering solution'

    def __copy__(self):
        new_obj = type(self)(None)
        new_obj.clustering_solutions = [sol.copy() for sol in self.clustering_solutions] \
            if self.clustering_solutions is not None else None
        new_obj.len_index = self.len_index.copy() if self.len_index is not None else None
        new_obj.n_clusters = self.n_clusters
        new_obj.n_features = self.n_features
        new_obj.n_solutions = self.n_solutions
        new_obj.cluster_sets = [this_set.copy() for this_set in self.cluster_sets] \
            if self.cluster_sets is not None else None

        return new_obj

    def __repr__(self):
        return f'BinaryFormatClusters({repr(self.clustering_solutions)})'

    def __str__(self):
        return f'BinaryFormatClusters({str(self.clustering_solutions)})'

    def __len__(self):
        return len(self.clustering_solutions)

    def __getitem__(self, item: int):
        assert isinstance(item, int)
        assert item >= 0
        assert item < self.n_clusters

        current_ind = 0
        for l_ind, l in enumerate(self.len_index):
            if l + current_ind <= item:
                current_ind += l
            else:
                return self.clustering_solutions[l_ind][item - current_ind, :]

    def __iter__(self):
        for solution in self.clustering_solutions:
            for i in range(solution.shape[0]):
                yield solution[i, :]


class ArbitraryClusterer:
    def __init__(self, labels, n_clusters, **kwargs):
        self.n_clusters_ = n_clusters
        self.labels_ = labels
        for key, val in kwargs.items():
            self.__dict__[key] = val


class CLICOM:
    def __init__(self, clustering_solutions: BinaryFormatClusters, threshold: float, cluster_wise_cliques: bool = True,
                 cluster_unclustered_features: bool = True, min_cluster_size: int = 15):
        self.clustering_solutions: BinaryFormatClusters = clustering_solutions
        self.threshold: float = threshold
        self.cluster_wise_cliques: bool = cluster_wise_cliques
        self.cluster_unclustered_features = cluster_unclustered_features
        self.min_cluster_size = min_cluster_size

        self.n_features: int = self.clustering_solutions.n_features

        if self.cluster_wise_cliques:
            self.adj_mat = self.get_cluster_similarity_matrix()
        else:
            self.adj_mat = self.get_evidence_accumulation_matrix()
        self.binary_mat = self.adj_mat >= self.threshold
        self.clique_set: Set[frozenset] = set()
        self.labels_ = None
        self.n_clusters_ = None

    def run(self):
        self.find_cliques()

        if self.cluster_wise_cliques:
            blocks = self.cliques_to_blocks()
            self.labels_ = self.majority_voter(blocks)
        else:
            clusters = self.cliques_to_clusters()
            self.labels_ = self.clusters_to_labels(clusters)
        self.n_clusters_ = int(self.labels_.max()) + 1

    def clusters_to_labels(self, clusters: List[Set[int]]) -> np.ndarray:
        labels = np.zeros((self.clustering_solutions.n_features,)) - 1
        for i, cluster in enumerate(clusters):
            for feature in cluster:
                labels[feature] = i
        return labels

    def find_cliques(self):
        # fast_cliquer algorithm:
        n_objs = self.binary_mat.shape[0]
        # build K0
        clique_mat = np.zeros_like(self.adj_mat, dtype='object')
        for i in range(n_objs):
            for j in range(n_objs):
                if self.binary_mat[i, j]:
                    clique_mat[i, j] = {i, j}
                else:
                    clique_mat[i, j] = set()
        # create neighbors set
        neighbor_sets = [set() for _ in range(n_objs)]
        for i in range(n_objs):
            for j in range(n_objs):
                if self.binary_mat[i, j] and i != j:
                    neighbor_sets[i].add(j)

        # sequentially construct K1..K|V|
        with tqdm(total=n_objs ** 2, desc='Finding cliques') as pbar:
            for pivot in range(n_objs):
                for i in range(n_objs):
                    for j in range(i, n_objs):
                        # empty-set entries stay the same; if Kij(p-1) is subset of neighbors[p], add p to Kij(p-1)
                        if len(clique_mat[i, j]) > 0 and clique_mat[i, j].issubset(neighbor_sets[pivot]):
                            clique_mat[i, j].add(pivot)
                    pbar.update(1)

        # extract cliques
        for i in range(self.binary_mat.shape[0]):
            for j in range(i + 1, self.binary_mat.shape[0]):
                if len(clique_mat[i, j]) > 1:  # only extract cliques with 2 or more members
                    self.clique_set.add(frozenset(clique_mat[i, j]))

    def cliques_to_clusters(self, allowed_overlap: float = 0.2) -> List[Set[int]]:
        sorted_cliques = [set(clique) for clique in
                          sorted(self.clique_set, reverse=True,
                                 key=lambda clique: (len(clique), sorted(clique)))]
        all_features = set(range(self.binary_mat.shape[0]))
        assigned = set()
        clusters = []
        for m in range(len(sorted_cliques)):
            if len(assigned.intersection(sorted_cliques[m])) / len(sorted_cliques[m]) <= allowed_overlap:
                features_to_assign = sorted_cliques[m] - assigned
                clusters.append(features_to_assign)
                assigned = assigned.union(features_to_assign)

        if self.cluster_unclustered_features:
            for feature in all_features - assigned:
                best_match = None
                best_score = 0
                for i, cluster in enumerate(clusters):
                    this_score = self.feature_cluster_similarity(feature, cluster)
                    if this_score > best_score:
                        best_score = this_score
                        best_match = i
                    if best_score >= self.threshold:
                        clusters[best_match].add(feature)

        filtered_clusters = []
        for cluster in clusters:
            if len(cluster) >= self.min_cluster_size:
                filtered_clusters.append(cluster)

        return filtered_clusters

    def feature_cluster_similarity(self, feature: int, cluster: Set[int]):
        assert not self.cluster_wise_cliques
        # the similarity between a features and a cluster is the mean similarity to every feature in the cluster
        # (equivalent to average linkage criterion)
        return np.mean([self.adj_mat[feature, i] for i in cluster])

    def cliques_to_blocks(self, allowed_overlap: float = 0.2) -> list:
        cluster_blocks = []
        sorted_cliques = sorted(self.clique_set, reverse=True, key=self.cumulative_cluster_similarity)
        all_clusters = set(range(self.binary_mat.shape[0]))
        assigned = set()
        for m in range(len(sorted_cliques)):
            if len(assigned.intersection(sorted_cliques[m])) / len(sorted_cliques[m]) <= allowed_overlap:
                cluster_blocks.append(sorted_cliques[m] - assigned)
                assigned = assigned.union(cluster_blocks[-1])

        if self.cluster_unclustered_features:
            for cluster in all_clusters - assigned:
                nearest_cluster = np.argmax(self.adj_mat[cluster, :])
                if not self.binary_mat[cluster, nearest_cluster]:
                    # only add object to a cluster block if it is sufficiently close
                    continue
                for i, block in enumerate(cluster_blocks):
                    if nearest_cluster in block:
                        cluster_blocks[i] = block.union({cluster})
        return cluster_blocks

    def majority_voter(self, cluster_blocks: list) -> np.ndarray:
        if len(cluster_blocks) == 0:  # if no cluster blocks were found, return an empty clustering result
            labels = np.zeros((self.n_features,)) - 1
            return labels
        voting_mat = np.zeros((len(cluster_blocks), self.n_features))
        for k, block in enumerate(cluster_blocks):
            for cluster in block:
                voting_mat[k, :] += self.clustering_solutions[cluster]
        labels = np.argmax(voting_mat, axis=0)
        # features that got 0 votes should not be clustered
        labels[np.max(voting_mat, axis=0) == 0] = -1
        this_cluster = 0
        for cluster in np.unique(labels):
            if cluster == -1:
                continue
            if sum(labels == cluster) >= self.min_cluster_size:
                labels[labels == cluster] = this_cluster
                this_cluster += 1
            else:
                labels[labels == cluster] = -1
        return labels

    def get_evidence_accumulation_matrix(self) -> np.ndarray:
        mat = np.zeros((self.n_features, self.n_features))
        for cluster in self.clustering_solutions:
            co_clustered = np.where(cluster)[0]
            for pair in itertools.combinations(co_clustered, 2):
                mat[pair] += 1
                mat[pair[::-1]] += 1
        mat /= self.clustering_solutions.n_solutions

        return mat

    def get_cluster_similarity_matrix(self) -> np.ndarray:
        n_clusters = self.clustering_solutions.n_clusters
        mat = np.zeros((n_clusters, n_clusters))
        # only calculate pairwise distance for the upper triangle of the array, to avoid doing double work
        indices = [(i, j) for i, j in zip(*np.triu_indices(n_clusters, 1))]

        n_calculations = (n_clusters ** 2) // 2
        batch_size = 1 + n_calculations // joblib.cpu_count()
        similarities = generic.ProgressParallel(n_jobs=-1, batch_size=batch_size,
                                                desc='Generating cluster similarity matrix')(
            joblib.delayed(CLICOM.inter_cluster_similarity)(*ind, self.clustering_solutions.cluster_sets,
                                                            self.n_features, len(self.clustering_solutions)) for ind in
            indices)
        # populate the upper triangle of the matrix
        np.put(mat, np.ravel_multi_index(np.triu_indices(n_clusters, 1), (n_clusters, n_clusters)), similarities)
        # copy the upper triangle of the matrix to the lower triangle of the matrix
        mat = mat + mat.T
        return mat

    @staticmethod
    def inter_cluster_similarity(a: int, b: int, cluster_sets: List[set], n_features: int, n_solutions: int) -> float:
        # ECS method
        cumsum = 0
        a_set = cluster_sets[a]
        b_set = cluster_sets[b]
        factor = n_solutions * len(a_set) * len(b_set)
        a_b_union = a_set.union(b_set)
        a_b_intersect = a_set.intersection(b_set)
        a_comp_b_intersect = a_set.intersection({i for i in range(n_features) if i not in b_set})
        comp_a_b_intersect = b_set.intersection({i for i in range(n_features) if i not in a_set})
        for cluster in cluster_sets:
            cumsum += generic.combination(len(a_b_union.intersection(cluster)), 2)
            cumsum -= generic.combination(len(a_comp_b_intersect.intersection(cluster)), 2)
            cumsum -= generic.combination(len(comp_a_b_intersect.intersection(cluster)), 2)
            cumsum += generic.combination(len(a_b_intersect.intersection(cluster)), 2)
        cumsum += n_solutions * len(a_b_intersect)
        return cumsum / factor

    def cumulative_cluster_similarity(self, clique: frozenset) -> float:
        # CECS method
        factor = len(clique) / generic.combination(len(clique), 2)
        cumsum = 0
        clique_members = list(clique)
        for i in range(len(clique)):
            for j in range(i + 1, len(clique)):
                cumsum += self.adj_mat[clique_members[i], clique_members[j]]
        return factor * cumsum


class KMedoidsIter:
    """
    A K-Medoids clusterer object with a similar API to other scikit-learn clusterers,
    with the added capability to run the algorithm n_init times and pick the highest-scoring result.
    """

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
        self.titles: list = []
        self.centers: list = []

        self.data_for_plot = None
        self.n_clusters_determine_plot = None

        self.data: pd.DataFrame = data
        self.power_transform: bool = power_transform
        transform = generic.standard_box_cox if power_transform else generic.standardize
        self.transform: Callable = transform

        if metric is not None:
            metric_parameter_name, metric_name = metric
            metric_name = metric_name.lower()
            validation.validate_clustering_parameters(self.legal_metrics.union(self.precomputed_metrics),
                                                      metric_name)
            if metric_name in self.precomputed_metrics:
                def precomputed_transform(x):
                    if isinstance(x, pd.DataFrame):
                        x = x.values
                    return self.precomputed_metrics[metric_name](transform(x))

                self.transform = precomputed_transform
                self.metric = 'precomputed'
            else:
                self.metric = metric_name.lower()
            self.clusterer_kwargs[metric_parameter_name] = self.metric
            self.metric_name = metric_name

        self.plot_style: str = plot_style.lower()
        self.split_plots: bool = split_plots

    def run(self):
        raise NotImplementedError

    @staticmethod
    def sort_clusters(clusterer, n_clusters: int, **arrays_to_sort: np.ndarray):
        updated_labels = clusterer.labels_.copy()
        sorted_arrays = {key: arr.copy() for key, arr in arrays_to_sort.items()}
        cluster_sizes = [sum(clusterer.labels_ == i) for i in np.unique(clusterer.labels_) if i >= 0]
        sorted_label_indices = [ind for _, ind in sorted(zip(cluster_sizes, range(len(cluster_sizes))), reverse=True)]
        for new_ind, orig_ind in enumerate(sorted_label_indices):
            updated_labels[clusterer.labels_ == orig_ind] = new_ind
            for sorted_arr, orig_arr in zip(sorted_arrays.values(), arrays_to_sort.values()):
                sorted_arr[new_ind] = orig_arr[orig_ind]
        return ArbitraryClusterer(updated_labels, n_clusters, **sorted_arrays)

    def _pca_plot(self, n_clusters: int, data: pd.DataFrame, labels: np.ndarray, title: str):
        data_standardized = generic.standardize(data)
        n_components = 2
        pca_obj = PCA(n_components=n_components)
        pcomps = pca_obj.fit_transform(data_standardized)

        columns = [f'Principal component {i + 1}' for i in range(n_components)]
        principal_df = pd.DataFrame(data=pcomps, columns=columns)
        final_df = principal_df
        final_df['labels'] = pd.Series(labels)

        pc_var = pca_obj.explained_variance_ratio_

        pc1_var = pc_var[0]
        pc2_var = pc_var[1]
        final_df = final_df[final_df['labels'] != -1]
        final_df = final_df[['Principal component 1', f'Principal component 2', 'labels']]
        fig = plt.figure(figsize=(9, 9))
        ax = fig.add_subplot(1, 1, 1)
        ax.grid(True)
        ax.set_xlabel(f'{final_df.columns[0]} (explained {pc1_var * 100 :.2f}%)', fontsize=15)
        ax.set_ylabel(f'{final_df.columns[1]} (explained {pc2_var * 100 :.2f}%)', fontsize=15)
        ax.set_title(title, fontsize=16)

        color_generator = generic.color_generator()
        color_opts = [next(color_generator) for _ in range(n_clusters)]
        for cluster in range(n_clusters):
            ax.scatter(final_df[final_df['labels'] == cluster].iloc[:, 0],
                       final_df[final_df['labels'] == cluster].iloc[:, 1],
                       label=f'Cluster {cluster + 1}', c=color_opts[cluster], s=20, alpha=0.4)
        ax.legend(title="Clusters")
        ax.grid(True)
        plt.tight_layout()
        plt.show()

    def _generate_plots(self, n_clusters: int, data: pd.DataFrame, labels: np.ndarray, centers: np.ndarray,
                        title: str) -> Union[Tuple[plt.Figure, List[plt.Axes]], None]:
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
            ylabel_fontsize = 10 if n_clusters > 15 else 8
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
                    vals = data.loc[labels == i, :].T.values
                    axes[-1].plot(x, vals, color=color, alpha=0.05, linewidth=0.35)

            axes[-1].set_title(f"Cluster number {i + 1} ({np.count_nonzero(labels == i)} genes)")

            ylabel = 'Standardized\nexpression' if n_clusters > 15 else 'Standardized expression'
            axes[-1].set_ylabel(ylabel, fontsize=ylabel_fontsize)
            axes[-1].set_xticks(x)
            axes[-1].set_xticklabels(data.columns, fontsize=ticklabels_fontsize)
        for ax in axes:
            ax.set_ylim((min_y, max_y))
        if not self.split_plots:
            fig.tight_layout(rect=[0, 0.03, 1, 0.92])
        plt.show()
        self._pca_plot(n_clusters, data, labels, title)
        return fig, axes

    def plot_clustering(self):
        data = self.data_for_plot
        for clusterer, centers, title in zip(self.clusterers, self.centers, self.titles):
            self._generate_plots(clusterer.n_clusters_, data, clusterer.labels_, centers, title)
        if self.n_clusters_determine_plot is not None:
            self.n_clusters_determine_plot()

    def clustering_labels_to_binary_format(self) -> np.ndarray:
        pass


class ClusteringRunnerWithNClusters(ClusteringRunner, ABC):
    def __init__(self, data: pd.DataFrame, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'auto', plot_style: str = 'none',
                 split_plots: bool = False, metric: Tuple[str, str] = None):
        super().__init__(data, power_transform, metric, plot_style, split_plots)
        self.max_n_clusters_estimate = min(20, data.shape[
            0] // 4) if max_n_clusters_estimate == 'auto' else max_n_clusters_estimate
        self.n_clusters: list = self.parse_n_clusters(n_clusters)

    def parse_n_clusters(self, n_clusters) -> List[int]:
        # get best n_clusters value using gap statistic/silhouette method, if requested by the user
        if isinstance(n_clusters, str) and n_clusters.lower() == 'silhouette':
            n_clusters = self.silhouette()
        elif isinstance(n_clusters, str) and n_clusters.lower() in {'gap', 'gap statistic', 'gap_statistic'}:
            n_clusters = self.gap_statistic()

        n_clusters = parsing.data_to_list(n_clusters)
        # make sure all values of n_clusters are in valid range
        n_clusters, n_clusters_copy = itertools.tee(n_clusters)
        n_clusters = parsing.data_to_list(n_clusters)
        assert np.all([isinstance(item, int) and 2 <= item <= len(self.data.index) for item in n_clusters_copy]), \
            f"Invalid value for n_clusters: '{n_clusters}'. n_clusters must be 'gap', 'silhouette', " \
            f"or an integer/Iterable of integers in range 2 <= n_clusters <= n_features."
        return n_clusters

    def gap_statistic(self, random_seed: int = None, n_refs: int = 10) -> int:
        raw_data = self.data.values
        # determine the range of k values to be tested
        n_clusters_range = np.arange(1, self.max_n_clusters_estimate + 1)
        print(f"Estimating the optimal number of clusters using the Gap Statistic method in range "
              f"{2}:{self.max_n_clusters_estimate}...")
        # set the random state, if one was supplied
        if random_seed is not None:
            np.random.seed(random_seed)
        # calculate the SVD of the observed data (X), and transform via X' = x_tag = dot(X, V)
        # note: the data is centered by sklearn.decomposition.PCA, no need to pre-center it.
        pca = PCA(random_state=random_seed).fit(raw_data)
        x_tag = pca.transform(raw_data)
        # determine the ranges of the columns of X' (x_tag)
        a, b = x_tag.min(axis=0, keepdims=True), x_tag.max(axis=0, keepdims=True)
        # transform the observed data using Box-Cox, and then standardize it
        data = self.transform(raw_data)
        # generate 'n_refs' random reference arrays:
        # draw uniform features Z' over the ranges of the columns of X', and back-transform via Z = dot(Z', V.T), then
        # transform the random reference data using Box-Cox, and then standardize it
        refs = [self.transform(generic.shift_to_baseline(
            pca.inverse_transform(np.random.random_sample(size=raw_data.shape) * (b - a) + a))) for _
            in range(n_refs)]
        # allocate empty arrays for observed/expected log(inertia), gap scores Gap(K) and gap error S(K)
        log_disp_obs = np.zeros((len(n_clusters_range)))
        log_disp_exp = np.zeros((len(n_clusters_range)))
        gap_scores = np.zeros((len(n_clusters_range)))
        gap_error = np.zeros((len(n_clusters_range)))

        # iterate over K values in the range
        for ind, this_n_clusters in enumerate(
            tqdm(n_clusters_range, desc="Testing 'n_clusters' values", unit='values')):
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
        # 1 is not a viable candidate
        if k_candidates[0] == 1:
            k_candidates = k_candidates[1:]
        # if no K fulfilled the condition (the gap function is monotonically increasing),
        # return the last value of K as the best K
        best_n_clusters = int(k_candidates[0]) if len(k_candidates) > 0 else int(n_clusters_range[-1])
        n_clusters_ind = np.argmax(n_clusters_range == best_n_clusters)
        print(f"Using the Gap Statistic method, {best_n_clusters} was chosen as the best number of clusters (K).")
        if len(k_candidates) >= 2:
            print(f"Other potentially good values of K that were found: {list(k_candidates[1:])}. ")

        # plot results
        self.n_clusters_determine_plot = functools.partial(self._plot_gap_statistic, n_clusters_range, log_disp_obs,
                                                           log_disp_exp, gap_scores, gap_error, best_n_clusters,
                                                           n_clusters_ind)

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
        plt.tight_layout()
        plt.show()
        return fig

    def silhouette(self) -> int:
        print(f"Estimating the optimal number of clusters using the Silhouette method in range "
              f"{2}:{self.max_n_clusters_estimate}...")
        data = self.transform(self.data.values)
        silhouette_scores = []
        n_clusters_range = np.arange(2, self.max_n_clusters_estimate + 1)
        for n_clusters in n_clusters_range:
            clusterer = self.clusterer_class(n_clusters=n_clusters, **self.clusterer_kwargs)
            silhouette_scores.append(silhouette_score(data, clusterer.fit_predict(data)))
        best_n_clusters = int(n_clusters_range[int(np.argmax(silhouette_scores))])

        self.n_clusters_determine_plot = functools.partial(self._plot_silhouette, best_n_clusters, n_clusters_range,
                                                           silhouette_scores)
        print(f"Using the Silhouette method, {best_n_clusters} was chosen as the best number of clusters (k).")
        return best_n_clusters

    @staticmethod
    def _plot_silhouette(best_n_clusters, n_clusters_range, silhouette_scores):
        fig = plt.figure(figsize=(7, 9))
        ax = fig.add_subplot(111)
        ax.plot(n_clusters_range, silhouette_scores, '--o', color='r')
        ax.set_ylim(-1.0, 1.0)
        ax.set_title("Silhouette method for estimating optimal number of clusters")
        ax.set_ylabel('Average silhouette score')
        ax.set_xlabel("Number of clusters (k)")
        ax.axhline(0, color='grey', linewidth=2, linestyle='--')
        ax.annotate(f'Best n_clusters={best_n_clusters}', xy=(best_n_clusters, max(silhouette_scores)),
                    xytext=(best_n_clusters + 1, 0.75),
                    arrowprops=dict(facecolor='black', shrink=0.15))
        ax.set_xticks(n_clusters_range)
        sns.despine()
        plt.tight_layout()
        return fig


class KMeansRunner(ClusteringRunnerWithNClusters):
    clusterer_class = KMeans

    def __init__(self, data, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'auto', random_seed: int = None, n_init: int = 3,
                 max_iter: int = 300, plot_style: str = 'none', split_plots: bool = False):
        assert isinstance(random_seed, int) or random_seed is None
        assert isinstance(n_init, int) and n_init >= 1
        assert isinstance(max_iter, int) and max_iter >= 1
        self.random_seed = random_seed
        self.n_init = n_init
        self.max_iter = max_iter
        self.clusterer_kwargs = dict(n_init=self.n_init, max_iter=self.max_iter, random_state=self.random_seed)

        super(KMeansRunner, self).__init__(data, power_transform, n_clusters, max_n_clusters_estimate, plot_style,
                                           split_plots)

    def run(self, plot: bool = True) -> List[ArbitraryClusterer]:
        self.clusterers = []
        # generate standardized data for plots
        data_for_plot = generic.standard_box_cox(self.data) if self.power_transform else generic.standardize(self.data)
        self.data_for_plot = data_for_plot
        # iterate over all K values, generating clustering results and plotting them
        data_for_algo = self.transform(self.data)
        for this_n_clusters in self.n_clusters:
            this_n_clusters = int(this_n_clusters)
            clusterer = self.clusterer_class(n_clusters=this_n_clusters, **self.clusterer_kwargs).fit(data_for_algo)
            clusterer = self.sort_clusters(clusterer, this_n_clusters, cluster_centers_=clusterer.cluster_centers_)
            # plot results
            centers = clusterer.cluster_centers_
            title = f"Results of K-Means Clustering for n_clusters={this_n_clusters} and " \
                    f"power_transform={self.power_transform}"
            self.centers.append(centers)
            self.titles.append(title)
            self.clusterers.append(clusterer)
        if plot:
            self.plot_clustering()
        return self.clusterers


class KMedoidsRunner(ClusteringRunnerWithNClusters):
    clusterer_class = KMedoidsIter
    legal_metrics = set(sklearn_pairwise._VALID_METRICS)

    def __init__(self, data, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'auto', metric: str = 'euclidean',
                 random_seed: int = None, n_init: int = 3, max_iter: int = 300, plot_style: str = 'none',
                 split_plots: bool = False):
        assert isinstance(random_seed, int) or random_seed is None
        assert isinstance(n_init, int) and n_init >= 1
        assert isinstance(max_iter, int) and max_iter >= 1
        assert isinstance(metric, str)

        self.random_seed = random_seed
        self.n_init = n_init
        self.max_iter = max_iter
        self.clusterer_kwargs = dict(n_init=self.n_init, max_iter=self.max_iter, random_state=self.random_seed)
        super(KMedoidsRunner, self).__init__(data, power_transform, n_clusters, max_n_clusters_estimate, plot_style,
                                             split_plots, ('metric', metric))

    def run(self, plot: bool = True) -> List[ArbitraryClusterer]:
        self.clusterers = []
        # generate standardized data for plots
        data_for_plot = generic.standard_box_cox(self.data) if self.power_transform else generic.standardize(self.data)
        self.data_for_plot = data_for_plot
        # iterate over all K values, generating clustering results and plotting them
        data_for_algo = self.transform(self.data)
        for this_n_clusters in self.n_clusters:
            this_n_clusters = int(this_n_clusters)
            clusterer = self.clusterer_class(n_clusters=this_n_clusters, **self.clusterer_kwargs).fit(data_for_algo)
            clusterer = self.sort_clusters(clusterer, this_n_clusters, medoid_indices_=clusterer.medoid_indices_)
            # get cluster centers
            centers = data_for_plot.values[clusterer.medoid_indices_, :]
            # plot results
            title = f"Results of K-Medoids Clustering for n_clusters={this_n_clusters}, " \
                    f"metric='{self.metric_name}', power_transform={self.power_transform}"
            self.centers.append(centers)
            self.titles.append(title)
            self.clusterers.append(clusterer)
        if plot:
            self.plot_clustering()
        return self.clusterers


class HierarchicalRunner(ClusteringRunnerWithNClusters):
    clusterer_class = AgglomerativeClustering
    legal_metrics = set(sklearn_pairwise.PAIRED_DISTANCES.keys())

    def __init__(self, data, power_transform: bool, n_clusters: Union[int, List[int], str],
                 max_n_clusters_estimate: Union[int, str] = 'auto', metric: str = 'euclidean',
                 linkage: str = 'average', distance_threshold: float = None, plot_style: str = 'none',
                 split_plots: bool = False):
        assert isinstance(linkage, str), f"'linkage' must be of type str, instead got {type(linkage)}."
        assert isinstance(metric, str), f"'metric' must be of type str, instead got {type(metric)}."
        assert linkage.lower() in {'ward', 'complete', 'average',
                                   'single'}, f"Invalid linkage method '{linkage.lower()}'."
        assert isinstance(distance_threshold, (int, float)) or distance_threshold is None
        if linkage.lower() == 'ward':
            assert metric.lower() == 'euclidean', "metric must be 'euclidean' if linkage is 'ward'."

        self.linkage = linkage.lower()
        self.distance_threshold = distance_threshold
        if n_clusters is None or n_clusters == 'distance':
            if distance_threshold is None:
                raise ValueError("Neither 'n_clusters' or 'distance_threshold' were provided. ")
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

    def run(self, plot: bool = True) -> List[ArbitraryClusterer]:
        self.clusterers = []
        # generate standardized data for plots
        data_for_plot = generic.standard_box_cox(self.data) if self.power_transform else generic.standardize(self.data)
        self.data_for_plot = data_for_plot
        data_for_algo = self.transform(self.data)
        if self.distance_threshold is None:
            # iterate over all K values, generating clustering results and plotting them
            for this_n_clusters in self.n_clusters:
                this_n_clusters = int(this_n_clusters)
                clusterer = self.clusterer_class(n_clusters=this_n_clusters, **self.clusterer_kwargs).fit(data_for_algo)
                clusterer = self.sort_clusters(clusterer, this_n_clusters)
                # get cluster centers
                centers = np.array(
                    [data_for_plot.values[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(this_n_clusters)])
                # plot results
                title = f"Results of Hierarchical Clustering for n_clusters={this_n_clusters}, " \
                        f"metric='{self.metric_name}', \nlinkage='{self.linkage}', power_transform={self.power_transform}"
                self.centers.append(centers)
                self.titles.append(title)
                self.clusterers.append(clusterer)
            if plot:
                self.plot_clustering()
        else:
            clusterer = self.clusterer_class(**self.clusterer_kwargs).fit(self.transform(self.data.values))
            clusterer = self.sort_clusters(clusterer, len(np.unique(clusterer.labels_)))
            # get cluster centers
            centers = np.array(
                [data_for_plot.values[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(clusterer.n_clusters_)])
            self.clusterers.append(clusterer)
            # plot results
            title = f"Results of Hierarchical Clustering for distance_threshold={self.distance_threshold}, \n" \
                    f"metric='{self.metric_name}', linkage='{self.linkage}', " f"power_transform={self.power_transform}"
            self.centers.append(centers)
            self.titles.append(title)
            if plot:
                self.plot_clustering()

        return self.clusterers


class HDBSCANRunner(ClusteringRunner):
    clusterer_class = hdbscan.HDBSCAN if HAS_HDBSCAN else None
    legal_metrics = set(hdbscan.dist_metrics.METRIC_MAPPING.keys()) if HAS_HDBSCAN else set()

    def __init__(self, data, power_transform: bool, min_cluster_size: int = 5, min_samples: int = 1,
                 metric: str = 'euclidean', cluster_selection_epsilon: float = 0, cluster_selection_method: str = 'eom',
                 return_probabilities: bool = False, plot_style: str = 'none', split_plots: bool = False):
        self.return_probabilities = return_probabilities
        if not HAS_HDBSCAN:
            return
        assert isinstance(metric, str)
        assert isinstance(min_cluster_size, int) and min_cluster_size > 1
        assert isinstance(min_samples, int) and min_samples >= 1
        assert isinstance(cluster_selection_epsilon, (int, float))
        assert isinstance(cluster_selection_method, str) and cluster_selection_method.lower() in {'eom', 'leaf'}

        self.min_cluster_size = min_cluster_size
        self.min_samples = min_samples
        self.cluster_selection_epsilon = cluster_selection_epsilon
        self.cluster_selection_method = cluster_selection_method
        self.clusterer_kwargs = dict(min_cluster_size=self.min_cluster_size, min_samples=self.min_samples,
                                     cluster_selection_epsilon=self.cluster_selection_epsilon,
                                     cluster_selection_method=self.cluster_selection_method)
        super().__init__(data, power_transform, ('metric', metric), plot_style, split_plots)

    @staticmethod
    def _missing_dependency_warning():
        warnings.warn(f"Package 'hdbscan' is not installed. \n"
                      f"If you want to use HDBSCAN clustering, please install package 'hdbscan' and try again. ")

    def run(self, plot: bool = True) -> Union[
        List[ArbitraryClusterer], Tuple[List[ArbitraryClusterer], np.ndarray], List[None]]:
        if not HAS_HDBSCAN:
            self._missing_dependency_warning()
            if self.return_probabilities:
                return [None], np.array([])
            return [None]
        self.clusterers = []

        # calculate clustering result
        clusterer = self.clusterer_class(**self.clusterer_kwargs)
        clusterer.fit(self.transform(self.data.values))
        clusterer = self.sort_clusters(clusterer, clusterer.labels_.max() + 1, probabilities_=clusterer.probabilities_)
        self.clusterers.append(clusterer)
        # extract clustering result
        n_clusters = clusterer.n_clusters_
        probabilities = clusterer.probabilities_
        unclustered = np.count_nonzero(clusterer.labels_ == -1)

        if n_clusters != 0:
            # generate standardized data for plots
            data_for_plot = generic.standard_box_cox(self.data) if self.power_transform else generic.standardize(
                self.data)
            self.data_for_plot = data_for_plot
            # get cluster centers
            centers = np.array([data_for_plot.loc[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(n_clusters)])
            # plot results
            title = f"Results of HDBSCAN Clustering for min_cluster_size={self.min_cluster_size}, " \
                    f"min_samples = {self.min_samples}, metric='{self.metric_name}', " \
                    f"\nepsilon={self.cluster_selection_epsilon}, " \
                    f"method='{self.cluster_selection_method}' and power_transform={self.power_transform}"
            self.centers.append(centers)
            self.titles.append(title)
            if plot:
                self.plot_clustering()

        if self.return_probabilities:
            return self.clusterers, probabilities
        return self.clusterers


class CLICOMRunner(ClusteringRunner):
    clusterer_class = CLICOM
    algorithm_mapper = {'kmeans': KMeansRunner, 'k means': KMeansRunner, 'k-means': KMeansRunner,
                        'kmedoids': KMedoidsRunner, 'k medoids': KMedoidsRunner, 'k-medoids': KMedoidsRunner,
                        'hierarchical': HierarchicalRunner, 'agglomerative': HierarchicalRunner,
                        'hdbscan': HDBSCANRunner}

    def __init__(self, data, power_transform: Union[bool, List[bool]], evidence_threshold: float,
                 cluster_unclustered_features: bool, min_cluster_size: int, *parameter_dicts: dict,
                 plot_style: str = 'none',
                 split_plots: bool = False, ):
        self.clustering_solutions: list = []
        self.evidence_threshold = evidence_threshold
        self.cluster_unclustered_features = cluster_unclustered_features
        self.min_cluster_size = min_cluster_size
        self.parameter_dicts = parameter_dicts
        self.clusterers = []
        super(CLICOMRunner, self).__init__(data, power_transform, plot_style=plot_style, split_plots=split_plots)

    def run(self, plot: bool = True) -> List[ArbitraryClusterer]:
        valid_setups = self.find_valid_clustering_setups()
        n_setups = len(valid_setups)
        print(f"Found {n_setups} legal clustering setups.")
        for setup in tqdm(valid_setups, desc='Running clustering setups', unit=' setup'):
            self.clustering_solutions.extend(self.run_clustering_setup(setup))
        solutions_binary_format = self.clusterers_to_binary_format()
        if solutions_binary_format.n_clusters >= 4 * np.sqrt(solutions_binary_format.n_features):
            print("Number of clusters found in all clustering solutions is large relatively to the number of features. "
                  "\nCLICOM clustering will be calculated based on a feature graph instead of a cluster graph "
                  "to reduce computation time.")
            cluster_wise_cliques = False
        else:
            cluster_wise_cliques = True

        clusterer = self.clusterer_class(solutions_binary_format, self.evidence_threshold,
                                         cluster_unclustered_features=self.cluster_unclustered_features,
                                         min_cluster_size=self.min_cluster_size,
                                         cluster_wise_cliques=cluster_wise_cliques)

        print("Calculating CLICOM clustering...", end='\t')
        clusterer.run()
        clusterer = self.sort_clusters(clusterer, clusterer.n_clusters_)
        self.clusterers.append(clusterer)
        print("Done")

        n_clusters = clusterer.n_clusters_
        if n_clusters != 0:
            # get cluster centers
            # generate standardized data for plots
            data_for_plot = generic.standard_box_cox(self.data)
            self.data_for_plot = data_for_plot
            centers = np.array([data_for_plot.loc[clusterer.labels_ == i, :].T.mean(axis=1) for i in range(n_clusters)])
            # plot results
            title = f"Results of Emsemble Clustering for {n_setups} clustering solutions, " \
                    f"\nevidence_threshold={self.evidence_threshold}, and min_cluster_size={self.min_cluster_size}"
            self.centers.append(centers)
            self.titles.append(title)
            if plot:
                self.plot_clustering()

        return self.clusterers

    def find_valid_clustering_setups(self) -> list:
        valid_setups = []
        for d in self.parameter_dicts:
            assert isinstance(d, dict)
            assert 'method' in d, \
                f"Each parameter_dict must contain a 'method' field. Failed to find 'method' field in: {d}"
            assert d['method'].lower() in self.algorithm_mapper, f"Unknown clustering method '{d['method']}'."
            runner_class = self.algorithm_mapper[d['method'].lower()]
            params = d.copy()
            params.pop('method')
            setups_to_test = itertools.product(parsing.data_to_list(self.power_transform),
                                               *[parsing.data_to_list(val) for val in params.values()])
            for setup in setups_to_test:
                kwargs = {key: val for key, val in zip(itertools.chain(['power_transform'], params.keys()), setup)}
                try:
                    runner_class(self.data, **kwargs)
                    valid_setups.append((runner_class, kwargs))
                except AssertionError:
                    continue

        return valid_setups

    def run_clustering_setup(self, setup) -> list:
        runner_class, kwargs = setup
        runner = runner_class(self.data, **kwargs)
        clusterers = runner.run(plot=False)
        return [clusterer.labels_ for clusterer in clusterers]

    def clusterers_to_binary_format(self) -> BinaryFormatClusters:
        lst = []
        for labels in self.clustering_solutions:
            mat = np.zeros((np.max(labels) + 1, len(labels)), dtype=bool)
            for i in np.unique(labels):
                if i == -1:
                    continue
                for j in np.where(labels == i):
                    mat[i, j] = True
            lst.append(mat)

        return BinaryFormatClusters(lst)
