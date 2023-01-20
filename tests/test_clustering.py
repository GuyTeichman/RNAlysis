from collections import namedtuple

import pytest

from rnalysis.utils.clustering import *
from rnalysis.utils.io import load_csv


@pytest.fixture
def valid_clustering_solutions():
    return [np.array([[1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 1]]),
            np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]), np.array(
            [[1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1]])]


@pytest.fixture
def valid_clustering_solutions_with_noise():
    return [np.array([[1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 1]]),
            np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]), np.array(
            [[1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1]])]


@pytest.fixture
def invalid_clustering_solutions():
    return [np.array([[1, 1, 1, 0, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 1]]),
            np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]), np.array(
            [[1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1]])]


def test_kmedoidsiter_api():
    truth = KMedoids(3, max_iter=300, init='k-medoids++', random_state=42)
    kmeds = KMedoidsIter(3, init='k-medoids++', max_iter=300, n_init=1, random_state=42)
    df = load_csv('tests/test_files/counted.csv', 0)
    truth.fit(df)
    kmeds.fit(df)
    assert np.all(truth.cluster_centers_ == kmeds.cluster_centers_)
    assert np.all(truth.inertia_ == kmeds.inertia_)

    assert np.all(truth.predict(df) == kmeds.predict(df))
    assert np.all(truth.fit_predict(df) == kmeds.fit_predict(df))

    kmeds_rand = KMedoidsIter(3, init='k-medoids++', max_iter=300, n_init=3)
    kmeds_rand.fit(df)
    kmeds_rand.predict(df)
    kmeds_rand.fit_predict(df)

    assert repr(kmeds) == repr(kmeds.clusterer)
    assert str(kmeds) == str(kmeds.clusterer)


def test_kmedoidsiter_iter():
    kmeds = KMedoidsIter(3, init='k-medoids++', max_iter=300, n_init=5, random_state=0)
    df = load_csv('tests/test_files/counted.csv', 0)
    kmeds.fit(df)

    inertias = []
    clusterers = []
    for i in range(5):
        clusterers.append(KMedoids(3, max_iter=300, init='k-medoids++', random_state=0).fit(df))
        inertias.append(clusterers[i].inertia_)
    truth_inertia = max(inertias)
    truth_kmeds = clusterers[np.argmax(inertias)]
    assert kmeds.inertia_ == truth_inertia
    assert np.all(kmeds.clusterer.predict(df) == truth_kmeds.predict(df))


@pytest.mark.parametrize("args,expected", [
    ((False, 'silhouette', 20), [6]),
    ((False, 'gap', 20), [7]),
    ((False, [7, 2, 5], 20), [7, 2, 5]),
    ((False, [10], 20), [10]),
    ((False, 10, 20), [10]),
    ((False, range(2, 9), 20), [2, 3, 4, 5, 6, 7, 8])
])
def test_parse_n_clusters(monkeypatch, args, expected):
    data = pd.read_csv('tests/test_files/counted.csv', index_col=0)
    monkeypatch.setattr(ClusteringRunnerWithNClusters, "silhouette", lambda self: 6)
    monkeypatch.setattr(ClusteringRunnerWithNClusters, "gap_statistic", lambda self: 7)

    runner = ClusteringRunnerWithNClusters(data, *args)
    assert runner.n_clusters == expected


@pytest.mark.parametrize("args", [
    (False, [5, 2, '3'], 20),
    (False, [3, 23, 2], 20),
    (False, '17', 20),
    (False, 1, 20),
    (False, [3, 5, 1], 20)
])
def test_parse_n_clusters_invalid_values(args: tuple):
    data = pd.read_csv('tests/test_files/counted.csv', index_col=0)
    with pytest.raises(AssertionError):
        _ = ClusteringRunnerWithNClusters(data, *args)


def test_compute_dispersion():
    clust_with_inertia = namedtuple('Clusterer', ['inertia_'])
    clust_without_inertia = namedtuple('Clusterer', ['labels_'])
    data = pd.read_csv('tests/test_files/counted.csv', index_col=0).values
    for k in [1, 3, data.shape[0]]:
        kmeans = KMeans(k, random_state=42).fit(data)
        assert ClusteringRunnerWithNClusters._compute_dispersion(clust_with_inertia(kmeans.inertia_), data,
                                                                 k) == kmeans.inertia_
        print(k, kmeans.inertia_, sorted(kmeans.labels_))
        assert np.isclose(
            ClusteringRunnerWithNClusters._compute_dispersion(clust_without_inertia(kmeans.labels_), data, k),
            kmeans.inertia_)


def _compare_setups(setups, truth):
    assert len(setups) == len(truth)
    for i in range(len(setups)):
        this_setup = setups[i]
        this_truth = truth[i]
        assert this_setup[0] == this_truth[0]
        assert this_setup[2] == this_truth[2]
        assert this_setup[1][0].equals(this_truth[1][0])


def test_clicomrunner_find_valid_clustering_setups():
    df = pd.DataFrame(np.zeros((10, 10)))
    truth = [(KMeansRunner, (df,), dict(power_transform=True, n_clusters=2, n_init=5)),
             (KMeansRunner, (df,), dict(power_transform=True, n_clusters=3, n_init=5)),
             (KMeansRunner, (df,), dict(power_transform=True, n_clusters=5, n_init=5)),
             (KMeansRunner, (df,), dict(power_transform=False, n_clusters=2, n_init=5)),
             (KMeansRunner, (df,), dict(power_transform=False, n_clusters=3, n_init=5)),
             (KMeansRunner, (df,), dict(power_transform=False, n_clusters=5, n_init=5)),
             (HierarchicalRunner, (df,), dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=5, metric='jackknife', linkage='average')),
             (HierarchicalRunner, (df,), dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=7, metric='jackknife', linkage='average')),
             (HierarchicalRunner, (df,), dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=10, metric='jackknife', linkage='average')),
             (HierarchicalRunner, (df,), dict(power_transform=False, n_clusters=5, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=5, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=5, metric='jackknife', linkage='average')),
             (HierarchicalRunner, (df,), dict(power_transform=False, n_clusters=7, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=7, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=7, metric='jackknife', linkage='average')),
             (
                 HierarchicalRunner, (df,),
                 dict(power_transform=False, n_clusters=10, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=10, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=10, metric='jackknife', linkage='average')),

             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=0.5)),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=1)),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=None, metric='euclidean', distance_threshold=0.5)),
             (HierarchicalRunner, (df,),
              dict(power_transform=False, n_clusters=None, metric='euclidean', distance_threshold=1))]

    runner = CLICOMRunner(pd.DataFrame(np.zeros((10, 10))), [list(range(10))], [True, False], 0.5,
                          False, 1,
                          dict(method='kmeans', n_clusters=[2, 3, 5], n_init=5),
                          dict(method='hierarchical', n_clusters=[5, 7, 10], metric=['euclidean', 'jackknife'],
                               linkage=['ward', 'average', 'doesnt exist']),
                          dict(method='hierarchical', n_clusters=None, distance_threshold=[0.5, 1], metric='euclidean'))
    setups = runner.find_valid_clustering_setups()
    _compare_setups(setups, truth)

    truth = [(KMeansRunner, (df,), dict(power_transform=True, n_clusters=2, n_init=5)),
             (KMeansRunner, (df,), dict(power_transform=True, n_clusters=3, n_init=5)),
             (KMeansRunner, (df,), dict(power_transform=True, n_clusters=5, n_init=5)),
             (HierarchicalRunner, (df,), dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=5, metric='jackknife', linkage='average')),
             (HierarchicalRunner, (df,), dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=7, metric='jackknife', linkage='average')),
             (HierarchicalRunner, (df,), dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=10, metric='jackknife', linkage='average')),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=0.5)),
             (HierarchicalRunner, (df,),
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=1))]

    runner = CLICOMRunner(pd.DataFrame(np.zeros((10, 10))), [list(range(10))], True, 0.5, False, 15,
                          dict(method='kmeans', n_clusters=[2, 3, 5], n_init=5),
                          dict(method='hierarchical', n_clusters=[5, 7, 10], metric=['euclidean', 'jackknife'],
                               linkage=['ward', 'average', 'doesnt exist']),
                          dict(method='hierarchical', n_clusters=None, distance_threshold=[0.5, 1], metric='euclidean'))
    setups = runner.find_valid_clustering_setups()
    _compare_setups(setups, truth)

    truth = [(KMeansRunner, (df[[0, 1, 2, 3, 4]],), dict(power_transform=True, n_clusters=6, n_init=5)),
             (KMeansRunner, (df[[5, 6, 7, 8, 9]],), dict(power_transform=True, n_clusters=6, n_init=5)),
             (KMedoidsRunner, (df[[0, 1, 2, 3, 4]],), dict(power_transform=True, n_clusters=7, n_init=5)),
             (KMedoidsRunner, (df[[5, 6, 7, 8, 9]],), dict(power_transform=True, n_clusters=7, n_init=5))]
    runner = CLICOMRunner(pd.DataFrame(np.zeros((10, 10))), [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]], True, 0.5, False, 15,
                          dict(method='kmeans', n_clusters=6, n_init=5),
                          dict(method='kmedoids', n_clusters=7, n_init=5))
    setups = runner.find_valid_clustering_setups()
    _compare_setups(setups, truth)


def test_find_cliques(monkeypatch):
    truth = {frozenset({2, 4, 7, 8}), frozenset({1, 7}), frozenset({1, 3, 6}), frozenset({0, 3, 5, 6})}
    binary_adj_mat = pd.read_csv('tests/test_files/clicom_binary_adj_matrix.csv', index_col=0).values
    monkeypatch.setattr(CLICOM, "get_cluster_similarity_matrix", lambda self: binary_adj_mat)
    clicom = CLICOM(BinaryFormatClusters([np.array([[0, 1]])]), 1)
    clicom.find_cliques()
    assert clicom.clique_set == truth


def test_clicom_cluster_wise_result(valid_clustering_solutions):
    bin_format = BinaryFormatClusters(valid_clustering_solutions)

    truth = {frozenset({0, 1, 2, 3}), frozenset({4, 5, 6, 7})}
    clusterer = CLICOM(bin_format, 1 / 3, min_cluster_size=1)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth == results

    truth = set()
    clusterer = CLICOM(bin_format, 1 / 3, min_cluster_size=5)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth == results

    truth_threshold_0 = {frozenset({0, 1, 2, 3, 4, 5, 6, 7})}
    clusterer = CLICOM(bin_format, 0, min_cluster_size=1)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth_threshold_0 == results

    truth_threshold_max = set()
    clusterer = CLICOM(bin_format, 1.01, min_cluster_size=1)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth_threshold_max == results


def test_clicom_feature_wise_result(valid_clustering_solutions):
    bin_format = BinaryFormatClusters(valid_clustering_solutions)

    truth = {frozenset({0, 1, 2, 3}), frozenset({4, 5, 6, 7})}
    clusterer = CLICOM(bin_format, 0.33, cluster_wise_cliques=False, min_cluster_size=1)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth == results

    truth = set()
    clusterer = CLICOM(bin_format, 0.33, cluster_wise_cliques=False, min_cluster_size=5)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth == results

    truth_threshold_0 = {frozenset({0, 1, 2, 3, 4, 5, 6, 7})}
    clusterer = CLICOM(bin_format, 0, cluster_wise_cliques=False, min_cluster_size=1)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth_threshold_0 == results

    truth_threshold_max = set()
    clusterer = CLICOM(bin_format, 1.01, cluster_wise_cliques=False, min_cluster_size=1)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth_threshold_max == results


def test_binary_format_clusters_init(valid_clustering_solutions):
    cl = BinaryFormatClusters()
    for attr in cl.__slots__:
        assert getattr(cl, attr) is None

    cl2 = BinaryFormatClusters(valid_clustering_solutions)
    assert cl2.n_solutions == 3
    assert cl2.n_features == 8
    assert cl2.n_clusters == 9
    assert cl2.len_index == [3, 2, 4]
    assert cl2.cluster_sets == [{0, 1}, {2, 3, 4}, {5, 6, 7}, {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1}, {2, 3}, {4, 5},
                                {6, 7}]
    assert cl2.clustering_solutions == valid_clustering_solutions


def test_binary_format_clusters_copy(valid_clustering_solutions):
    cl = BinaryFormatClusters(valid_clustering_solutions)
    cl2 = cl.__copy__()
    for attr in cl.__slots__:
        if attr == 'clustering_solutions':
            for sol1, sol2 in zip(getattr(cl, attr), getattr(cl2, attr)):
                assert np.all(sol1 == sol2)
        else:
            assert np.all(getattr(cl, attr) == getattr(cl2, attr))
        if not isinstance(getattr(cl, attr), (int, float)):
            assert getattr(cl, attr) is not getattr(cl2, attr)

    _ = BinaryFormatClusters(None).__copy__()


def test_binary_format_clusters_validate_clustering_solutions(valid_clustering_solutions,
                                                              valid_clustering_solutions_with_noise,
                                                              invalid_clustering_solutions):
    BinaryFormatClusters._validate_clustering_solutions(valid_clustering_solutions)
    BinaryFormatClusters._validate_clustering_solutions(valid_clustering_solutions_with_noise)
    with pytest.raises(AssertionError):
        BinaryFormatClusters._validate_clustering_solutions(invalid_clustering_solutions)
    with pytest.raises(AssertionError):
        BinaryFormatClusters._validate_clustering_solutions(tuple(valid_clustering_solutions))
    with pytest.raises(AssertionError):
        BinaryFormatClusters._validate_clustering_solutions([])
    with pytest.raises(AssertionError):
        BinaryFormatClusters._validate_clustering_solutions(valid_clustering_solutions + [[0, 0, 0]])


def _sort_clusters(clusters: np.ndarray):
    return sorted([np.where(clusters == i)[0] for i in np.unique(clusters) if i >= 0], key=lambda x: x.sum())


def test_clicom_majority_voter(valid_clustering_solutions):
    truth = _sort_clusters(np.array([0, 0, 0, 0, 1, 1, 1, 1]))
    n_features = 8
    clusterer = CLICOM.__new__(CLICOM)
    clusterer.n_features = n_features
    clusterer.min_cluster_size = 1
    clusterer.clustering_solutions = BinaryFormatClusters(valid_clustering_solutions)
    cluster_blocks = [frozenset({2, 4, 7, 8}), frozenset({0, 1, 3, 5, 6})]

    res = _sort_clusters(clusterer.majority_voter(cluster_blocks))
    assert np.all([np.all(arr_res == arr_truth) for arr_res, arr_truth in zip(res, truth)])


def test_clicom_majority_voter_no_solutions():
    truth = np.array([-1, -1, -1, -1, -1])
    n_features = 5
    clusterer = CLICOM.__new__(CLICOM)
    clusterer.n_features = n_features
    clusterer.min_cluster_size = 1

    assert np.all(clusterer.majority_voter([]) == truth)


def test_clicom_cliques_to_blocks():
    truth = sorted([frozenset({0, 1, 3, 5, 6}), frozenset({2, 4, 7, 8})], key=len)
    clusterer = CLICOM.__new__(CLICOM)
    clusterer.cluster_unclustered_features = True
    clusterer.clique_set = {frozenset({2, 4, 7, 8}),
                            frozenset({1, 7}),
                            frozenset({1, 3, 6}),
                            frozenset({0, 3, 5, 6})}
    clusterer.binary_mat = np.array(
        [[0, 0, 0, 1, 0, 1, 1, 0, 0],
         [0, 0, 0, 1, 0, 0, 1, 1, 0],
         [0, 0, 0, 0, 1, 0, 0, 1, 1],
         [1, 1, 0, 0, 0, 1, 1, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 1, 1],
         [1, 0, 0, 1, 0, 0, 1, 0, 0],
         [1, 1, 0, 1, 0, 1, 0, 0, 0],
         [0, 1, 1, 0, 1, 0, 0, 0, 1],
         [0, 0, 1, 0, 1, 0, 0, 1, 0]])
    clusterer.adj_mat = clusterer.binary_mat

    assert sorted(clusterer.cliques_to_blocks(), key=len) == truth


def test_clicom_get_cluster_similarity_matrix(valid_clustering_solutions):
    truth = np.array([
        [0, 48, 0, 144, 0, 216, 72, 0, 0],
        [48, 0, 32, 108, 54, 48, 168, 84, 24],
        [0, 32, 0, 0, 162, 0, 0, 132, 192],
        [144, 108, 0, 0, 9, 144, 144, 18, 0],
        [0, 54, 162, 9, 0, 0, 18, 144, 162],
        [216, 48, 0, 144, 0, 0, 72, 0, 0],
        [72, 168, 0, 144, 18, 72, 0, 36, 0],
        [0, 84, 132, 18, 144, 0, 36, 0, 108],
        [0, 24, 192, 0, 162, 0, 0, 108, 0], ]) / 216

    clusterer = CLICOM.__new__(CLICOM)
    clusterer.parallel_backend = 'loky'
    clusterer.clustering_solutions = BinaryFormatClusters(valid_clustering_solutions)
    clusterer.n_features = 8

    res = clusterer.get_cluster_similarity_matrix()
    assert np.isclose(res, truth).all()


def test_clicom_get_evidence_accumulation_matrix(valid_clustering_solutions):
    n_solutions = 3
    truth = np.array([
        [0, 3, 1, 1, 0, 0, 0, 0],
        [3, 0, 1, 1, 0, 0, 0, 0],
        [1, 1, 0, 3, 1, 0, 0, 0],
        [1, 1, 3, 0, 1, 0, 0, 0],
        [0, 0, 1, 1, 0, 2, 1, 1],
        [0, 0, 0, 0, 2, 0, 2, 2],
        [0, 0, 0, 0, 1, 2, 0, 3],
        [0, 0, 0, 0, 1, 2, 3, 0]]) / n_solutions

    clusterer = CLICOM.__new__(CLICOM)
    clusterer.clustering_solutions = BinaryFormatClusters(valid_clustering_solutions)
    clusterer.n_features = 8
    assert np.all(clusterer.get_evidence_accumulation_matrix() == truth)


def test_clicom_clusters_to_labels():
    n_features = 9
    clusters = [{0, 1, 2, 3}, {4, 5, 6, 7}]
    binary_format_clusters = BinaryFormatClusters.__new__(BinaryFormatClusters)
    binary_format_clusters.n_features = n_features
    truth = np.array([[0, 0, 0, 0, 1, 1, 1, 1, -1]])

    clusterer = CLICOM.__new__(CLICOM)
    clusterer.clustering_solutions = binary_format_clusters
    assert np.all(clusterer.clusters_to_labels(clusters) == truth)


@pytest.mark.parametrize("cluster_unclustered_features,allowed_overlap, min_cluster_size,expected_clusters", [
    (True, 1 / 3, 1, [{0, 1, 2, 3, 12}, {6, 7, 8, 9}, {10, 11}, {4, 5}]),
    (False, 1 / 3, 1, [{0, 1, 2, 3}, {6, 7, 8, 9}, {10, 11}, {4, 5}]),
    (True, 0, 1, [{0, 1, 2, 3, 12}, {6, 7, 8, 9, 5}, {10, 11}]),
    (False, 0, 1, [{0, 1, 2, 3}, {6, 7, 8, 9}, {10, 11}]),
    (False, 0, 3, [{0, 1, 2, 3}, {6, 7, 8, 9}])
])
def test_clicom_cliques_to_clusters(monkeypatch, cluster_unclustered_features, allowed_overlap, min_cluster_size,
                                    expected_clusters):
    threshold = 0.8

    def mock_feature_cluster_similarity(self, feature, cluster):
        if feature == 4:
            best_cluster = None  # feature 4 should remain unclustered
        elif feature == 5:
            best_cluster = {6, 7, 8, 9}
        elif feature == 12:
            best_cluster = {0, 1, 2, 3}
        else:
            raise ValueError("an unexpected feature remained unclustered. ")

        if cluster == best_cluster:
            return threshold
        else:
            return threshold - 0.01
            # return a similarity score that's lower than the best

    monkeypatch.setattr(CLICOM, 'feature_cluster_similarity', mock_feature_cluster_similarity)

    n_features = 13
    cliques_lst = [{0, 1, 2, 3}, {4, 5, 6}, {2, 3, 4}, {6, 7, 8, 9}, {10, 11}, {2, 3, 12}]
    cliques_set = {frozenset(this_clique) for this_clique in cliques_lst}
    clusterer = CLICOM.__new__(CLICOM)
    clusterer.binary_mat = np.zeros((n_features, 1))
    clusterer.cluster_unclustered_features = cluster_unclustered_features
    clusterer.min_cluster_size = min_cluster_size
    clusterer.clique_set = cliques_set
    if cluster_unclustered_features:
        clusterer.cluster_wise_cliques = False  # we are not in cluster-wise cliques mode
        clusterer.threshold = threshold
    res = clusterer.cliques_to_clusters(allowed_overlap)

    assert isinstance(res, list)
    for clstr in res:
        assert isinstance(clstr, set)

    assert {frozenset(clstr) for clstr in res} == {frozenset(clstr) for clstr in expected_clusters}


@pytest.mark.parametrize("feature,cluster,expected", [
    (0, {4, 5, 6}, np.mean([0, 216, 72])),
    (0, {7, 8, 1}, np.mean([48, 0, 0])),
    (3, {4, 5, 6}, np.mean([9, 144, 144])),
    (3, {7, 8, 1}, np.mean([108, 18, 0])),
])
def test_clicom_feature_cluster_similarity(feature, cluster, expected):
    adj_mat = np.array([
        [0, 48, 0, 144, 0, 216, 72, 0, 0],
        [48, 0, 32, 108, 54, 48, 168, 84, 24],
        [0, 32, 0, 0, 162, 0, 0, 132, 192],
        [144, 108, 0, 0, 9, 144, 144, 18, 0],
        [0, 54, 162, 9, 0, 0, 18, 144, 162],
        [216, 48, 0, 144, 0, 0, 72, 0, 0],
        [72, 168, 0, 144, 18, 72, 0, 36, 0],
        [0, 84, 132, 18, 144, 0, 36, 0, 108],
        [0, 24, 192, 0, 162, 0, 0, 108, 0], ])
    clusterer = CLICOM.__new__(CLICOM)
    clusterer.adj_mat = adj_mat
    clusterer.cluster_wise_cliques = False
    assert clusterer.feature_cluster_similarity(feature, cluster) == expected

    clusterer.cluster_wise_cliques = True
    with pytest.raises(AssertionError):
        _ = clusterer.feature_cluster_similarity(feature, cluster)


@pytest.mark.parametrize("a,b,expected", [
    (0, 1, 4 / 18),
    (0, 3, 16 / 24),
    (0, 5, 1),
    (3, 5, 16 / 24)])
def test_clicom_inter_cluster_similarity(a, b, expected, valid_clustering_solutions):
    clusterer = CLICOM(BinaryFormatClusters(valid_clustering_solutions), 0.5)
    assert clusterer.inter_cluster_similarity(a, b, clusterer.clustering_solutions.cluster_sets,
                                              clusterer.clustering_solutions.n_features,
                                              clusterer.clustering_solutions.n_solutions) == expected


@pytest.mark.parametrize("clique,expected", [
    (frozenset({0, 1}), 8 / 18),
    (frozenset({0, 5}), 2),
    (frozenset({0, 3, 5}), (32 / 24) + 1)])
def test_clicom_cumulative_cluster_similarity(valid_clustering_solutions, clique, expected):
    clusterer = CLICOM.__new__(CLICOM)
    clusterer.adj_mat = np.array([
        [0, 48, 0, 144, 0, 216, 72, 0, 0],
        [48, 0, 32, 108, 54, 48, 168, 84, 24],
        [0, 32, 0, 0, 162, 0, 0, 132, 192],
        [144, 108, 0, 0, 9, 144, 144, 18, 0],
        [0, 54, 162, 9, 0, 0, 18, 144, 162],
        [216, 48, 0, 144, 0, 0, 72, 0, 0],
        [72, 168, 0, 144, 18, 72, 0, 36, 0],
        [0, 84, 132, 18, 144, 0, 36, 0, 108],
        [0, 24, 192, 0, 162, 0, 0, 108, 0], ]) / 216

    assert clusterer.cumulative_cluster_similarity(clique) == expected
