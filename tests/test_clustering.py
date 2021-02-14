import pytest
from rnalysis.utils.clustering import *
from rnalysis.utils.io import load_csv
from collections import namedtuple


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


def test_parse_n_clusters(monkeypatch):
    data = pd.read_csv('tests/test_files/counted.csv', index_col=0)
    monkeypatch.setattr(ClusteringRunnerWithNClusters, "silhouette", lambda self: 6)
    monkeypatch.setattr(ClusteringRunnerWithNClusters, "gap_statistic", lambda self: 7)

    runner = ClusteringRunnerWithNClusters(data, False, 'silhouette', 20)
    assert runner.n_clusters == [6]
    runner = ClusteringRunnerWithNClusters(data, False, 'gap', 20)
    assert runner.n_clusters == [7]

    runner = ClusteringRunnerWithNClusters(data, False, [7, 2, 5], 20)
    assert runner.n_clusters == [7, 2, 5]

    runner = ClusteringRunnerWithNClusters(data, False, [10], 20)
    assert runner.n_clusters == [10]

    runner = ClusteringRunnerWithNClusters(data, False, 10, 20)
    assert runner.n_clusters == [10]

    runner = ClusteringRunnerWithNClusters(data, False, range(2, 9), 20)
    assert runner.n_clusters == list(range(2, 9))

    with pytest.raises(AssertionError):
        runner = ClusteringRunnerWithNClusters(data, False, [5, 2, '3'], 20)
    with pytest.raises(AssertionError):
        runner = ClusteringRunnerWithNClusters(data, False, '17', 20)
    with pytest.raises(AssertionError):
        runner = ClusteringRunnerWithNClusters(data, False, 1, 20)
    with pytest.raises(AssertionError):
        runner = ClusteringRunnerWithNClusters(data, False, [3, 5, 1], 20)
    with pytest.raises(AssertionError):
        runner = ClusteringRunnerWithNClusters(data, False, [3, 23, 2], 20)


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


def test_clicomrunner_find_valid_clustering_setups():
    truth = [(KMeansRunner, dict(power_transform=True, n_clusters=2, n_init=5)),
             (KMeansRunner, dict(power_transform=True, n_clusters=3, n_init=5)),
             (KMeansRunner, dict(power_transform=True, n_clusters=5, n_init=5)),
             (KMeansRunner, dict(power_transform=False, n_clusters=2, n_init=5)),
             (KMeansRunner, dict(power_transform=False, n_clusters=3, n_init=5)),
             (KMeansRunner, dict(power_transform=False, n_clusters=5, n_init=5)),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=5, metric='jackknife', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=7, metric='jackknife', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=10, metric='jackknife', linkage='average')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=5, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=5, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=5, metric='jackknife', linkage='average')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=7, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=7, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=7, metric='jackknife', linkage='average')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=10, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=10, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=False, n_clusters=10, metric='jackknife', linkage='average')),

             (HierarchicalRunner,
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=0.5)),
             (HierarchicalRunner,
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=1)),
             (HierarchicalRunner,
              dict(power_transform=False, n_clusters=None, metric='euclidean', distance_threshold=0.5)),
             (HierarchicalRunner,
              dict(power_transform=False, n_clusters=None, metric='euclidean', distance_threshold=1))]

    runner = CLICOMRunner(pd.DataFrame(np.zeros((10, 10))), [True, False], 0.5,
                          dict(method='kmeans', n_clusters=[2, 3, 5], n_init=5),
                          dict(method='hierarchical', n_clusters=[5, 7, 10], metric=['euclidean', 'jackknife', 'l1'],
                               linkage=['ward', 'average', 'doesnt exist']),
                          dict(method='hierarchical', n_clusters=None, distance_threshold=[0.5, 1], metric='euclidean'))
    setups = runner.find_valid_clustering_setups()
    assert setups == truth

    truth = [(KMeansRunner, dict(power_transform=True, n_clusters=2, n_init=5)),
             (KMeansRunner, dict(power_transform=True, n_clusters=3, n_init=5)),
             (KMeansRunner, dict(power_transform=True, n_clusters=5, n_init=5)),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=5, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=5, metric='jackknife', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=7, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=7, metric='jackknife', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='ward')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=10, metric='euclidean', linkage='average')),
             (HierarchicalRunner, dict(power_transform=True, n_clusters=10, metric='jackknife', linkage='average')),
             (HierarchicalRunner,
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=0.5)),
             (HierarchicalRunner,
              dict(power_transform=True, n_clusters=None, metric='euclidean', distance_threshold=1))]

    runner = CLICOMRunner(pd.DataFrame(np.zeros((10, 10))), True, 0.5,
                          dict(method='kmeans', n_clusters=[2, 3, 5], n_init=5),
                          dict(method='hierarchical', n_clusters=[5, 7, 10], metric=['euclidean', 'jackknife', 'l1'],
                               linkage=['ward', 'average','doesnt exist']),
                          dict(method='hierarchical', n_clusters=None, distance_threshold=[0.5, 1], metric='euclidean'))
    setups = runner.find_valid_clustering_setups()
    assert setups == truth


def test_find_cliques(monkeypatch):
    binary_adj_mat = pd.read_csv('tests/test_files/clicom_binary_adj_matrix.csv', index_col=0).values
    truth = {frozenset({2, 4, 7, 8}), frozenset({1, 7}), frozenset({1, 3, 6}), frozenset({0, 3, 5, 6})}
    monkeypatch.setattr(CLICOM, "get_cluster_similarity_matrix", lambda self: binary_adj_mat)
    clicom = CLICOM(BinaryFormatClusters([np.array([[0, 1]])]), 1)
    clicom.find_cliques()
    assert clicom.clique_set == truth


def test_clicom_cluster_wise_result():
    clustering_collection = [np.array([[1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0, 0],
                                       [0, 0, 0, 0, 0, 1, 1, 1]]),
                             np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]),
                             np.array([[1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0],
                                       [0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1]])]
    bin_format = BinaryFormatClusters(clustering_collection)

    truth = {frozenset({0, 1, 2, 3}), frozenset({4, 5, 6, 7})}
    clusterer = CLICOM(bin_format, 1 / 3)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth == results

    truth_threshold_0 = {frozenset({0, 1, 2, 3, 4, 5, 6, 7})}
    clusterer = CLICOM(bin_format, 0)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth_threshold_0 == results

    truth_threshold_max = set()
    clusterer = CLICOM(bin_format, 1.1)
    clusterer.run()
    results = {frozenset(np.where(clusterer.labels_ == i)[0]) for i in np.unique(clusterer.labels_) if i >= 0}
    assert truth_threshold_max == results


def test_clicom_majority_voter():
    assert False


def test_clicom_cliques_to_clusters():
    assert False


def test_clicom_cliques_to_blocks():
    assert False


def test_clicom_get_cluster_similarity_matrix():
    assert False


def test_clicom_get_evidence_accumulation_matrix():
    assert False


def test_clicom_clusters_to_labels():
    assert False
