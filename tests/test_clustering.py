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
