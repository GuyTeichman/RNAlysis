import numpy as np
from sklearn_extra.cluster import KMedoids


class KMedoidsIter:
    def __init__(self, n_clusters: int, metric: str = 'euclidean', init: str = 'k-medoids++', max_iter: int = 300,
                 random_state: int = None, n_init: int = 10):
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
