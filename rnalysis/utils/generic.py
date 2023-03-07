import inspect
import itertools
import math
import typing
import warnings
from functools import lru_cache
from typing import Union, Callable, Tuple

import joblib
import matplotlib.collections
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import figure
from matplotlib.widgets import MultiCursor
from scipy.special import comb
from sklearn.preprocessing import PowerTransformer, StandardScaler
from tqdm.auto import tqdm

try:
    import numba
except ImportError:  # pragma: no cover
    warnings.warn("RNAlysis can perform faster when package 'numba' is installed. \n"
                  "If you want to improve the performance of slow operations on RNAlysis, "
                  "please install package 'numba'. ")


    class numba:  # pragma: no cover
        @staticmethod
        def jit(*args, **kwargs):
            return lambda f: f

        @staticmethod
        def njit(*args, **kwargs):
            return lambda f: f


class ProgressParallel(joblib.Parallel):
    # tqdm progress bar for parallel tasks based on:
    # https://stackoverflow.com/questions/37804279/how-can-we-use-tqdm-in-a-parallel-execution-with-joblib/50925708
    # answer by 'user394430'
    def __init__(self, use_tqdm=True, total=None, desc: str = '', unit: str = 'it', backend='loky', *args, **kwargs):
        self._use_tqdm = use_tqdm
        self._total = total
        self._desc = desc
        self._unit = unit
        kwargs['n_jobs'] = -2
        super().__init__(*args, backend=backend, **kwargs, verbose=100)

    def __call__(self, *args, **kwargs):
        fmt = '{desc}: {percentage:3.0f}%|{bar}| [{elapsed}<{remaining}, {rate_fmt}{postfix}]'
        with tqdm(disable=not self._use_tqdm, total=self._total, desc=self._desc, unit=self._unit,
                  bar_format=fmt) as self._pbar:
            return joblib.Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        if self._total is None:
            self._pbar.total = self.n_dispatched_tasks
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()


def standard_box_cox(data: Union[np.ndarray, pd.DataFrame]) -> Union[np.ndarray, pd.DataFrame]:
    """

    :param data:
    :type data:
    :return:
    :rtype:
    """
    res_array = StandardScaler().fit_transform(PowerTransformer(method='box-cox').fit_transform(data + 1))
    if isinstance(data, pd.DataFrame):
        return pd.DataFrame(res_array, index=data.index, columns=data.columns)
    return res_array


def shift_to_baseline(data: Union[np.ndarray, pd.DataFrame], baseline: float = 0) -> Union[np.ndarray, pd.DataFrame]:
    """

    :param data:
    :type data:
    :param baseline:
    :type baseline:
    :return:
    :rtype:
    """
    min_val = data.min()
    while isinstance(min_val, (pd.DataFrame, np.ndarray, pd.Series)):
        min_val = min_val.min()
    diff = baseline - min_val
    return data + diff


def standardize(data: Union[np.ndarray, pd.DataFrame]) -> Union[np.ndarray, pd.DataFrame]:
    """

    :param data:
    :type data:
    :return:
    :rtype:
    """
    res_array = StandardScaler().fit_transform(data)
    if isinstance(data, pd.DataFrame):
        return pd.DataFrame(res_array, index=data.index, columns=data.columns)
    return res_array


def intersection_nonempty(*objs: Union[list, set, tuple]):
    return set.intersection(*[set(item) for item in objs if len(item) > 0])


def color_generator():
    """
    A generator that yields distinct colors up to a certain limit, and then yields randomized RGB values.

    :return: a color name string (like 'black', \
    or a numpy.ndarray of size (3,) containing three random values each between 0 and 1.

    """
    preset_colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink',
                     'tab:gray', 'tab:olive', 'tab:cyan', 'gold', 'maroon', 'mediumslateblue', 'fuchsia',
                     'lawngreen', 'moccasin', 'thistle']
    for color in preset_colors:
        yield color
    while True:
        yield np.random.random(3)


@lru_cache(maxsize=2 ** 16)
def combination(a: int, b: int) -> int:
    return int(comb(a, b))


class SetWithMajorityVote(set):
    def majority_vote_intersection(self, *others, majority_threshold: float = 0.5):
        counts = {}
        n_sets = len(others) + 1
        for set_obj in itertools.chain([self], others):
            for obj in set_obj:
                counts[obj] = counts.get(obj, 0) + (1 / n_sets)
        result = {key for key, val in counts.items() if val >= majority_threshold}
        return result


def get_method_readable_name(method: Union[str, Callable], obj: object = None):
    try:
        if isinstance(method, str):
            func = getattr(obj, method)
        else:
            func = method
        if hasattr(func, 'readable_name'):
            name = func.readable_name
        else:
            name = func.__name__
        return name
    except AttributeError:
        return 'NAME NOT FOUND'


def get_method_signature(method: Union[str, Callable], obj: object = None):
    try:
        if isinstance(method, str):
            func = getattr(obj, method)
        else:
            func = method
        signature = inspect.signature(func)
        return signature.parameters
    except AttributeError:
        return {}


def despine(ax):
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)


def mix_colors(*colors: Tuple[float, float, float]):
    n_colors = len(colors)
    if n_colors == 1:
        return colors[0]

    colors = np.array(colors)
    multiplier = (1 / n_colors) + (1.6 / (2 ** n_colors)) / n_colors

    mix_color = multiplier * np.sum(colors, axis=0)
    mix_color = np.min([mix_color, [1.0, 1.0, 1.0]], 0)
    return mix_color


def sum_intervals_inclusive(intervals: typing.List[typing.Tuple[int, int]]) -> int:
    sorted_intervals = sorted(intervals)
    total = 0
    prev_interval = (-np.inf, -np.inf)
    for i in range(len(sorted_intervals)):
        this_interval = sorted_intervals[i]
        if this_interval[0] <= prev_interval[1]:
            this_interval = (prev_interval[1] + 1, this_interval[1])
        total += this_interval[1] - this_interval[0] + 1
        prev_interval = this_interval

    return total


def bic_score(X: np.ndarray, labels: np.ndarray):
    """
      BIC score for the goodness of fit of clusters.
      This Python function is translated from the Golang implementation by the author of the paper.
      The original code is available here:
      https://github.com/bobhancock/goxmeans/blob/a78e909e374c6f97ddd04a239658c7c5b7365e5c/km.go#L778
      """

    n_points = len(labels)
    n_clusters = len(set(labels))
    n_dimensions = X.shape[1]

    n_parameters = (n_clusters - 1) + (n_dimensions * n_clusters) + 1

    loglikelihood = 0
    for label_name in set(labels):
        X_cluster = X[labels == label_name]
        n_points_cluster = len(X_cluster)
        centroid = np.mean(X_cluster, axis=0)
        variance = np.sum((X_cluster - centroid) ** 2) / (len(X_cluster) - 1)
        loglikelihood += \
            n_points_cluster * np.log(n_points_cluster) \
            - n_points_cluster * np.log(n_points) \
            - n_points_cluster * n_dimensions / 2 * np.log(2 * math.pi * variance) \
            - (n_points_cluster - 1) / 2

    bic = loglikelihood - (n_parameters / 2) * np.log(n_points)

    return bic


def format_time(seconds: float):
    seconds = int(seconds)
    n_minutes = seconds // 60
    n_seconds = seconds % 60
    return f'{n_minutes:02d}:{n_seconds:02d}'


def sanitize_variable_name(name: str) -> str:
    """
    Sanitize a string to turn it into a legal variable name in R/Python.
    :param name: name to sanitize
    :type name: str
    :return: sanitizeed name
    :rtype: str
    """
    new_name = name.rstrip().replace(' ', '_')
    if new_name[0].isdigit():
        new_name = 'var_' + new_name

    if not new_name.isalnum():
        new_name = ''.join([char if char.isalnum() else '_' for char in new_name])

    return new_name


class InteractiveScatterFigure(figure.Figure):
    def __init__(self, labels: typing.List[str], annotation_fontsize: float = 10, show_cursor: bool = True, *args,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.ax = self.add_subplot()
        self.labels = labels
        self.annotation_fontsize = annotation_fontsize
        self.is_labeled: typing.Dict[int, plt.Line2D] = {}
        self.canvas.mpl_connect('pick_event', self.on_pick)
        if show_cursor:
            self.cursor = MultiCursor(self.canvas, self.axes, color='k', lw=0.5, horizOn=True, vertOn=True)
            self.canvas.mpl_connect('axes_leave_event', self.on_exit)
            self.canvas.mpl_connect('motion_notify_event', self.on_move)

    def on_exit(self, event):
        self.cursor.clear(event)
        self.canvas.draw()

    def on_move(self, event):
        self.canvas.draw()

    def on_pick(self, event):
        for this_ind in event.ind:
            if isinstance(event.artist, matplotlib.collections.PathCollection):
                data = event.artist.get_offsets().data
                xdata, ydata = data[:, 0], data[:, 1]
            elif isinstance(event.artist, plt.Line2D):
                thisline = event.artist
                xdata = thisline.get_xdata()
                ydata = thisline.get_ydata()
            else:
                return

            if this_ind in self.is_labeled:
                ann = self.is_labeled.pop(this_ind)
                ann.remove()
            else:
                self.is_labeled[this_ind] = plt.annotate(self.labels[this_ind],
                                                         (np.take(xdata, this_ind), np.take(ydata, this_ind)),
                                                         xytext=(3, 3), textcoords='offset points',
                                                         fontsize=self.annotation_fontsize)

            self.canvas.draw()
