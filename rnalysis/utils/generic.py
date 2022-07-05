import itertools
import inspect
from functools import lru_cache
from typing import Union, Callable

import numpy as np
import pandas as pd
from joblib import Parallel
from scipy.special import comb
from sklearn.preprocessing import PowerTransformer, StandardScaler
from tqdm.auto import tqdm


class ProgressParallel(Parallel):
    # tqdm progress bar for parallel tasks based on:
    # https://stackoverflow.com/questions/37804279/how-can-we-use-tqdm-in-a-parallel-execution-with-joblib/50925708
    # answer by 'user394430'
    def __init__(self, use_tqdm=True, total=None, desc: str = '', unit: str = 'it', *args, **kwargs):
        self._use_tqdm = use_tqdm
        self._total = total
        self._desc = desc
        self._unit = unit
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        fmt = '{desc}: {percentage:3.0f}%|{bar}| [{elapsed}<{remaining}, {rate_fmt}{postfix}]'
        with tqdm(disable=not self._use_tqdm, total=self._total, desc=self._desc, unit=self._unit,
                  bar_format=fmt) as self._pbar:
            return Parallel.__call__(self, *args, **kwargs)

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
                     'mediumblue', 'black', 'lawngreen']
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
