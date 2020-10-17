import pandas as pd
import numpy as np
from sklearn.preprocessing import PowerTransformer, StandardScaler
from typing import Union


def standard_box_cox(data: np.ndarray):
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


def standardize(data: Union[np.ndarray, pd.DataFrame]):
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
