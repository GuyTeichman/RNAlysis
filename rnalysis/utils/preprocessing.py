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
    return StandardScaler().fit_transform(PowerTransformer(method='box-cox').fit_transform(data + 1))


def standardize(data: np.ndarray):
    """

    :param data:
    :type data:
    :return:
    :rtype:
    """
    return StandardScaler().fit_transform(data)


def intersection_nonempty(*objs: Union[list, set, tuple]):
    return set.intersection(*[set(item) for item in objs if len(item) > 0])
