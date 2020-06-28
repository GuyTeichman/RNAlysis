import numpy as np
from sklearn.preprocessing import PowerTransformer, StandardScaler


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
