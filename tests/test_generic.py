import pytest
from rnalysis.utils.generic import *
import numpy as np


def test_intersection_nonempty():
    assert intersection_nonempty({1, 2, 3}, {2, 5, 7, 1}) == {1, 2}
    assert intersection_nonempty({1, 3, 7}, set()) == {1, 3, 7}
    assert intersection_nonempty({1, 2, 4, 5}, set(), {1, 3, 4, 7}) == {1, 4}


def test_standardize():
    np.random.seed(42)
    data = np.random.randint(-200, 100000, (100, 5))
    res = standardize(data)
    assert res.shape == data.shape
    assert np.isclose(res.mean(axis=0), 0).all()
    assert np.isclose(res.std(axis=0), 1).all()


def test_standard_box_cox():
    np.random.seed(42)
    data = np.random.randint(-200, 100000, (100, 5))
    res = standard_box_cox(data)
    assert res.shape == data.shape
    assert np.isclose(res.mean(axis=0), 0).all()
    assert np.isclose(res.std(axis=0), 1).all()
    assert not np.isclose(res, standardize(data)).all()

    data_df = pd.DataFrame(data, index=[f'ind{i}' for i in range(100)], columns=[f'col{j}' for j in range(5)])
    res_df = standard_box_cox(data_df)
    assert isinstance(res_df, pd.DataFrame)
    assert res_df.shape == data_df.shape
    assert np.all(res_df.index == data_df.index)
    assert np.all(res_df.columns == data_df.columns)


def test_color_generator():
    gen = color_generator()
    preset_colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink',
                     'tab:gray', 'tab:olive', 'tab:cyan', 'gold', 'maroon', 'mediumslateblue', 'fuchsia',
                     'mediumblue', 'black', 'lawngreen']
    for i in range(150):
        color = next(gen)
        assert color in preset_colors or (isinstance(color, np.ndarray) and len(color) == 3 and
                                          np.max(color) <= 1 and np.min(color) >= 0)
