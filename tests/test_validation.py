import pytest
from rnalysis.utils.validation import *


class DummyClass:
    def __init__(self):
        pass


def test_is_df_dataframe():
    my_df = pd.DataFrame()
    assert (check_is_df_like(my_df))


def test_is_df_str():
    correct_path = "myfile.csv"
    correct_path_2 = r"myfolder\anotherfile.csv"
    assert (not check_is_df_like(correct_path))
    assert (not check_is_df_like(correct_path_2))


def test_is_df_pathobj():
    correct_pathobj = Path("all_feature_96_new.csv")
    assert (not check_is_df_like(correct_pathobj))


def test_is_df_str_notcsv():
    incorrect_pth = "myfile.xlsx"
    incorrect_pth2 = "all_feature_96_new"
    with pytest.raises(ValueError):
        check_is_df_like(incorrect_pth)
    with pytest.raises(ValueError):
        check_is_df_like(incorrect_pth2)


def test_is_df_pathobj_notcsv():
    incorrect_pathobj = Path("test_general.py")
    with pytest.raises(ValueError):
        check_is_df_like(incorrect_pathobj)


def test_is_df_invalid_type():
    invalid_type = 67
    with pytest.raises(ValueError):
        check_is_df_like(invalid_type)


def test_is_iterable():
    assert isiterable(range(10))
    assert isiterable([])
    assert isiterable('fortytwo')
    assert not isiterable(42)
    assert not isiterable(3.14)
    assert not isiterable(True)
    assert isiterable((i for i in range(10)))
    assert isiterable({})
    assert isiterable({1, 2, 4})


def test_is_iterable_not_emptying_generator():
    gen = (i for i in range(10) if i != 7)
    assert isiterable(gen)
    assert list(gen) == [i for i in range(10) if i != 7]


def test_validate_uniprot_dataset_name():
    validate_uniprot_dataset_name({'one': 1, 'two': 2, 'three': 3}, 'one', 'two', 'three')
    with pytest.raises(AssertionError):
        validate_uniprot_dataset_name({'one': 1, 'two': 2, 'three': 3}, 'one', 2, 'three')
    with pytest.raises(AssertionError):
        validate_uniprot_dataset_name({'one': 1, 'two': 2, 'three': 3}, 'one', 'Two', 'three')


def test_validate_hdbscan_parameters():
    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(1, 'euclidean', 'eom', 1)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(2.0, 'euclidean', 'eom', 13)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(14, 'euclidean', 'eom', 13)

    validate_hdbscan_parameters(13, 'euclidean', 'eom', 13)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(2, 'euclidean', 5, 13)

    with pytest.raises(AssertionError):
        validate_hdbscan_parameters(2, 5, 'EOM', 13)


@pytest.mark.parametrize("test_input,expected_type,expected", [
    (["3", "hello world", "", ], str, True),
    ({5, 3, 7, 8}, int, True),
    ((5, 3, '7', 8), int, False),
    ([[], [1], [1, 2]], int, False),
    ([], DummyClass, True),
    ({'one': DummyClass(), 'two': DummyClass()}.values(), DummyClass, True)
])
def test_isinstanceiter(test_input, expected_type, expected):
    assert isinstanceiter(test_input, expected_type) == expected


class TestClass:
    def cls_method(self):
        pass

    def other_cls_method(self):
        pass


@pytest.mark.parametrize('mthd,cls,truth',
                         [(str.lower, str, True),
                          (str.lower, int, False),
                          (lambda x: x + 1, str, False),
                          (TestClass.cls_method, TestClass, True),
                          (TestClass().cls_method, TestClass, True),
                          (str.lower, TestClass, False)])
def test_is_method_of_class(mthd, cls, truth):
    assert is_method_of_class(mthd, cls) == truth


def test_validate_clustering_parameters():
    assert False


def test_validate_attr_table():
    assert False


def test_validate_biotype_table():
    assert False
