from pathlib import Path
from typing import Union, Iterable
import pandas as pd
import warnings


def check_is_df_like(inp):
    """
    checks whether an input file is a pandas DataFrame, a string that represent a path of a .csv file, a Path object \
    of a .csv file, or an invalid input.

    :param inp: the input we wish to test
    :return: True if pandas DataFrame, False if a string/Path object that leads to a .csv file.\
    Raises ValueError otherwise.
    """
    if isinstance(inp, pd.DataFrame):
        return True
    elif isinstance(inp, str):
        if inp[-4::] == '.csv':
            return False
    elif isinstance(inp, Path):
        if inp.suffix == '.csv':
            return False
    raise ValueError("The input is neither a pandas DataFrame or a csv file")


def isinstanceinh(obj, parent_class):
    return True if issubclass(obj.__class__, parent_class) else False


def isinstanceiter(iterable: Iterable, object_class: type):
    assert isiterable(iterable), f"Object of type {type(iterable)} is not iterable."
    return all([isinstance(i, object_class) for i in iterable])


def isiterable(obj):
    try:
        _ = iter(obj)
    except TypeError:
        return False
    else:
        return True


def validate_biotype_table(ref_df: pd.DataFrame):
    """
    Assert legality of Biotype Reference Table, and rename column names to standard names ('gene' and 'biotype').
    :param ref_df: the loaded Biotype Reference Table
    :type ref_df: pandas DataFrame

    """
    assert ref_df.shape[
               1] == 2, f"Invalid number of columns in Biotype Reference Table: found {ref_df.shape[1]} columns instead of 2!"
    assert ref_df.shape[
               0] >= 2, f"Biotype Reference Table must have at least two rows, found only  {ref_df.shape[0]}!"
    ref_df.rename(columns={ref_df.columns[0]: 'gene', ref_df.columns[1]: 'biotype'}, inplace=True)


def validate_attr_table(ref_df: pd.DataFrame):
    """
    Assert legality of Attribute Reference Table, and renames the first column to standard name ('gene').
    :param ref_df:
    :type ref_df: pandas DataFrame

    """
    assert ref_df.shape[
               1] >= 2, f"Attribute Reference Table must have at least two columns, found only  {ref_df.shape[1]}!"
    assert ref_df.shape[
               0] >= 2, f"Attribute Reference Table must have at least two rows, found only  {ref_df.shape[0]}!"
    ref_df.rename(columns={ref_df.columns[0]: 'gene'}, inplace=True)


def validate_uniprot_dataset_name(dataset_dict: dict, *names: str):
    for name in names:
        assert name in dataset_dict, f"Dataset '{name}' is not a valid Uniprot Dataset for mapping gene names/IDs. " \
                                     f"Valid Uniprot Datasets are: {', '.join(dataset_dict.keys())}."


def validate_threshold(threshold: float = 1):
    """
    Assertions for functions that filter normalized values by some threshold, \
    or are meant to be used on normalized values.

    :param threshold: optional. A threshold value for filter_low_rpm to be asserted.

    """
    assert isinstance(threshold, (float, int)), "Threshold must be a number!"
    assert threshold >= 0, "Threshold must be zero or larger!"


def validate_is_normalized(filter_obj_fname: Union[str, Path]):
    if '_norm' not in str(filter_obj_fname):
        warnings.warn("This function is meant for normalized values, and your values may not be normalized. ")


def validate_clustering_parameters(metric: str, linkage: str = None):
    legal_metrics = {'euclidean', 'spearman', 'pearson', 'manhattan', 'cosine', 'ys1', 'yr1', 'jackknife'}
    legal_linkages = {'single', 'average', 'complete', 'ward'}
    assert isinstance(metric, str), f"'metric' must be a string. Instead got '{type(metric)}'."
    metric = metric.lower()
    assert metric in legal_metrics, f"'metric' must be one of {legal_metrics}. Instead got '{metric}'."

    if linkage is not None:
        assert isinstance(linkage, str), f"'linkage' must be a string. Instead got '{type(linkage)}'."
        linkage = linkage.lower()
        assert linkage in legal_linkages, f"'linkage' must be in {legal_linkages}. Instead got '{linkage}'."
        return metric, linkage
    return metric


def validate_hdbscan_parameters(min_cluster_size: int, metric: str, cluster_selection_method: str, n_features: int):
    assert isinstance(min_cluster_size, int) and min_cluster_size >= 2, \
        f"'min_cluster_size' must be an integer >=2. Instead got {min_cluster_size}, type={type(min_cluster_size)}."
    assert min_cluster_size <= n_features, \
        "'min_cluster_size' cannot be larger than the number of features in the CountFilter object. "
    assert isinstance(cluster_selection_method, str), \
        f"'cluster_selection_method' must be a string. Instead got {type(cluster_selection_method)}."
    assert isinstance(metric, str), f"'metric' must be a string. Instead got {type(metric)}."
