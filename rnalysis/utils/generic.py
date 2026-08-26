import abc
import enum
import inspect
import itertools
import math
import types
import typing
import warnings
from datetime import date, datetime
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple, Union

import joblib
import lazy_loader as lazy
import matplotlib.collections
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import polars.selectors as cs
import yaml
from matplotlib import figure
from matplotlib.widgets import MultiCursor
from scipy.special import comb
from tqdm.auto import tqdm

from rnalysis import FROZEN_ENV, __version__
from rnalysis.exceptions import InternalError, InvalidTypeError, InvalidValueError
from rnalysis.utils.param_typing import POWER_TRANSFORMS

# scikit-learn costs ~1s to import and is only needed once a transform actually runs, so it is loaded
# lazily (SPEC 1 / https://scientific-python.org/specs/spec-0001/). Nothing may touch an attribute of
# `_sklearn` at import time -- doing so imports the package and defeats the whole point.
_sklearn = lazy.load('sklearn')

try:
    import numba
except ImportError:  # pragma: no cover
    warnings.warn(
        "RNAlysis can perform faster when package 'numba' is installed. \n"
        'If you want to improve the performance of slow operations on RNAlysis, '
        "please install package 'numba'. "
    )

    class numba:  # pragma: no cover
        @staticmethod
        def jit(*args, **kwargs):
            return lambda f: f

        @staticmethod
        def njit(*args, **kwargs):
            return lambda f: f


# Cache numba's compiled kernels on disk, so only the first process ever pays JIT compilation and
# every later one (including each Windows `spawn` multiprocessing worker) loads the compiled code
# instead of recompiling it. Pass this as the `cache=` argument of every `numba.jit` in RNAlysis.
#
# It is switched off in the frozen (PyInstaller) app: numba picks the cache directory from the
# *source file* of the jitted function (numba.core.caching.CacheImpl), and if none of its locators
# matches it raises RuntimeError from the decorator -- that is, at import time, which would stop the
# standalone app from starting at all. numba does ship a frozen-aware locator
# (UserWideCacheLocator, which caches into the user-wide cache dir and stamps the cache against
# sys.executable), so this gate can be lifted once it has been verified against a real
# standalone build.
NUMBA_CACHE = not FROZEN_ENV


def readable_name(name: str):
    def decorator(item):
        item.readable_name = name
        return item

    return decorator


def param_readable_names(names: Dict[str, str]):
    """
    Decorator that attaches human-readable display labels for specific parameters of a function.

    The mapping is stored on the function's ``readable_param_names`` attribute (mirroring how
    :func:`readable_name` sets ``readable_name``) and is used by :func:`get_param_readable_name`
    to override the auto-humanized label for parameters that do not humanize well. This is
    display-only: it never changes the actual parameter names used when calling the function.

    :param names: mapping of raw (snake_case) parameter names to their human-readable labels.
    :type names: dict of str to str
    """

    def decorator(item):
        item.readable_param_names = names
        return item

    return decorator


# Domain acronyms / tokens preserved verbatim when auto-humanizing parameter labels.
# Keyed by the lowercase snake_case token; the value is the exact display form.
PARAM_LABEL_ACRONYMS = {
    'go': 'GO',
    'kegg': 'KEGG',
    'pca': 'PCA',
    'fdr': 'FDR',
    'id': 'ID',
    'ids': 'IDs',
    'rna': 'RNA',
    'mrna': 'mRNA',
    'dna': 'DNA',
    'deseq2': 'DESeq2',
    'sam': 'SAM',
    'bam': 'BAM',
    'fastq': 'FASTQ',
    'gtf': 'GTF',
    'gff': 'GFF',
    'utr': 'UTR',
    'cpm': 'CPM',
    'tpm': 'TPM',
    'rpkm': 'RPKM',
    'fpkm': 'FPKM',
    'padj': 'padj',
    'clicom': 'CLICOM',
    'hdbscan': 'HDBSCAN',
    'umap': 'UMAP',
    'ncbi': 'NCBI',
    'url': 'URL',
    'csv': 'CSV',
}

# Whole parameter-name overrides that auto-humanization handles poorly but which apply
# regardless of the containing function (per-function overrides take precedence over these).
PARAM_LABEL_SPECIAL_NAMES = {
    'n_components': 'Number of components',
}


def _humanize_param_name(param_name: str) -> str:
    words = []
    first = True
    for token in param_name.split('_'):
        if token == '':
            continue
        lower = token.lower()
        if lower in PARAM_LABEL_ACRONYMS:
            words.append(PARAM_LABEL_ACRONYMS[lower])
        elif first:
            words.append(token.capitalize())
        else:
            words.append(lower)
        first = False
    if not words:
        return param_name
    return ' '.join(words)


def get_param_readable_name(param_name: str, func: Union[Callable, None] = None) -> str:
    """
    Return a human-readable display label for a function parameter (display-only).

    Resolution order:

    1. a per-function override declared via :func:`param_readable_names`, if ``func`` is given
       and declares a label for ``param_name``;
    2. a global whole-name special case (:data:`PARAM_LABEL_SPECIAL_NAMES`);
    3. auto-humanization: split on ``_``, sentence-case the words and join with spaces, while
       preserving known domain acronyms (:data:`PARAM_LABEL_ACRONYMS`).

    The raw ``param_name`` is never modified for use as an actual keyword argument; this only
    affects what the GUI shows the user.

    :param param_name: the raw (snake_case) parameter name from the function signature.
    :type param_name: str
    :param func: the function the parameter belongs to, used to look up per-function overrides.
    :type func: callable or None
    """
    if func is not None:
        overrides = getattr(func, 'readable_param_names', None)
        if overrides is not None and param_name in overrides:
            return overrides[param_name]
    if param_name in PARAM_LABEL_SPECIAL_NAMES:
        return PARAM_LABEL_SPECIAL_NAMES[param_name]
    return _humanize_param_name(param_name)


class TemporaryMatplotlibBackend:
    def __init__(self, backend):
        self.backend = backend
        self.original_backend = matplotlib.get_backend()

    def __enter__(self):
        matplotlib.use(self.backend, force=True)

    def __exit__(self, exc_type, exc_val, exc_tb):
        matplotlib.use(self.original_backend, force=True)


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
        with tqdm(
            disable=not self._use_tqdm, total=self._total, desc=self._desc, unit=self._unit, bar_format=fmt
        ) as self._pbar:
            return joblib.Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        if self._total is None:
            self._pbar.total = self.n_dispatched_tasks
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()


# The Box-Cox power transform fits one lambda per column via a separate optimization. When the data is
# transposed so that genes are columns (as in the PCA/clustering paths), that is thousands of independent
# optimizations run in a Python loop, which dominates the runtime of every clustering/PCA call. Splitting the
# columns across workers is embarrassingly parallel; the result matches the single-process one up to
# floating-point precision (~1e-15 on the real count-matrix fixture -- transforming a column-slice of a
# C-contiguous matrix rather than the whole matrix can shift numpy's reduction blocking, far below anything
# that affects PCA loadings or cluster assignments). Only pay the process-spawn cost when there are enough
# columns to make it worthwhile.
BOX_COX_PARALLEL_MIN_COLUMNS = 500

# The Box-Cox transform standardizes each column after transforming it, so every value it returns for a
# healthy column is a z-score: the largest one a column of n values can possibly hold is sqrt(n). Anything
# beyond this threshold therefore cannot be a properly standardized column of any real table (it would need
# ~10^12 rows) -- it can only mean the fitted lambda exploded and the transform blew up without quite
# overflowing to infinity. Such a column would silently dominate every principal component, so it is treated
# exactly like a non-finite one.
BOX_COX_MAX_ABS_VALUE = 1e6

#: how many offending gene names to spell out in the instability error before summarizing the rest
_MAX_NAMES_IN_ERROR = 10


def _box_cox_fit_transform(array: np.ndarray) -> np.ndarray:
    # sklearn's PowerTransformer(box-cox) fits + applies a separate lambda to each column independently
    # (and, with the default standardize=True, standardizes each column afterwards).
    return _sklearn.preprocessing.PowerTransformer(method='box-cox').fit_transform(array + 1)


def _parallel_box_cox(array: np.ndarray, backend: str, n_jobs: int) -> np.ndarray:
    # split the columns into contiguous chunks, Box-Cox each chunk in its own worker, then reassemble in the
    # original column order. Because every column is transformed independently, this equals the single-process
    # result exactly.
    n_columns = array.shape[1]
    column_chunks = np.array_split(np.arange(n_columns), min(n_jobs, n_columns))
    transformed_chunks = joblib.Parallel(n_jobs=len(column_chunks), backend=backend)(
        joblib.delayed(_box_cox_fit_transform)(array[:, chunk]) for chunk in column_chunks
    )
    return np.concatenate(transformed_chunks, axis=1)


def _validate_box_cox_stability(box_cox_array: np.ndarray, feature_names: Optional[typing.Sequence[str]]):
    """
    Raise a clear, actionable error if the Box-Cox transform blew up on any column.

    Box-Cox fits one lambda per column. When a column is near-constant at a high magnitude and is measured
    across very few values (the classic case: a reporter transgene such as GFP, sitting at ~8,400 counts in
    a handful of grouped samples, once the table is transposed so genes become columns), the maximum-likelihood
    lambda explodes and the transform overflows float64. The result is either non-finite -- which used to
    reach sklearn's PCA as ``ValueError: Input X contains NaN``, a message that names nothing the user can act
    on -- or a finite but astronomical value that silently dominates every principal component.
    """
    with warnings.catch_warnings():  # comparing NaNs is exactly what we are here to detect
        warnings.simplefilter('ignore')
        unstable = ~np.isfinite(box_cox_array).all(axis=0) | (
            np.abs(np.nan_to_num(box_cox_array)).max(axis=0) > BOX_COX_MAX_ABS_VALUE
        )
    if not unstable.any():
        return
    unstable_indices = np.nonzero(unstable)[0]
    if feature_names is None:
        names = [f'#{ind + 1}' for ind in unstable_indices]
    else:
        names = [str(feature_names[ind]) for ind in unstable_indices]
    n_unstable = len(names)
    listed = ', '.join(f"'{name}'" for name in names[:_MAX_NAMES_IN_ERROR])
    if n_unstable > _MAX_NAMES_IN_ERROR:
        listed += f', and {n_unstable - _MAX_NAMES_IN_ERROR} more'
    noun, these = ('gene', 'this gene') if n_unstable == 1 else ('genes', 'these genes')
    raise InvalidValueError(
        f'Box-Cox transformation is numerically unstable for {n_unstable} {noun} in this '
        f'table ({listed}): their values are near-constant at high magnitude. '
        f"Filter {these} out, or choose a different transform (e.g. 'log')."
    )


def standard_box_cox(
    data: Union[np.ndarray, pl.DataFrame],
    parallel_backend: str = 'sequential',
    feature_names: Optional[typing.Sequence[str]] = None,
) -> Union[np.ndarray, pl.DataFrame]:
    """
    Apply a per-column Box-Cox power transform followed by standardization.

    :param data: the data to transform (columns are transformed independently).
    :type data: np.ndarray or pl.DataFrame
    :param parallel_backend: joblib backend used to parallelize the per-column Box-Cox across columns. \
    The default ('sequential') keeps the original single-process behavior; any other backend parallelizes \
    once the number of columns reaches ``BOX_COX_PARALLEL_MIN_COLUMNS``. The result is identical either way \
    up to floating-point precision.
    :type parallel_backend: str (default='sequential')
    :param feature_names: names of the transformed columns, used only to name the offending genes if the \
    transform turns out to be numerically unstable. Defaults to the numeric column names of a DataFrame, \
    or to 1-based column positions for a bare array.
    :type feature_names: sequence of str or None (default=None)
    :raises InvalidValueError: if the Box-Cox transform is numerically unstable for one or more columns.
    """
    if isinstance(data, pl.DataFrame):
        numeric_columns = data.select(cs.numeric()).columns
        array = data.select(cs.numeric()).to_numpy()
        if feature_names is None:
            feature_names = numeric_columns
    else:
        array = data
    if parallel_backend not in (None, 'sequential') and array.shape[1] >= BOX_COX_PARALLEL_MIN_COLUMNS:
        box_cox_array = _parallel_box_cox(array, backend=parallel_backend, n_jobs=joblib.cpu_count())
    else:
        box_cox_array = _box_cox_fit_transform(array)
    _validate_box_cox_stability(box_cox_array, feature_names)
    res_array = _sklearn.preprocessing.StandardScaler().fit_transform(box_cox_array)
    if isinstance(data, pl.DataFrame):
        return data.select(~cs.numeric()).with_columns(pl.DataFrame(res_array, schema=numeric_columns))
    return res_array


def standard_log(data: Union[np.ndarray, pl.DataFrame]) -> Union[np.ndarray, pl.DataFrame]:
    """
    Apply a log2(x+1) transform followed by standardization.

    This is the robust alternative to the Box-Cox power transform: it uses one fixed formula for the whole
    table instead of fitting a separate exponent per gene, so it cannot blow up on a near-constant,
    high-magnitude gene the way Box-Cox can.

    :param data: the data to transform (must not contain values below -1).
    :type data: np.ndarray or pl.DataFrame
    :raises InvalidValueError: if the data contains values below -1, for which log2(x+1) is undefined.
    """
    if isinstance(data, pl.DataFrame):
        array = data.select(cs.numeric()).to_numpy()
    else:
        array = data
    if np.nanmin(array) <= -1:
        raise InvalidValueError(
            "The 'log' transform is only defined for values greater than -1 "
            '(it computes log2(x+1)), but this table contains smaller values. '
            'Choose a different transform, or filter out the offending values.'
        )
    log_array = np.log2(array + 1)
    res_array = _sklearn.preprocessing.StandardScaler().fit_transform(log_array)
    if isinstance(data, pl.DataFrame):
        return data.select(~cs.numeric()).with_columns(
            pl.DataFrame(res_array, schema=data.select(cs.numeric()).columns)
        )
    return res_array


def parse_power_transform(power_transform: Union[bool, str]) -> str:
    """
    Normalize a ``power_transform`` argument into one of the names in ``param_typing.POWER_TRANSFORMS``.

    ``power_transform`` used to be a boolean. ``True``/``False`` (in their Python and numpy spellings alike)
    are therefore still accepted, and always will be, so that Pipelines and exported parameter files saved by
    older versions of *RNAlysis* keep running with exactly the behavior they had when they were saved.

    :param power_transform: the transform to apply: 'box-cox', 'log', 'none', or the legacy True/False.
    :type power_transform: 'box-cox', 'log', 'none', or bool
    :return: the normalized transform name.
    :rtype: str
    """
    # numpy's boolean scalar is accepted alongside Python's, because that is what a comparison over a numpy or
    # polars column yields, and passing one straight in worked before 4.3 (np.bool_ stopped being a bool in numpy
    # 2, and the pre-4.3 code tested truthiness rather than the type). Other truthy/falsey values are not accepted
    # -- 1, 0 and 1.0 are not spellings of a transform, and are rejected below.
    if isinstance(power_transform, (bool, np.bool_)):
        return 'box-cox' if power_transform else 'none'
    if isinstance(power_transform, str) and power_transform.lower() in POWER_TRANSFORMS:
        return power_transform.lower()
    raise InvalidValueError(
        f"Invalid value for 'power_transform': {power_transform!r}. "
        f"'power_transform' must be one of {list(POWER_TRANSFORMS)}."
    )


def get_transform_function(power_transform: Union[bool, str]) -> Callable:
    """
    Return the transform+standardization function matching a ``power_transform`` argument.

    :param power_transform: the transform to apply: 'box-cox', 'log', 'none', or the legacy True/False.
    :type power_transform: 'box-cox', 'log', 'none', or bool
    """
    method = parse_power_transform(power_transform)
    if method == 'box-cox':
        return standard_box_cox
    elif method == 'log':
        return standard_log
    return standardize


def transform_and_standardize(
    data: Union[np.ndarray, pl.DataFrame],
    power_transform: Union[bool, str] = 'box-cox',
    parallel_backend: str = 'sequential',
    feature_names: Optional[typing.Sequence[str]] = None,
) -> Union[np.ndarray, pl.DataFrame]:
    """
    Standardize the data, optionally applying a power/log transform to it first.

    :param data: the data to transform (columns are transformed independently).
    :type data: np.ndarray or pl.DataFrame
    :param power_transform: the transform to apply: 'box-cox', 'log', 'none', or the legacy True/False.
    :type power_transform: 'box-cox', 'log', 'none', or bool (default='box-cox')
    :param parallel_backend: joblib backend used to parallelize the per-column Box-Cox across columns. \
    Ignored by the other transforms.
    :type parallel_backend: str (default='sequential')
    :param feature_names: names of the transformed columns, used only to name the offending genes if the \
    Box-Cox transform turns out to be numerically unstable.
    :type feature_names: sequence of str or None (default=None)
    """
    if parse_power_transform(power_transform) == 'box-cox':
        return standard_box_cox(data, parallel_backend=parallel_backend, feature_names=feature_names)
    return get_transform_function(power_transform)(data)


def box_cox_parallel_backend() -> str:
    """Return a joblib backend suitable for parallelizing the Box-Cox transform in the current environment.

    Frozen (PyInstaller) builds cannot use the ``loky`` backend, so they fall back to ``multiprocessing``
    (see ``rnalysis.utils.param_typing.PARALLEL_BACKENDS``). This is only ever used from top-level, non-nested
    call sites (the PCA methods), so ``multiprocessing``'s no-nested-children limitation does not apply.
    """
    return 'multiprocessing' if FROZEN_ENV else 'loky'


def shift_to_baseline(data: Union[np.ndarray, pl.DataFrame], baseline: float = 0) -> Union[np.ndarray, pl.DataFrame]:
    """

    :param data:
    :type data:
    :param baseline:
    :type baseline:
    :return:
    :rtype:
    """
    if isinstance(data, pl.DataFrame):
        array = data.select(cs.numeric()).to_numpy()
    elif isinstance(data, pl.Series):
        array = data.to_numpy()
    else:
        array = data
    min_val = array.min()

    while isinstance(min_val, np.ndarray):
        min_val = min_val.min()
    diff = baseline - min_val
    return data + diff


def standardize(data: Union[np.ndarray, pl.DataFrame]) -> Union[np.ndarray, pl.DataFrame]:
    """

    :param data:
    :type data:
    :return:
    :rtype:
    """
    if isinstance(data, pl.DataFrame):
        numeric_columns = data.select(cs.numeric()).columns
        array = data.select(cs.numeric()).to_numpy()
    else:
        array = data
    res_array = _sklearn.preprocessing.StandardScaler().fit_transform(array)
    if isinstance(data, pl.DataFrame):
        # only the numeric columns were transformed, so only their names may label the result -- passing every
        # column name (the gene-ID one included) raised a shape mismatch on any real RNAlysis table
        return data.select(~cs.numeric()).with_columns(pl.DataFrame(res_array, schema=numeric_columns))
    return res_array


def intersection_nonempty(*objs: Union[list, set, tuple]):
    return set.intersection(*[set(item) for item in objs if len(item) > 0])


def color_generator():
    """
    A generator that yields distinct colors up to a certain limit, and then yields randomized RGB values.

    :return: a color name string (like 'black', \
    or a numpy.ndarray of size (3,) containing three random values each between 0 and 1.

    """
    preset_colors = [
        'tab:blue',
        'tab:red',
        'tab:green',
        'tab:orange',
        'tab:purple',
        'tab:brown',
        'tab:pink',
        'tab:gray',
        'tab:olive',
        'tab:cyan',
        'gold',
        'maroon',
        'mediumslateblue',
        'fuchsia',
        'lawngreen',
        'moccasin',
        'thistle',
    ]
    for color in preset_colors:
        yield color
    while True:
        yield np.random.random(3)


@lru_cache(maxsize=2**16)
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
    multiplier = (1 / n_colors) + (1.6 / (2**n_colors)) / n_colors

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
        loglikelihood += n_points_cluster * (
            np.log(n_points_cluster)
            - np.log(n_points)
            - (n_dimensions / 2) * np.log(2 * math.pi * variance)
            - (n_points_cluster - 1) / 2
        )

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
    :return: sanitized name
    :rtype: str
    """
    new_name = name.rstrip().replace(' ', '_')
    if not new_name:
        return new_name
    if new_name[0].isdigit():
        new_name = 'var_' + new_name

    if not new_name.isalnum():
        new_name = ''.join([char if char.isalnum() else '_' for char in new_name])

    return new_name


class InteractiveScatterFigure(figure.Figure):
    def __init__(
        self, labels: typing.List[str], annotation_fontsize: float = 10, show_cursor: bool = True, *args, **kwargs
    ):
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
                self.is_labeled[this_ind] = plt.annotate(
                    self.labels[this_ind],
                    (np.take(xdata, this_ind), np.take(ydata, this_ind)),
                    xytext=(3, 3),
                    textcoords='offset points',
                    fontsize=self.annotation_fontsize,
                )

            self.canvas.draw()


#: types that ``yaml.safe_dump`` can always represent, and that ``yaml.safe_load`` returns as-is.
#: membership must be tested by exact type - PyYAML's SafeRepresenter dispatches on ``type(value)``,
#: so a *subclass* of any of these (an Enum mixin, a pandas Timestamp) is NOT representable.
_YAML_SAFE_SCALARS = (str, bool, int, float, bytes, type(None), datetime, date)


def _coerce_yaml_scalar_subclass(value):
    """
    Convert a subclass of a YAML-safe scalar down to the plain built-in PyYAML can represent.

    :param value: a value that is an instance - but not an exact instance - of a YAML-safe scalar
    :return: the equivalent plain built-in value
    """
    if isinstance(value, datetime):  # datetime subclasses date, so it must be checked first
        return datetime(
            value.year, value.month, value.day, value.hour, value.minute, value.second, value.microsecond, value.tzinfo
        )
    if isinstance(value, date):
        return date(value.year, value.month, value.day)
    for base in (bool, int, float, str, bytes):  # bool subclasses int, so it must be checked first
        if isinstance(value, base):
            return base(value)
    raise InternalError(f"'{type(value).__name__}' is not a subclass of a YAML-safe scalar type.")


def _sanitize_for_yaml(value, context: str):
    """
    Recursively convert a Pipeline parameter into an object ``yaml.safe_dump`` can represent.

    numpy scalars become their Python equivalents, numpy arrays and tuples become lists, and \
    pathlib Paths become POSIX-style strings. Values that cannot be represented at all raise a \
    typed error naming the offending function and parameter, instead of a bare ``RepresenterError``.

    :param value: the parameter value to sanitize
    :param context: human-readable description of the parameter, used in error messages
    :type context: str
    :return: a YAML-representable version of ``value``
    """
    if isinstance(value, np.generic):
        value = value.item()
    elif isinstance(value, np.ndarray):
        value = value.tolist()
    elif isinstance(value, Path):
        return value.as_posix()
    elif isinstance(value, enum.Enum):
        # an Enum member stands for its value; str()ing it would yield 'Mode.UNION', not 'union'
        return _sanitize_for_yaml(value.value, context)

    if isinstance(value, dict):
        return {
            _sanitize_hashable(key, context, 'dictionary key'): _sanitize_for_yaml(val, context)
            for key, val in value.items()
        }
    if isinstance(value, (list, tuple)):
        # YAML has no tuple type - tuples are stored (and re-loaded) as lists
        return [_sanitize_for_yaml(item, context) for item in value]
    if isinstance(value, (set, frozenset)):
        return {_sanitize_hashable(item, context, 'set item') for item in value}
    if type(value) in _YAML_SAFE_SCALARS:
        return value
    if isinstance(value, _YAML_SAFE_SCALARS):
        return _coerce_yaml_scalar_subclass(value)

    try:
        yaml.safe_dump(value)
    except yaml.YAMLError as err:
        raise InvalidTypeError(
            f'Cannot export Pipeline: the value of {context} (of type '
            f"'{type(value).__name__}') cannot be saved to a Pipeline YAML file. "
            f'Please replace it with a simple value, such as a number, a string, or a list. '
        ) from err
    return value


def _sanitize_hashable(value, context: str, kind: str):
    """
    Sanitize a value that must stay hashable (a dictionary key or a set item).

    :param value: the value to sanitize
    :param context: human-readable description of the parameter it belongs to, used in error messages
    :type context: str
    :param kind: what the value is ('dictionary key' or 'set item'), used in error messages
    :type kind: str
    :return: a YAML-representable, hashable version of ``value``
    """
    sanitized = _sanitize_for_yaml(value, context)
    try:
        hash(sanitized)
    except TypeError as err:
        raise InvalidTypeError(
            f'Cannot export Pipeline: the {kind} {value!r} of {context} '
            f'cannot be saved to a Pipeline YAML file. '
            f'Please use a simple value such as a string or a number. '
        ) from err
    return sanitized


def _parse_pipeline_yaml_string(content: str) -> Tuple[bool, typing.Any, Optional[yaml.YAMLError]]:
    """
    Determine whether the given string plausibly *is* a Pipeline YAML document, and if so, parse it.

    This is what keeps a mistyped filename ('C:/no/such/file.yaml') from being silently parsed as \
    YAML and failing later with an unrelated Python error.

    :param content: the string to examine
    :type content: str
    :return: a tuple of (whether the string was meant as a YAML document, the parsed document, \
    and the error raised while parsing it, if any)
    """
    # a multi-line string was clearly meant as file *contents* rather than as a file name
    is_document = '\n' in content
    try:
        parsed = yaml.safe_load(content)
    except yaml.YAMLError as err:
        return is_document, None, err
    # a single-line string is only treated as YAML if it parses into a Pipeline-shaped mapping
    return is_document or isinstance(parsed, dict), parsed, None


class GenericPipeline(abc.ABC):
    __slots__ = {'functions': 'list of functions to perform', 'params': 'list of function parameters'}

    def __init__(self):
        self.functions = []
        self.params = []

    def __len__(self):
        return len(self.functions)

    def __eq__(self, other):
        if self.functions == other.functions and self.params == other.params:
            return True
        return False

    @staticmethod
    def _param_string(args: tuple, kwargs: dict):
        """
        Returns a formatted string of the given arguments and keyworded arguments.

        :param args: arguments to format as string
        :type args: tuple
        :param kwargs: keyworded arguments to format as string
        :type kwargs: dict
        :return: a formatted string of arguments and keyworded argumentss
        :rtype: str

        """
        args_str = ', '.join([f"'{arg}'" if isinstance(arg, str) else f'{arg}' for arg in args])
        kwargs_str = ', '.join(
            [f"{key}='{arg}'" if isinstance(arg, str) else f'{key}={arg}' for key, arg in kwargs.items()]
        )
        if len(args_str) == 0:
            return kwargs_str
        elif len(kwargs_str) == 0:
            return args_str
        else:
            return f'{args_str}, {kwargs_str}'

    def remove_last_function(self):
        """
        Removes from the Pipeline the last function that was added to it. Removal is in-place.

        :Examples:
            >>> from rnalysis import filtering
            >>> pipe = filtering.Pipeline()
            >>> pipe.add_function(filtering.Filter.filter_missing_values)
            Added function 'Filter.filter_missing_values()' to the pipeline.
            >>> pipe.remove_last_function()
            Removed function filter_missing_values with parameters [] from the pipeline.

        """
        if not (len(self.functions) > 0 and len(self.params) > 0):
            raise InvalidValueError('Pipeline is empty, no functions to remove!')
        func = self.functions.pop(-1)
        args, kwargs = self.params.pop(-1)
        print(
            f'Removed function {func.__name__} with parameters [{self._param_string(args, kwargs)}] from the pipeline.'
        )

    def _readable_func_signature(self, func: types.FunctionType, args: tuple, kwargs: dict):
        """
        Returns a human-readable string functions signature for the given function and arguments.

        :param func: the function or method to generate signature for
        :type func: function
        :param args: arguments given for the function
        :type args: tuple
        :param kwargs: keyworded arguments given for the function
        :type kwargs: dict
        :return: function signature string
        :rtype: str
        """
        return f'{get_method_readable_name(func)}: ({self._param_string(args, kwargs)})'

    def export_pipeline(self, filename: Union[str, Path, None]) -> Union[None, str]:
        """
        Export a Pipeline to a Pipeline YAML file or YAML-like string.

        Function parameters are converted into their closest YAML equivalent before saving: \
        numpy scalars and arrays become plain Python numbers and lists, and pathlib Paths become \
        strings. Note that YAML has no tuple type - tuples are therefore saved (and re-loaded) as \
        lists, except for the top-level tuple of positional arguments, which is restored on import.

        :param filename: filename to save the Pipeline YAML to, or None to return a YAML-like string instead.
        :type filename: str, pathlib.Path, or None
        :return: if filename is None, returns the Pipeline YAML-like string.
        """
        pipeline_dict = self._get_pipeline_dict()
        for func, (args, kwargs) in zip(self.functions, self.params):
            pipeline_dict['functions'].append(func.__name__)
            pipeline_dict['params'].append(self._sanitize_params(func.__name__, args, kwargs))
        # the whole Pipeline is serialized *before* the target file is opened: opening it for
        # writing truncates it, so a Pipeline that cannot be serialized would otherwise destroy the
        # Pipeline file the user was overwriting.
        pipeline_yaml = yaml.safe_dump(pipeline_dict)
        if filename is None:
            return pipeline_yaml
        with open(filename, 'w') as f:
            f.write(pipeline_yaml)

    @staticmethod
    def _sanitize_params(func_name: str, args: tuple, kwargs: dict) -> list:
        """
        Return a YAML-representable version of a single function's arguments.

        :param func_name: name of the function the arguments belong to (used in error messages)
        :type func_name: str
        :param args: unkeyworded arguments of the function
        :type args: tuple
        :param kwargs: keyworded arguments of the function
        :type kwargs: dict
        :return: a [args, kwargs] list that ``yaml.safe_dump`` can represent
        """
        sanitized_args = [
            _sanitize_for_yaml(arg, f"positional argument #{i + 1} of function '{func_name}'")
            for i, arg in enumerate(args)
        ]
        sanitized_kwargs = {
            key: _sanitize_for_yaml(val, f"parameter '{key}' of function '{func_name}'") for key, val in kwargs.items()
        }
        return [sanitized_args, sanitized_kwargs]

    def _get_pipeline_dict(self):
        d = dict(
            functions=[],
            params=[],
            metadata={
                'rnalysis_version': f'{__version__}',
                'export_time': datetime.today().strftime('%Y/%m/%d, %H:%M:%S'),
            },
        )
        return d

    @classmethod
    def import_pipeline(cls, filename: Union[str, Path]) -> 'GenericPipeline':
        """
        Import a Pipeline from a Pipeline YAML file or YAML-like string.

        Note that YAML has no tuple type: tuples within the Pipeline's function parameters were \
        saved as lists, and are therefore returned as lists. The only exception is the top-level \
        tuple of positional arguments of each function, which is restored on import.

        :param filename: name of the YAML file containing the Pipeline, or a YAML-like string.
        :type filename: str or pathlib.Path
        :return: the imported Pipeline
        :rtype: Pipeline
        """
        pipeline_dict = cls._read_pipeline_dict(filename)
        pipeline = cls.__new__(cls)
        pipeline._init_from_dict(pipeline_dict)
        return pipeline

    @classmethod
    def _read_pipeline_dict(cls, filename: Union[str, Path]):
        """
        Read a Pipeline YAML document from either a file or a YAML-like string.

        :param filename: name of the YAML file containing the Pipeline, or a YAML-like string.
        :type filename: str or pathlib.Path
        :return: the parsed Pipeline YAML document
        """
        try:
            with open(filename) as f:
                try:
                    return yaml.safe_load(f)
                except yaml.YAMLError as err:
                    raise cls._pipeline_load_error(
                        None, f"the file '{filename}' is not a valid YAML document ({err})."
                    ) from err
        except OSError as err:
            if isinstance(filename, str):
                is_yaml_string, pipeline_dict, parse_error = _parse_pipeline_yaml_string(filename)
                if is_yaml_string:
                    if parse_error is not None:
                        raise cls._pipeline_load_error(
                            None, f'the given Pipeline is not a valid YAML document ({parse_error}).'
                        ) from parse_error
                    return pipeline_dict
            if isinstance(err, FileNotFoundError):
                raise FileNotFoundError(
                    f"Could not find the Pipeline file '{filename}'. "
                    f'Please make sure the file exists, '
                    f'and that its path is spelled correctly. '
                ) from err
            raise

    @staticmethod
    def _exported_version_string(pipeline_dict) -> str:
        """
        Describe the RNAlysis version that exported the given Pipeline, based on its version stamp.

        :param pipeline_dict: a parsed Pipeline YAML document (not necessarily a valid one)
        :return: a human-readable description of the exporting RNAlysis version
        :rtype: str
        """
        version = None
        if isinstance(pipeline_dict, dict):
            metadata = pipeline_dict.get('metadata')
            if isinstance(metadata, dict):
                version = metadata.get('rnalysis_version')
        if version is None:
            return 'this Pipeline does not record which version of RNAlysis exported it'
        return f'this Pipeline was exported by RNAlysis {version}'

    @classmethod
    def _pipeline_load_error(cls, pipeline_dict, reason: str) -> InvalidValueError:
        """
        Build a Pipeline load-failure error that names the exporting and the current RNAlysis version.

        :param pipeline_dict: a parsed Pipeline YAML document (not necessarily a valid one)
        :param reason: description of what went wrong
        :type reason: str
        :return: the error to raise
        :rtype: InvalidValueError
        """
        return InvalidValueError(
            f'Failed to load Pipeline: {reason} '
            f'({cls._exported_version_string(pipeline_dict)}, '
            f'and the current version is RNAlysis {__version__})'
        )

    @classmethod
    def _validate_pipeline_dict(cls, pipeline_dict):
        """
        Verify that the given parsed YAML document has the shape of a Pipeline.

        :param pipeline_dict: a parsed Pipeline YAML document
        """
        if not isinstance(pipeline_dict, dict):
            raise cls._pipeline_load_error(
                pipeline_dict,
                f"expected a Pipeline YAML file or YAML-like string, but got '{type(pipeline_dict).__name__}' instead.",
            )
        for key in ('functions', 'params'):
            if key not in pipeline_dict:
                raise cls._pipeline_load_error(pipeline_dict, f"the Pipeline is missing the mandatory '{key}' field.")
            if not isinstance(pipeline_dict[key], list):
                raise cls._pipeline_load_error(
                    pipeline_dict,
                    f"the Pipeline's '{key}' field must be a list, "
                    f"but is '{type(pipeline_dict[key]).__name__}' instead.",
                )
        if len(pipeline_dict['functions']) != len(pipeline_dict['params']):
            raise cls._pipeline_load_error(
                pipeline_dict,
                f'the Pipeline lists {len(pipeline_dict["functions"])} functions, '
                f'but {len(pipeline_dict["params"])} sets of parameters.',
            )

    @classmethod
    def _params_from_dict(cls, pipeline_dict: dict) -> list:
        """
        Parse the 'params' field of a Pipeline YAML document into a list of (args, kwargs) pairs.

        :param pipeline_dict: a parsed Pipeline YAML document
        :return: the Pipeline's function parameters
        :rtype: list
        """
        params = []
        for entry in pipeline_dict['params']:
            if not isinstance(entry, (list, tuple)) or len(entry) != 2:
                raise cls._pipeline_load_error(
                    pipeline_dict,
                    f'the function parameters {entry} are malformed - '
                    f'expected a pair of [arguments, keyword arguments].',
                )
            args, kwargs = entry
            if kwargs is None:
                kwargs = {}
            if not isinstance(kwargs, dict):
                raise cls._pipeline_load_error(
                    pipeline_dict, f'the keyword arguments {kwargs} are malformed - expected a dictionary.'
                )
            args = tuple(args) if isinstance(args, (list, tuple, set)) else (args,)
            params.append((args, kwargs))
        return params

    @classmethod
    def _resolve_function(cls, func_name, namespace, namespace_name: str, pipeline_dict: dict) -> types.FunctionType:
        """
        Look up a Pipeline function by the name recorded in a Pipeline YAML document.

        :param func_name: the function name recorded in the Pipeline YAML document
        :param namespace: the class or module the function should belong to
        :param namespace_name: human-readable name of the namespace, used in error messages
        :type namespace_name: str
        :param pipeline_dict: the parsed Pipeline YAML document, used to report the exporting version
        :return: the function matching the given name
        """
        # this check is deliberately no stricter than the one in add_function (a function object of
        # the right namespace): anything that could be exported must import back (hard rule 4).
        if not isinstance(func_name, str):
            raise cls._pipeline_load_error(pipeline_dict, f"'{func_name}' is not a valid RNAlysis function name.")
        func = getattr(namespace, func_name, None)
        if not isinstance(func, types.FunctionType):
            raise cls._pipeline_load_error(
                pipeline_dict,
                f"function '{func_name}' does not exist in {namespace_name}. "
                f'It may have been renamed or removed since this Pipeline was exported.',
            )
        return func

    def _init_from_dict(self, pipeline_dict: dict):
        raise NotImplementedError

    def _func_signature(self, func: types.FunctionType, args: tuple, kwargs: dict):
        """
        Returns a string functions signature for the given function and arguments.

        :param func: the function or method to generate signature for
        :type func: function
        :param args: arguments given for the function
        :type args: tuple
        :param kwargs: keyworded arguments given for the function
        :type kwargs: dict
        :return: function signature string
        :rtype: str
        """
        return f'{func.__name__}({self._param_string(args, kwargs)})'

    def add_function(self, func: types.FunctionType, *args, **kwargs):
        if not isinstance(func, types.FunctionType):
            raise InvalidTypeError(f"'func' must be a function, is {type(func)} instead.")

        self.functions.append(func)
        self.params.append((args, kwargs))
        print(f"Added function '{self._func_signature(func, args, kwargs)}' to the pipeline.")

    def _validate_pipeline(self):
        if not (len(self.functions) > 0 and len(self.params) > 0):
            raise InvalidValueError('Cannot apply an empty pipeline!')
        if len(self.functions) != len(self.params):
            raise InvalidValueError("Cannot apply Pipeline: length of 'functions' different from length of 'params'!")

    @abc.abstractmethod
    def apply_to(self, *args, **kwargs):
        raise NotImplementedError


def jitter(n: int, jitter_range: float = 0.2) -> np.ndarray:
    if n == 0:
        return np.array([])
    if jitter_range == 0:
        return np.zeros(n)

    jittered = np.random.uniform(-jitter_range, jitter_range, n)
    jittered -= np.mean(jittered)  # center the jittered points
    scaling_factor = jitter_range / max(np.max(np.abs(jittered)), 1e-6)
    jittered *= scaling_factor  # scale the jitter points to the jitter range
    return jittered
