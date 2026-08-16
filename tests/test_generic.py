from types import ModuleType

import pytest
from matplotlib.backend_bases import PickEvent

# scikit-learn is imported lazily by rnalysis.utils.generic (see tests/test_imports.py), so the
# star-import below no longer re-exports it -- the tests import it directly instead.
from sklearn.preprocessing import PowerTransformer, StandardScaler

from rnalysis import FROZEN_ENV
from rnalysis.exceptions import InvalidValueError
from rnalysis.utils import generic
from rnalysis.utils.generic import *


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


def test_standardize_dataframe_with_a_gene_id_column():
    # get_transform_function hands 'box-cox', 'log' and 'none' out as interchangeable, and every RNAlysis table
    # carries a non-numeric gene-ID column -- so 'none' has to survive the same input the other two already handle.
    # standardize() used to label its result frame with *all* the column names while holding only the numeric ones.
    array = np.array([[10.0, 3.0], [50.0, 90.0], [200.0, 40.0], [1000.0, 7.0]])
    data_df = pl.DataFrame([f'gene{i}' for i in range(4)]).with_columns(pl.DataFrame(array, schema=['cond1', 'cond2']))

    res_df = standardize(data_df)

    assert isinstance(res_df, pl.DataFrame)
    assert res_df.columns == data_df.columns
    assert np.all(res_df.select(pl.first()) == data_df.select(pl.first()))
    assert np.array_equal(res_df.select(cs.numeric()).to_numpy(), StandardScaler().fit_transform(array))


def test_standard_box_cox():
    np.random.seed(42)
    data = np.random.randint(-200, 100000, (100, 5))
    res = standard_box_cox(data)
    assert res.shape == data.shape
    assert np.isclose(res.mean(axis=0), 0).all()
    assert np.isclose(res.std(axis=0), 1).all()
    assert not np.isclose(res, standardize(data)).all()

    data_df = pl.DataFrame([f'ind{i}' for i in range(100)]).with_columns(
        pl.DataFrame(data, schema=[f'col{j}' for j in range(5)])
    )
    res_df = standard_box_cox(data_df)
    assert isinstance(res_df, pl.DataFrame)
    assert res_df.shape == data_df.shape
    assert np.all(res_df.select(pl.first()) == data_df.select(pl.first()))
    assert np.all(res_df.columns == data_df.columns)


# Box-Cox each column independently, so splitting the columns across workers gives the same result as the
# single-process path. It is not guaranteed bit-for-bit: transforming a column-slice of a C-contiguous matrix
# vs the full matrix can shift numpy's reduction blocking, perturbing some data by ~1e-15 (measured on the
# real count-matrix fixture; far below anything that affects PCA loadings, cluster assignments, or sort order).
# The end-to-end loky path is covered against a committed truth file by test_split_by_principal_components.
_BOX_COX_PARALLEL_TOL = 1e-10


@pytest.mark.parametrize('backend', ['threading'])
def test_parallel_box_cox_matches_serial(backend):
    rng = np.random.default_rng(0)
    array = rng.integers(1, 5000, (6, 40)).astype(float)
    serial = generic._box_cox_fit_transform(array)
    parallel = generic._parallel_box_cox(array, backend=backend, n_jobs=4)
    assert parallel.shape == serial.shape
    assert np.allclose(parallel, serial, rtol=0, atol=_BOX_COX_PARALLEL_TOL)


def test_standard_box_cox_parallel_equals_sequential(monkeypatch):
    # lower the parallelization threshold so a small (fast) array still exercises the parallel branch
    monkeypatch.setattr(generic, 'BOX_COX_PARALLEL_MIN_COLUMNS', 4)
    rng = np.random.default_rng(2)
    array = rng.integers(1, 5000, (6, 30)).astype(float)

    sequential = standard_box_cox(array, parallel_backend='sequential')
    parallel = standard_box_cox(array, parallel_backend='threading')
    assert np.allclose(parallel, sequential, rtol=0, atol=_BOX_COX_PARALLEL_TOL)

    data_df = pl.DataFrame([f'ind{i}' for i in range(6)]).with_columns(
        pl.DataFrame(array, schema=[f'col{j}' for j in range(30)])
    )
    seq_df = standard_box_cox(data_df, parallel_backend='sequential')
    par_df = standard_box_cox(data_df, parallel_backend='threading')
    assert np.all(par_df.select(pl.first()) == seq_df.select(pl.first()))
    assert np.allclose(
        par_df.select(cs.numeric()).to_numpy(),
        seq_df.select(cs.numeric()).to_numpy(),
        rtol=0,
        atol=_BOX_COX_PARALLEL_TOL,
    )


def test_standard_box_cox_default_is_sequential(monkeypatch):
    # the default call (no parallel_backend) must be byte-identical to the pre-existing serial path
    monkeypatch.setattr(generic, 'BOX_COX_PARALLEL_MIN_COLUMNS', 4)
    rng = np.random.default_rng(3)
    array = rng.integers(1, 5000, (6, 30)).astype(float)
    reference = StandardScaler().fit_transform(PowerTransformer(method='box-cox').fit_transform(array + 1))
    assert np.array_equal(standard_box_cox(array), reference)


# The real-world offending column (a GFP reporter at ~8,400 counts across 4 grouped samples), salvaged from the
# closed PR #238 (branch fix/boxcox-overflow-pca-nan) together with its non-finite detection logic. Box-Cox fits
# one lambda per column, and across only 4 near-constant high-magnitude values the MLE lambda explodes (~78), so
# the transform overflows float64. Instead of handing sklearn a matrix full of NaNs (which used to crash PCA with
# "Input X contains NaN"), RNAlysis now fails loudly and names the offending genes.
UNSTABLE_BOX_COX_COLUMN = np.array([8482.73, 8459.93, 8423.67, 8288.78])
WELL_BEHAVED_COLUMNS = np.array([[10.0, 3.0], [50.0, 90.0], [200.0, 40.0], [1000.0, 7.0]])


def _array_with_unstable_column():
    return np.column_stack([WELL_BEHAVED_COLUMNS[:, 0], UNSTABLE_BOX_COX_COLUMN, WELL_BEHAVED_COLUMNS[:, 1]])


@pytest.mark.parametrize('parallel_backend', ['sequential', 'threading'])
def test_standard_box_cox_raises_on_non_finite_column(parallel_backend, monkeypatch):
    monkeypatch.setattr(generic, 'BOX_COX_PARALLEL_MIN_COLUMNS', 2)
    array = _array_with_unstable_column()
    with pytest.raises(InvalidValueError) as err:
        standard_box_cox(array, parallel_backend=parallel_backend, feature_names=['good1', 'GFP', 'good2'])
    msg = str(err.value)
    assert 'GFP' in msg
    assert 'good1' not in msg and 'good2' not in msg
    assert "'log'" in msg


def test_standard_box_cox_unstable_column_error_names_dataframe_columns():
    array = _array_with_unstable_column()
    data_df = pl.DataFrame(['sample1', 'sample2', 'sample3', 'sample4']).with_columns(
        pl.DataFrame(array, schema=['good1', 'GFP', 'good2'])
    )
    with pytest.raises(InvalidValueError) as err:
        standard_box_cox(data_df)
    assert 'GFP' in str(err.value)


def test_standard_box_cox_raises_on_absurdly_large_output(monkeypatch):
    # a Box-Cox output can also come back *finite* but astronomically large (an exploded lambda that did not
    # quite overflow). Such a column silently dominates every principal component instead of crashing, so it
    # must be caught too. Which real inputs land there depends on the exact scipy/sklearn optimizer, so the
    # magnitude branch is pinned here by faking the transform output rather than by a brittle input fixture.
    array = _array_with_unstable_column()

    def fake_box_cox(arr):
        out = np.ones_like(arr)
        out[:, 1] = np.array([1e200, -1e200, 1e199, -1e199])
        return out

    monkeypatch.setattr(generic, '_box_cox_fit_transform', fake_box_cox)
    with pytest.raises(InvalidValueError) as err:
        standard_box_cox(array, feature_names=['good1', 'GFP', 'good2'])
    msg = str(err.value)
    assert 'GFP' in msg
    assert 'good1' not in msg and 'good2' not in msg  # only the offending gene is named


def test_standard_box_cox_well_behaved_data_is_unaffected():
    # rule 5: currently-working inputs must keep producing bit-identical output
    reference = StandardScaler().fit_transform(
        PowerTransformer(method='box-cox').fit_transform(WELL_BEHAVED_COLUMNS + 1)
    )
    assert np.array_equal(standard_box_cox(WELL_BEHAVED_COLUMNS, feature_names=['good1', 'good2']), reference)


@pytest.mark.parametrize(
    'value,truth',
    [
        (True, 'box-cox'),
        (False, 'none'),
        ('box-cox', 'box-cox'),
        ('Box-Cox', 'box-cox'),
        ('log', 'log'),
        ('none', 'none'),
    ],
)
def test_parse_power_transform(value, truth):
    assert parse_power_transform(value) == truth


@pytest.mark.parametrize('value,truth', [(np.True_, 'box-cox'), (np.False_, 'none')])
def test_parse_power_transform_accepts_numpy_booleans(value, truth):
    # `isinstance(np.True_, bool)` is False on numpy >= 2, but the pre-4.3 code was a plain truthiness test, so a
    # numpy boolean (which is what a comparison over a numpy/polars column yields) used to work. The promise that
    # True/False stay valid indefinitely has to keep covering the numpy flavour too.
    assert parse_power_transform(value) == truth


@pytest.mark.parametrize('value', ['boxcox', 'log2', '', None, 5])
def test_parse_power_transform_invalid(value):
    with pytest.raises(InvalidValueError):
        parse_power_transform(value)


def test_standard_log():
    rng = np.random.default_rng(4)
    array = rng.integers(0, 5000, (6, 8)).astype(float)
    reference = StandardScaler().fit_transform(np.log2(array + 1))
    assert np.array_equal(standard_log(array), reference)

    data_df = pl.DataFrame([f'ind{i}' for i in range(6)]).with_columns(
        pl.DataFrame(array, schema=[f'col{j}' for j in range(8)])
    )
    res_df = standard_log(data_df)
    assert isinstance(res_df, pl.DataFrame)
    assert res_df.columns == data_df.columns
    assert np.all(res_df.select(pl.first()) == data_df.select(pl.first()))
    assert np.array_equal(res_df.select(cs.numeric()).to_numpy(), reference)


def test_standard_log_survives_unstable_box_cox_data():
    # 'log' is the escape hatch offered by the Box-Cox error message -- it must never fail on data that
    # Box-Cox chokes on
    res = standard_log(_array_with_unstable_column())
    assert np.isfinite(res).all()


def test_standard_log_rejects_negative_values():
    with pytest.raises(InvalidValueError):
        standard_log(np.array([[1.0, -5.0], [2.0, 3.0]]))


@pytest.mark.parametrize(
    'power_transform,truth',
    [
        (True, 'standard_box_cox'),
        ('box-cox', 'standard_box_cox'),
        ('log', 'standard_log'),
        (False, 'standardize'),
        ('none', 'standardize'),
    ],
)
def test_get_transform_function(power_transform, truth):
    assert get_transform_function(power_transform).__name__ == truth


def test_color_generator():
    gen = color_generator()
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
    for i in range(150):
        color = next(gen)
        assert (isinstance(color, str) and color in preset_colors) or (
            isinstance(color, np.ndarray) and len(color) == 3 and np.max(color) <= 1 and np.min(color) >= 0
        )


@pytest.mark.parametrize(
    'this_set,other_sets,majority_threshold,truth',
    [
        ({1, 2, 3, 4}, [{1, 2, 3, 6}, {4, 5, 6}], 2 / 3, {1, 2, 3, 4, 6}),
        ({'a', 'ab', 'aab'}, [{'ba', 'b'}], 0.501, set()),
        ({'a', 'ab', 'aab'}, [{'ba', 'b'}], 0.5, {'a', 'ab', 'aab', 'ba', 'b'}),
        ({1, 2, 3}, [{2, 3, 4}, {3, 4, 5}], 0.5, {2, 3, 4}),
        ({1, 2, 3}, [{2, 3, 4}, {3, 4, 5}], 1, {3}),
    ],
)
def test_majority_vote_intersection(this_set, other_sets, majority_threshold, truth):
    result = SetWithMajorityVote.majority_vote_intersection(
        this_set, *other_sets, majority_threshold=majority_threshold
    )
    assert result == truth


@pytest.mark.parametrize('is_df', [True, False])
@pytest.mark.parametrize(
    'data,baseline,truth',
    [
        (np.array([1, 2, 3, 4, 5]), 0, np.array([0, 1, 2, 3, 4])),
        (
            np.array([[1, 2, 3], [-2, 4, 5], [0, 0, -1], [3, -2, 1]]),
            1,
            np.array([[4, 5, 6], [1, 7, 8], [3, 3, 2], [6, 1, 4]]),
        ),
        (np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]), -1, np.array([[[-1, 0], [1, 2]], [[3, 4], [5, 6]]])),
    ],
)
def test_shift_to_baseline(data, baseline, is_df, truth):
    if is_df and len(data.shape) <= 2:
        assert shift_to_baseline(pl.DataFrame(data), baseline).equals(pl.DataFrame(truth))
    else:
        assert np.all(shift_to_baseline(data, baseline) == truth)


def first_test_func():
    pass


def second_test_func(a: str, b: bool, c):
    pass


def third_test_func(a: typing.List[str], b: typing.Callable, c: int = 3, d=None):
    pass


class TestObj:
    def __init__(self):
        pass

    def fourth_test_func(self, a, b: None, c: float = 5.2):
        pass


@pytest.mark.parametrize(
    'func,obj,truth',
    [
        (first_test_func, None, {}),
        (
            second_test_func,
            None,
            {
                'a': {'annotation': str, 'default': inspect._empty},
                'b': {'annotation': bool, 'default': inspect._empty},
                'c': {'annotation': inspect._empty, 'default': inspect._empty},
            },
        ),
        (
            third_test_func,
            None,
            {
                'a': {'annotation': typing.List[str], 'default': inspect._empty},
                'b': {'annotation': typing.Callable, 'default': inspect._empty},
                'c': {'annotation': int, 'default': 3},
                'd': {'annotation': inspect._empty, 'default': None},
            },
        ),
        (
            'fourth_test_func',
            TestObj(),
            {
                'a': {'annotation': inspect._empty, 'default': inspect._empty},
                'b': {'annotation': None, 'default': inspect._empty},
                'c': {'annotation': float, 'default': 5.2},
            },
        ),
    ],
)
def test_get_signature(func, obj, truth):
    this_signature = get_method_signature(func, obj)
    assert len(this_signature) == len(truth)
    for key, val in truth.items():
        assert key in this_signature
        param = this_signature[key]
        assert param.name == key
        assert param.annotation == val['annotation']
        assert param.default == val['default']


@pytest.mark.parametrize(
    'intervals,expected',
    [
        ([], 0),
        ([(1, 3)], 3),
        ([(1, 3), (4, 6), (7, 10)], 10),
        ([(4, 6), (1, 3), (7, 10)], 10),
        ([(1, 4), (3, 6), (3, 5), (3, 6), (4, 9)], 9),
        ([(7, 10), (2, 5), (1, 4), (2, 5)], 9),
    ],
)
def test_sum_intervals_inclusive(intervals, expected):
    res = sum_intervals_inclusive(intervals)
    assert res == expected


@pytest.mark.parametrize(
    'seconds,expected',
    [
        (13, '00:13'),
        (0, '00:00'),
        (59.34, '00:59'),
        (60.0, '01:00'),
        (60.2, '01:00'),
        (192.17, '03:12'),
        (13 * 60 + 1.5, '13:01'),
        (100 * 60, '100:00'),
    ],
)
def test_format_time(seconds, expected):
    res = format_time(seconds)
    assert res == expected


@pytest.mark.parametrize(
    'name,expected',
    [
        ('aname', 'aname'),
        ('name123n', 'name123n'),
        ('camelCase', 'camelCase'),
        ('snake_case', 'snake_case'),
        ('name with spaces ', 'name_with_spaces'),
        ('1name of var', 'var_1name_of_var'),
        ('12345', 'var_12345'),
        ('%asdf&sign', '_asdf_sign'),
        (' ^^more things123 ', '___more_things123'),
    ],
)
def test_sanitize_variable_name(name, expected):
    res = sanitize_variable_name(name)
    assert res == expected


def test_get_method_readable_name():
    def func(x):
        return x + 1

    func.readable_name = 'readable name'
    assert get_method_readable_name(func) == 'readable name'

    def func2(a, b):
        return a + b

    assert get_method_readable_name(func2) == 'func2'

    func2.readable_name = 'readable name'
    assert get_method_readable_name(func2) == 'readable name'


def test_param_readable_names_decorator_sets_attribute():
    mapping = {'trim_n': 'Trim ambiguous (N) bases', 'drop_columns': 'Columns to drop'}

    @param_readable_names(mapping)
    def func(trim_n=True, drop_columns=None):
        return trim_n, drop_columns

    assert func.readable_param_names == mapping
    # the decorator must not alter the function's behavior
    assert func(False, ['a']) == (False, ['a'])


@pytest.mark.parametrize(
    'param_name,expected',
    [
        ('fastq_folder', 'FASTQ folder'),
        ('three_prime_adapters', 'Three prime adapters'),
        ('minimum_read_length', 'Minimum read length'),
        ('power_transform', 'Power transform'),
        ('drop_columns', 'Drop columns'),
        ('trim_n', 'Trim n'),
        ('single', 'Single'),
        ('go_id', 'GO ID'),
        ('kegg_pathways', 'KEGG pathways'),
        ('gene_ids', 'Gene IDs'),
        ('rna_type', 'RNA type'),
        ('mrna_length', 'mRNA length'),
        ('use_pca', 'Use PCA'),
        ('fdr_level', 'FDR level'),
        ('deseq2_normalization', 'DESeq2 normalization'),
        ('sam_path', 'SAM path'),
        ('bam_index', 'BAM index'),
        ('fastq_to_sam', 'FASTQ to SAM'),
        ('padj_threshold', 'padj threshold'),
        ('gtf_file', 'GTF file'),
        ('n_components', 'Number of components'),
        ('', ''),
    ],
)
def test_get_param_readable_name_auto(param_name, expected):
    assert get_param_readable_name(param_name) == expected


def test_get_param_readable_name_override():
    @param_readable_names({'trim_n': 'Trim ambiguous (N) bases'})
    def func(trim_n=True, other_param=1):
        pass

    # overridden parameter uses the declared label
    assert get_param_readable_name('trim_n', func) == 'Trim ambiguous (N) bases'
    # non-overridden parameters fall back to auto-humanization
    assert get_param_readable_name('other_param', func) == 'Other param'


def test_get_param_readable_name_func_without_overrides():
    def func(trim_n=True):
        pass

    # a function without a readable_param_names attribute auto-humanizes
    assert get_param_readable_name('trim_n', func) == 'Trim n'
    assert get_param_readable_name('trim_n', None) == 'Trim n'


def test_get_param_readable_name_override_beats_special_name():
    @param_readable_names({'n_components': 'Num comps'})
    def func(n_components=2):
        pass

    assert get_param_readable_name('n_components', func) == 'Num comps'


@pytest.mark.parametrize(
    'X, labels, expected_bic',
    [
        (np.array([[1, 2], [1.5, 2.5], [3, 4], [4, 5]]), np.array([0, 0, 1, 1]), -13.510391348997054),
        (np.array([[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]]), np.array([0, 0, 0, 1, 1]), -23.462198946075148),
        (
            np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4], [5, 5, 5]]),
            np.array([0, 0, 0, 1, 1]),
            -33.747038606183764,
        ),
    ],
)
def test_bic_score(X, labels, expected_bic):
    assert np.isclose(bic_score(X, labels), expected_bic)


def test_bic_score_single_cluster():
    X = np.array([[1, 2], [1.5, 2.5], [2, 3]])
    labels = np.array([0, 0, 0])
    m = np.mean(X, axis=0)
    var = np.sum((X[labels == 0] - m) ** 2) / float(3 - 1)
    log_likelihood = 3 * (-np.log(2 * np.pi * var) - 1)
    n_params = 3
    expected_bic = log_likelihood - 0.5 * n_params * np.log(3)
    assert np.isclose(bic_score(X, labels), expected_bic)


def test_bic_score_empty_input():
    X = np.empty((0, 2))
    labels = np.array([])
    assert np.isnan(bic_score(X, labels))


class MockMouseEvent:
    def __init__(self):
        self.guiEvent = None


class TestInteractiveScatterFigure:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.labels = ['A', 'B', 'C']
        self.fig = InteractiveScatterFigure(self.labels, annotation_fontsize=12, show_cursor=True)
        self.fig.ax.scatter([1, 2, 3], [4, 5, 6], picker=True)

    def test_initialization(self):
        assert isinstance(self.fig.ax, plt.Axes)
        assert self.fig.labels == ['A', 'B', 'C']
        assert self.fig.annotation_fontsize == 12
        assert self.fig.is_labeled == {}

    def test_on_exit(self):
        event = PickEvent('axes_leave_event', self.fig.canvas, mouseevent=MockMouseEvent(), artist=self.fig.ax)
        self.fig.on_exit(event)

    def test_on_move(self):
        event = PickEvent('motion_notify_event', self.fig.canvas, mouseevent=MockMouseEvent(), artist=self.fig.ax)
        self.fig.on_move(event)
        # No specific assertions, but ensure the canvas is redrawn

    def test_on_pick_add_annotation(self):
        event = PickEvent(
            'pick_event', self.fig.canvas, mouseevent=MockMouseEvent(), artist=self.fig.ax.collections[0], ind=[0]
        )
        self.fig.on_pick(event)
        assert 0 in self.fig.is_labeled
        assert isinstance(self.fig.is_labeled[0], plt.Annotation)
        assert self.fig.is_labeled[0].get_text() == 'A'

    def test_on_pick_remove_annotation(self):
        event = PickEvent(
            'pick_event', self.fig.canvas, mouseevent=MockMouseEvent(), artist=self.fig.ax.collections[0], ind=[0]
        )
        self.fig.on_pick(event)
        self.fig.on_pick(event)
        assert 0 not in self.fig.is_labeled

    def test_on_pick_invalid_artist(self):
        event = PickEvent(
            'pick_event', self.fig.canvas, mouseevent=MockMouseEvent(), artist=plt.Rectangle((0, 0), 1, 1), ind=[0]
        )
        self.fig.on_pick(event)
        assert self.fig.is_labeled == {}


@pytest.mark.parametrize('n', [1, 2, 3, 4, 5, 10, 100])
@pytest.mark.parametrize('jitter_range', [0, 0.01, 0.1, 0.5, 1, 5])
def test_jitter(n, jitter_range):
    res = jitter(n, jitter_range)
    assert len(res) == n
    assert isinstance(res, np.ndarray)
    assert np.isclose(np.mean(res), 0, atol=0.001)
    assert max(res) <= jitter_range or np.isclose(max(res), jitter_range, atol=0.001)
    assert min(res) >= -jitter_range or np.isclose(min(res), -jitter_range, atol=0.001)


# --- Pipeline YAML parameter sanitization -------------------------------------------------------
# Pipeline parameters are dumped with yaml.safe_dump, which can only represent a handful of types.
# Anything a user can plausibly pass (a numpy scalar plucked out of a result table, a pathlib Path)
# must be converted first, and anything that genuinely cannot be represented must fail with a
# message naming the offending function and parameter - not with a bare RepresenterError.


@pytest.mark.parametrize(
    'value,expected,expected_type',
    [
        (np.float64(5.0), 5.0, float),
        (np.int64(7), 7, int),
        (np.bool_(True), True, bool),
        (np.str_('abc'), 'abc', str),
        (Path('some/dir/table.csv'), 'some/dir/table.csv', str),
        (np.array([1, 2, 3]), [1, 2, 3], list),
        (np.array([[1.5], [2.5]]), [[1.5], [2.5]], list),
        ((1, 2), [1, 2], list),
        ('plain string', 'plain string', str),
        (None, None, type(None)),
    ],
)
def test_sanitize_for_yaml_converts_values(value, expected, expected_type):
    res = generic._sanitize_for_yaml(value, 'context')
    assert res == expected
    assert type(res) is expected_type
    assert yaml.safe_load(yaml.safe_dump(res)) == expected


def test_sanitize_for_yaml_recurses_into_containers():
    value = {'a': [np.float64(1.5), (Path('x/y'), {'b': np.array([1, 2])})], np.str_('c'): {np.int64(3)}}
    res = generic._sanitize_for_yaml(value, 'context')
    assert res == {'a': [1.5, ['x/y', {'b': [1, 2]}]], 'c': {3}}
    assert yaml.safe_load(yaml.safe_dump(res)) == res


def test_sanitize_for_yaml_unrepresentable_value_raises_typed_error():
    class Unrepresentable:
        pass

    with pytest.raises(InvalidTypeError) as err:
        generic._sanitize_for_yaml(Unrepresentable(), "parameter 'threshold' of function 'foo'")
    assert "parameter 'threshold' of function 'foo'" in str(err.value)
    assert 'Unrepresentable' in str(err.value)


def test_sanitize_for_yaml_unhashable_dict_key_raises_typed_error():
    with pytest.raises(InvalidTypeError) as err:
        generic._sanitize_for_yaml({(1, 2): 'val'}, "parameter 'grouping' of function 'foo'")
    assert "parameter 'grouping' of function 'foo'" in str(err.value)


@pytest.mark.parametrize(
    'content,is_yaml',
    [
        ('tests/test_files/no_such_file.yaml', False),
        ('C:/no/such/file.yaml', False),
        ('a_bare_word', False),
        ('', False),
        ('filter_type: countfilter', True),
        ('functions:\n- describe\nparams:\n- - []\n  - {}\n', True),
    ],
)
def test_parse_pipeline_yaml_string_only_accepts_plausible_yaml(content, is_yaml):
    assert generic._parse_pipeline_yaml_string(content)[0] == is_yaml


# --- numba disk cache (issue #257) ---------------------------------------------------------------
# Without `cache=True` every fresh process -- including every Windows `spawn` multiprocessing worker
# -- recompiles the jitted kernels from scratch. These guards keep that from silently regressing.


def test_numba_cache_flag_follows_frozen_env():
    """Caching must be on from source, and off in the frozen app (see generic.NUMBA_CACHE's comment)."""
    assert generic.NUMBA_CACHE is not FROZEN_ENV


@pytest.mark.skipif(not isinstance(generic.numba, ModuleType), reason='numba is not installed')
def test_jitted_kernels_use_the_disk_cache():
    """A NullCache on either kernel means it gets recompiled from scratch in every process."""
    from rnalysis.filtering import FoldChangeFilter
    from rnalysis.utils.enrichment_runner import PermutationTest

    for kernel in (FoldChangeFilter._foldchange_randomization, PermutationTest._calc_permutation_pval):
        assert type(kernel._cache).__name__ != 'NullCache', f'{kernel.__name__} is not disk-cached'
