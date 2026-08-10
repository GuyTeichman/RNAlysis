import pytest
from matplotlib.backend_bases import PickEvent

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


def test_standard_box_cox():
    np.random.seed(42)
    data = np.random.randint(-200, 100000, (100, 5))
    res = standard_box_cox(data)
    assert res.shape == data.shape
    assert np.isclose(res.mean(axis=0), 0).all()
    assert np.isclose(res.std(axis=0), 1).all()
    assert not np.isclose(res, standardize(data)).all()

    data_df = pl.DataFrame([f'ind{i}' for i in range(100)]).with_columns(
        pl.DataFrame(data, schema=[f'col{j}' for j in range(5)]))
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
        pl.DataFrame(array, schema=[f'col{j}' for j in range(30)]))
    seq_df = standard_box_cox(data_df, parallel_backend='sequential')
    par_df = standard_box_cox(data_df, parallel_backend='threading')
    assert np.all(par_df.select(pl.first()) == seq_df.select(pl.first()))
    assert np.allclose(par_df.select(cs.numeric()).to_numpy(),
                       seq_df.select(cs.numeric()).to_numpy(), rtol=0, atol=_BOX_COX_PARALLEL_TOL)


def test_standard_box_cox_default_is_sequential(monkeypatch):
    # the default call (no parallel_backend) must be byte-identical to the pre-existing serial path
    monkeypatch.setattr(generic, 'BOX_COX_PARALLEL_MIN_COLUMNS', 4)
    rng = np.random.default_rng(3)
    array = rng.integers(1, 5000, (6, 30)).astype(float)
    reference = StandardScaler().fit_transform(PowerTransformer(method='box-cox').fit_transform(array + 1))
    assert np.array_equal(standard_box_cox(array), reference)


@pytest.mark.parametrize('parallel_backend', ['sequential', 'threading'])
def test_standard_box_cox_survives_boxcox_overflow(parallel_backend, monkeypatch):
    # A near-constant, high-magnitude column with very few samples drives the per-column Box-Cox MLE
    # lambda to an extreme value whose transform overflows float64 -> non-finite output. That used to
    # crash PCA/clustering downstream (sklearn: "Input X contains NaN"). Such a column must be handled
    # gracefully (standardized without the power transform), while well-behaved columns stay untouched.
    monkeypatch.setattr(generic, 'BOX_COX_PARALLEL_MIN_COLUMNS', 2)
    overflow_col = np.array([8482.73, 8459.93, 8423.67, 8288.78])  # the real GFP-like offender
    good = np.array([[10.0, 3.0], [50.0, 90.0], [200.0, 40.0], [1000.0, 7.0]])
    array = np.column_stack([good[:, 0], overflow_col, good[:, 1]])

    res = standard_box_cox(array, parallel_backend=parallel_backend)
    assert res.shape == array.shape
    assert np.isfinite(res).all(), 'the power transform must not leak NaN/inf that crash PCA'

    # the well-behaved columns must be identical to their own Box-Cox+standardize result -- the
    # overflow fallback only touches the offending column.
    good_ref = StandardScaler().fit_transform(PowerTransformer(method='box-cox').fit_transform(good + 1))
    assert np.allclose(res[:, [0, 2]], good_ref, rtol=0, atol=1e-10)


def test_color_generator():
    gen = color_generator()
    preset_colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink',
                     'tab:gray', 'tab:olive', 'tab:cyan', 'gold', 'maroon', 'mediumslateblue', 'fuchsia',
                     'lawngreen', 'moccasin', 'thistle']
    for i in range(150):
        color = next(gen)
        assert (isinstance(color, str) and color in preset_colors) or (
            isinstance(color, np.ndarray) and len(color) == 3 and
            np.max(color) <= 1 and np.min(color) >= 0)


@pytest.mark.parametrize("this_set,other_sets,majority_threshold,truth",
                         [({1, 2, 3, 4}, [{1, 2, 3, 6}, {4, 5, 6}], 2 / 3, {1, 2, 3, 4, 6}),
                          ({'a', 'ab', 'aab'}, [{'ba', 'b'}], 0.501, set()),
                          ({'a', 'ab', 'aab'}, [{'ba', 'b'}], 0.5, {'a', 'ab', 'aab', 'ba', 'b'}),
                          ({1, 2, 3}, [{2, 3, 4}, {3, 4, 5}], 0.5, {2, 3, 4}),
                          ({1, 2, 3}, [{2, 3, 4}, {3, 4, 5}], 1, {3})])
def test_majority_vote_intersection(this_set, other_sets, majority_threshold, truth):
    result = SetWithMajorityVote.majority_vote_intersection(this_set, *other_sets,
                                                            majority_threshold=majority_threshold)
    assert result == truth


@pytest.mark.parametrize("is_df", [True, False])
@pytest.mark.parametrize("data,baseline,truth", [
    (np.array([1, 2, 3, 4, 5]), 0, np.array([0, 1, 2, 3, 4])),
    (np.array([[1, 2, 3], [-2, 4, 5], [0, 0, -1], [3, -2, 1]]), 1,
     np.array([[4, 5, 6], [1, 7, 8], [3, 3, 2], [6, 1, 4]])),
    (np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]), -1, np.array([[[-1, 0], [1, 2]], [[3, 4], [5, 6]]]))
])
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


@pytest.mark.parametrize("func,obj,truth", [
    (first_test_func, None, {}),
    (second_test_func, None, {'a': {'annotation': str, 'default': inspect._empty},
                              'b': {'annotation': bool, 'default': inspect._empty},
                              'c': {'annotation': inspect._empty, 'default': inspect._empty}}),
    (third_test_func, None, {'a': {'annotation': typing.List[str], 'default': inspect._empty},
                             'b': {'annotation': typing.Callable, 'default': inspect._empty},
                             'c': {'annotation': int, 'default': 3},
                             'd': {'annotation': inspect._empty, 'default': None}}),
    ('fourth_test_func', TestObj(), {'a': {'annotation': inspect._empty, 'default': inspect._empty},
                                     'b': {'annotation': None, 'default': inspect._empty},
                                     'c': {'annotation': float, 'default': 5.2}})])
def test_get_signature(func, obj, truth):
    this_signature = get_method_signature(func, obj)
    assert len(this_signature) == len(truth)
    for key, val in truth.items():
        assert key in this_signature
        param = this_signature[key]
        assert param.name == key
        assert param.annotation == val['annotation']
        assert param.default == val['default']


@pytest.mark.parametrize('intervals,expected', [
    ([], 0),
    ([(1, 3)], 3),
    ([(1, 3), (4, 6), (7, 10)], 10),
    ([(4, 6), (1, 3), (7, 10)], 10),
    ([(1, 4), (3, 6), (3, 5), (3, 6), (4, 9)], 9),
    ([(7, 10), (2, 5), (1, 4), (2, 5)], 9),

])
def test_sum_intervals_inclusive(intervals, expected):
    res = sum_intervals_inclusive(intervals)
    assert res == expected


@pytest.mark.parametrize('seconds,expected', [
    (13, '00:13'),
    (0, '00:00'),
    (59.34, '00:59'),
    (60.0, '01:00'),
    (60.2, '01:00'),
    (192.17, '03:12'),
    (13 * 60 + 1.5, '13:01'),
    (100 * 60, '100:00')])
def test_format_time(seconds, expected):
    res = format_time(seconds)
    assert res == expected


@pytest.mark.parametrize('name,expected', [
    ('aname', 'aname'),
    ('name123n', 'name123n'),
    ('camelCase', 'camelCase'),
    ('snake_case', 'snake_case'),
    ('name with spaces ', 'name_with_spaces'),
    ('1name of var', 'var_1name_of_var'),
    ('12345', 'var_12345'),
    ('%asdf&sign', '_asdf_sign'),
    (' ^^more things123 ', '___more_things123')
])
def test_sanitize_variable_name(name, expected):
    res = sanitize_variable_name(name)
    assert res == expected


def test_get_method_readable_name():
    def func(x):
        return x + 1
    func.readable_name = "readable name"
    assert get_method_readable_name(func) == "readable name"

    def func2(a, b):
        return a + b

    assert get_method_readable_name(func2) == 'func2'

    func2.readable_name = "readable name"
    assert get_method_readable_name(func2) == "readable name"


def test_param_readable_names_decorator_sets_attribute():
    mapping = {'trim_n': 'Trim ambiguous (N) bases', 'drop_columns': 'Columns to drop'}

    @param_readable_names(mapping)
    def func(trim_n=True, drop_columns=None):
        return trim_n, drop_columns

    assert func.readable_param_names == mapping
    # the decorator must not alter the function's behavior
    assert func(False, ['a']) == (False, ['a'])


@pytest.mark.parametrize("param_name,expected", [
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
])
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
    "X, labels, expected_bic",
    [
        (
            np.array([[1, 2], [1.5, 2.5], [3, 4], [4, 5]]),
            np.array([0, 0, 1, 1]),
            -13.510391348997054
        ),
        (
            np.array([[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]]),
            np.array([0, 0, 0, 1, 1]),
            -23.462198946075148
        ),
        (
            np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4], [5, 5, 5]]),
            np.array([0, 0, 0, 1, 1]),
            -33.747038606183764
        ),
    ]
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
        event = PickEvent('pick_event', self.fig.canvas, mouseevent=MockMouseEvent(), artist=self.fig.ax.collections[0],
                          ind=[0])
        self.fig.on_pick(event)
        assert 0 in self.fig.is_labeled
        assert isinstance(self.fig.is_labeled[0], plt.Annotation)
        assert self.fig.is_labeled[0].get_text() == 'A'

    def test_on_pick_remove_annotation(self):
        event = PickEvent('pick_event', self.fig.canvas, mouseevent=MockMouseEvent(), artist=self.fig.ax.collections[0],
                          ind=[0])
        self.fig.on_pick(event)
        self.fig.on_pick(event)
        assert 0 not in self.fig.is_labeled

    def test_on_pick_invalid_artist(self):
        event = PickEvent('pick_event', self.fig.canvas, mouseevent=MockMouseEvent(),
                          artist=plt.Rectangle((0, 0), 1, 1), ind=[0])
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
