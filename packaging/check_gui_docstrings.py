"""
Pre-commit hook (``check-gui-docstrings``): every parameter the GUI renders must have help text.

The GUI is generated from the public API by reflection, and each parameter's help text/tooltip comes from its
``:param name: ...`` docstring line, parsed by ``rnalysis/utils/docstring_parsing.py`` (see CLAUDE.md, "the API<->GUI
reflection contract"). A typo in that line silently leaves the parameter with no help in the GUI. This hook catches
that at commit time, for ``rnalysis/filtering.py``, ``rnalysis/enrichment.py`` and ``rnalysis/fastq.py``, by checking
each GUI-exposed function for:

- ``missing-help``: a signature parameter (other than self/cls/*args/**kwargs) whose parsed help text is missing/empty;
- ``unknown-param``: a ``:param name:`` line naming a parameter that is not in the signature (usually a typo);
- ``malformed-field``: a line that looks like a ``:param`` field but does not parse as one (e.g. ``param x:``);
- ``no-blank-line``: a ``:param`` field inside the description paragraph (no blank line before the ``:param`` block),
  which the parser swallows into the description;
- ``no-docstring``: a ``@readable_name``-decorated function with parameters but no docstring at all.

pre-commit runs this in an isolated environment without RNAlysis's dependencies (polars, numpy, PyQt6...), so it
never imports the ``rnalysis`` package: modules are read with ``ast``, and the shared, stdlib-only parser is loaded
by file path. Being static, the check approximates what the GUI does at runtime:

- GUI-exposed functions mirror ``tests/test_parsing.py::_gui_exposed_functions``: public methods of ``Filter`` /
  ``FeatureSet`` and their subclasses (``gui.TabPage.get_all_actions``), and public module-level functions of
  ``fastq`` and ``enrichment``. ``EXCLUDED_FUNCS`` are still checked (a stricter superset; they are public API whose
  docstrings also feed the API docs). Properties are skipped (``get_all_actions`` only keeps callables).
- For methods, parameters in ``GUI_HIDDEN_METHOD_PARAMS`` are skipped: the tab pages never render them
  (``gui.FuncTypeStack.EXCLUDED_PARAMS``). Module-level functions are rendered by ``FuncExternalWindow`` subclasses,
  most of which render every parameter, so only ``GUI_HIDDEN_FUNCTION_PARAMS`` is skipped for them. Other
  per-window exclusions are not mirrored; a parameter they hide but that lacks help text is listed in
  ``KNOWN_ISSUES`` instead.
- Inherited docstrings (``inspect.getdoc``) are resolved through base classes defined in the same file. A method
  with no docstring whose class derives from a class defined elsewhere (other than a builtin such as ``set``, which
  is checked directly) is skipped when that base may define the method, since its docstring may be inherited from
  there.
- ``no-docstring`` is only reported for ``@readable_name``-decorated functions. The GUI also shows undecorated
  public methods under their function name; an undocumented one is not reported.
- ``KNOWN_ISSUES`` is a baseline of problems that predate this hook. They are tolerated, but an entry that no longer
  occurs is reported, so the baseline can only shrink.

Usage: ``python packaging/check_gui_docstrings.py [FILE ...]`` (pre-commit passes the changed files). If this script
or the shared parser is among them, all checked modules are checked.
"""

import ast
import builtins
import importlib.util
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
PARSER_PATH = 'rnalysis/utils/docstring_parsing.py'
SCRIPT_PATH = 'packaging/check_gui_docstrings.py'

# module path -> (check public module-level functions?, root classes whose public methods (and subclasses') are checked)
CHECKED_MODULES: Dict[str, Tuple[bool, Tuple[str, ...]]] = {
    'rnalysis/filtering.py': (False, ('Filter',)),
    'rnalysis/enrichment.py': (True, ('FeatureSet',)),
    'rnalysis/fastq.py': (True, ()),
}

# parameters the GUI never renders for Filter/FeatureSet methods (see gui.FuncTypeStack.EXCLUDED_PARAMS)
GUI_HIDDEN_METHOD_PARAMS = frozenset({'self', 'cls', 'backend', 'gui_mode', 'legacy_args', 'function_kwargs'})
# ... and for module-level fastq/enrichment functions (excluded by every FuncExternalWindow that could show it)
GUI_HIDDEN_FUNCTION_PARAMS = frozenset({'legacy_args'})

MISSING = 'missing-help'
UNKNOWN = 'unknown-param'
MALFORMED = 'malformed-field'
NO_BLANK_LINE = 'no-blank-line'
NO_DOCSTRING = 'no-docstring'

MESSAGES = {
    MISSING: 'has no help text in the docstring (add a ":param {param}: ..." line)',
    UNKNOWN: 'is documented with ":param {param}:" but is not in the signature (typo?)',
    MALFORMED: 'has a malformed ":param" line that the GUI cannot parse (expected ":param {param}: ...")',
    NO_BLANK_LINE: 'docstring has no blank line between the description and the ":param" block, so the GUI shows no '
    'parameter help',
    NO_DOCSTRING: 'is shown in the GUI but has no docstring, so its parameters have no help text',
}

# Problems that predate this hook: (module, qualname, param, kind). Fix them and delete the entry.
KNOWN_ISSUES: Set[Tuple[str, str, Optional[str], str]] = {
    ('rnalysis/enrichment.py', 'FeatureSet.biotypes_from_gtf', 'long_format', UNKNOWN),
    ('rnalysis/enrichment.py', 'FeatureSet.go_enrichment', 'biotype', UNKNOWN),
    ('rnalysis/enrichment.py', 'FeatureSet.go_enrichment', 'biotype_ref_path', UNKNOWN),
    ('rnalysis/enrichment.py', 'FeatureSet.go_enrichment', 'random_seed', MISSING),
    ('rnalysis/enrichment.py', 'FeatureSet.kegg_enrichment', 'random_seed', MISSING),
    ('rnalysis/enrichment.py', 'FeatureSet.map_orthologs_phylomedb', 'consistency_score_threshold', MISSING),
    ('rnalysis/enrichment.py', 'FeatureSet.translate_gene_ids', 'inplace', MISSING),
    ('rnalysis/enrichment.py', 'FeatureSet.user_defined_enrichment', 'biotype', UNKNOWN),
    ('rnalysis/enrichment.py', 'FeatureSet.user_defined_enrichment', 'biotype_ref_path', UNKNOWN),
    ('rnalysis/enrichment.py', 'FeatureSet.user_defined_enrichment', 'random_seed', MISSING),
    ('rnalysis/fastq.py', 'create_bam_index', None, NO_DOCSTRING),
    ('rnalysis/fastq.py', 'fastq_to_sam_paired', 'input_folder', UNKNOWN),
    ('rnalysis/fastq.py', 'fastq_to_sam_paired', 'output_format', MISSING),
    ('rnalysis/fastq.py', 'fastq_to_sam_paired', 'quality_score_type', MISSING),
    ('rnalysis/fastq.py', 'fastq_to_sam_paired', 'r1_files', MISSING),
    ('rnalysis/fastq.py', 'fastq_to_sam_paired', 'r2_files', MISSING),
    ('rnalysis/fastq.py', 'fastq_to_sam_single', 'output_format', MISSING),
    ('rnalysis/fastq.py', 'fastq_to_sam_single', 'quality_score_type', MISSING),
    ('rnalysis/fastq.py', 'kallisto_quantify_paired_end', 'learn_bias', UNKNOWN),
    ('rnalysis/fastq.py', 'kallisto_quantify_paired_end', 'seek_fusion_genes', UNKNOWN),
    ('rnalysis/fastq.py', 'kallisto_quantify_single_end', 'learn_bias', UNKNOWN),
    ('rnalysis/fastq.py', 'kallisto_quantify_single_end', 'seek_fusion_genes', UNKNOWN),
    ('rnalysis/fastq.py', 'sam_to_fastq_paired', 're_reverse_reads', MISSING),
    ('rnalysis/fastq.py', 'sam_to_fastq_paired', 'return_new_filenames', MISSING),
    ('rnalysis/fastq.py', 'trim_adapters_paired_end', 'return_new_filenames', MISSING),
    ('rnalysis/fastq.py', 'validate_sam', 'new_sample_names', UNKNOWN),
    ('rnalysis/filtering.py', 'CountFilter.differential_expression_deseq2', 'cooks_cutoff', MISSING),
    ('rnalysis/filtering.py', 'CountFilter.differential_expression_deseq2', 'scaling_factors', MISSING),
    ('rnalysis/filtering.py', 'CountFilter.differential_expression_deseq2_simplified', 'return_log', MISSING),
    ('rnalysis/filtering.py', 'CountFilter.differential_expression_limma_voom_simplified', 'return_log', MISSING),
    ('rnalysis/filtering.py', 'Filter.concatenate', None, NO_DOCSTRING),
    ('rnalysis/filtering.py', 'Filter.histogram', None, NO_DOCSTRING),
    ('rnalysis/filtering.py', 'Filter.map_orthologs_phylomedb', 'consistency_score_threshold', MISSING),
    ('rnalysis/filtering.py', 'Filter.split_by_percentile', 'interpolate', MISSING),
    ('rnalysis/filtering.py', 'FoldChangeFilter.randomization_test', 'random_seed', MISSING),
}

_STRICT_PARAM = re.compile(r'^[ \t]*:param ([a-zA-Z_0-9]+):', re.MULTILINE)
_PARAM_LIKE = re.compile(r'^[ \t]*(?::param\b|param\s+\w+\s*:)')
_PARAM_LIKE_NAME = re.compile(r'param\s+(\w+)')
_FIELD_HEAD = re.compile(r'^[ \t]*:param\s+([^:]*?)\s*:')


def _load_parser():
    # loaded by path, so that rnalysis/__init__.py (which imports the heavy dependencies) is never executed
    spec = importlib.util.spec_from_file_location('rnalysis_docstring_parsing', REPO_ROOT / PARSER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


parse_docstring = _load_parser().parse_docstring


@dataclass(frozen=True)
class Problem:
    path: str
    line: int
    qualname: str
    param: Optional[str]
    kind: str

    def key(self) -> Tuple[str, str, Optional[str], str]:
        return self.path, self.qualname, self.param, self.kind

    def __str__(self) -> str:
        subject = f'parameter {self.param!r}' if self.param is not None else 'function'
        return f'{self.path}:{self.line}: {self.qualname}: {subject} {MESSAGES[self.kind].format(param=self.param)}'


def _decorator_names(func: ast.AST) -> Set[str]:
    names = set()
    for dec in func.decorator_list:
        node = dec.func if isinstance(dec, ast.Call) else dec
        if isinstance(node, ast.Attribute):
            names.add(node.attr)
        elif isinstance(node, ast.Name):
            names.add(node.id)
    return names


def _is_function(node: ast.AST) -> bool:
    return isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))


def _signature_params(func: ast.FunctionDef, is_method: bool) -> Tuple[List[str], Set[str]]:
    """Return the parameters to check, and every name in the signature (including self, *args and **kwargs)."""
    args = func.args
    named = [arg.arg for arg in args.posonlyargs + args.args + args.kwonlyargs]
    if is_method and 'staticmethod' not in _decorator_names(func) and named:
        named = named[1:]  # self / cls
    all_names = set(named) | {arg.arg for arg in (args.vararg, args.kwarg) if arg is not None}
    hidden = GUI_HIDDEN_METHOD_PARAMS if is_method else GUI_HIDDEN_FUNCTION_PARAMS
    return [name for name in named if name not in hidden], all_names


def _find_line(lines: Sequence[str], func: ast.FunctionDef, pattern: str) -> int:
    for lineno in range(func.lineno, (func.end_lineno or func.lineno) + 1):
        if re.search(pattern, lines[lineno - 1]):
            return lineno
    return func.lineno


class _ClassIndex:
    def __init__(self, tree: ast.Module):
        self.classes: Dict[str, ast.ClassDef] = {
            node.name: node for node in tree.body if isinstance(node, ast.ClassDef)
        }

    def _bases(self, cls: ast.ClassDef) -> List[ast.expr]:
        return [base for base in cls.bases if not (isinstance(base, ast.Name) and base.id == 'object')]

    def local_base(self, base: ast.expr) -> Optional[ast.ClassDef]:
        return self.classes.get(base.id) if isinstance(base, ast.Name) else None

    def mro(self, cls: ast.ClassDef) -> Tuple[List[ast.ClassDef], List[ast.expr]]:
        """Same-file ancestors (depth-first, including cls), and the bases defined elsewhere."""
        order, external, stack = [], [], [cls]
        while stack:
            current = stack.pop(0)
            if current in order:
                continue
            order.append(current)
            for base in self._bases(current):
                local = self.local_base(base)
                if local is None:
                    external.append(base)
                else:
                    stack.append(local)
        return order, external

    def is_gui_class(self, cls: ast.ClassDef, roots: Iterable[str]) -> bool:
        return any(ancestor.name in roots for ancestor in self.mro(cls)[0])

    def inherited_docstring(self, cls: ast.ClassDef, name: str) -> Tuple[Optional[str], bool]:
        """
        Mimic inspect.getdoc(): the first docstring of method `name` along the same-file MRO. Also return whether a
        base defined elsewhere may define `name` (and so supply a docstring this static check cannot see):
        builtin bases (e.g. `set`) are checked directly, any other external base is assumed to.
        """
        order, external = self.mro(cls)
        for ancestor in order:
            for node in ancestor.body:
                if _is_function(node) and node.name == name:
                    doc = ast.get_docstring(node)
                    if doc is not None:
                        return doc, False
        may_inherit = False
        for base in external:
            builtin = getattr(builtins, base.id, None) if isinstance(base, ast.Name) else None
            if not isinstance(builtin, type) or hasattr(builtin, name):
                may_inherit = True
        return None, may_inherit


def _check_function(
    func: ast.FunctionDef, qualname: str, doc: Optional[str], is_method: bool, path: str, lines: Sequence[str]
) -> List[Problem]:
    params, all_names = _signature_params(func, is_method)
    if doc is None:
        if params and 'readable_name' in _decorator_names(func):
            return [Problem(path, func.lineno, qualname, None, NO_DOCSTRING)]
        return []

    problems = []
    desc, parsed = parse_docstring(doc)
    if _STRICT_PARAM.search(desc):
        return [Problem(path, _find_line(lines, func, r':param '), qualname, None, NO_BLANK_LINE)]

    for line in doc.splitlines():
        collapsed = re.sub(' +', ' ', line)
        if _PARAM_LIKE.match(collapsed) and not _STRICT_PARAM.match(collapsed):
            head = _FIELD_HEAD.match(collapsed)
            if head:  # ':param int x: ...' / ':param x : ...' -> the name is the last word before the colon
                words = re.findall(r'\w+', head.group(1))
                name = words[-1] if words else None
            else:
                match = _PARAM_LIKE_NAME.search(collapsed)
                name = match.group(1) if match else None
                if not collapsed.lstrip().startswith(':') and name not in all_names:
                    continue  # a prose line such as 'param values: ...', not a broken field
            problems.append(Problem(path, _find_line(lines, func, re.escape(line.strip())), qualname, name, MALFORMED))
    malformed = {problem.param for problem in problems}

    for name in params:
        if not parsed.get(name, '').strip() and name not in malformed:
            line = _find_line(lines, func, rf':param {name}:|\b{name}\b')
            problems.append(Problem(path, line, qualname, name, MISSING))

    for name in _STRICT_PARAM.findall(doc):
        if name not in all_names:
            problems.append(Problem(path, _find_line(lines, func, rf':param {name}:'), qualname, name, UNKNOWN))
    return problems


def check_source(source: str, path: str, module_functions: bool, root_classes: Sequence[str]) -> List[Problem]:
    """Check the GUI-exposed functions of one module's source code."""
    tree = ast.parse(source)
    lines = source.splitlines()
    index = _ClassIndex(tree)
    problems = []
    for node in tree.body:
        if module_functions and _is_function(node) and not node.name.startswith('_'):
            problems.extend(_check_function(node, node.name, ast.get_docstring(node), False, path, lines))
        elif (
            isinstance(node, ast.ClassDef) and not node.name.startswith('_') and index.is_gui_class(node, root_classes)
        ):
            for method in node.body:
                if not _is_function(method) or method.name.startswith('_'):
                    continue
                if _decorator_names(method) & {'property', 'setter', 'getter', 'deleter', 'cached_property'}:
                    continue
                doc = ast.get_docstring(method)
                if doc is None:
                    doc, may_inherit = index.inherited_docstring(node, method.name)
                    if doc is None and may_inherit:
                        continue  # the docstring may be inherited from a class the static check cannot see
                qualname = f'{node.name}.{method.name}'
                problems.extend(_check_function(method, qualname, doc, True, path, lines))
    return problems


def check_file(path: str) -> List[Problem]:
    module_functions, root_classes = CHECKED_MODULES[path]
    source = (REPO_ROOT / path).read_text(encoding='utf-8')
    return check_source(source, path, module_functions, root_classes)


def _normalize(filename: str) -> str:
    path = Path(filename)
    if path.is_absolute():
        try:
            path = path.resolve().relative_to(REPO_ROOT)
        except ValueError:
            pass
    return path.as_posix()


def main(argv: Sequence[str]) -> int:
    changed = {_normalize(filename) for filename in argv}
    if changed & {PARSER_PATH, SCRIPT_PATH}:
        to_check = sorted(CHECKED_MODULES)
    else:
        to_check = sorted(changed & set(CHECKED_MODULES))

    failures = []
    for path in to_check:
        problems = check_file(path)
        found = {problem.key() for problem in problems}
        failures.extend(str(problem) for problem in problems if problem.key() not in KNOWN_ISSUES)
        for stale in sorted(key for key in KNOWN_ISSUES if key[0] == path and key not in found):
            failures.append(
                f'{path}: {stale[1]}: {stale[2]!r} ({stale[3]}) is fixed - remove it from KNOWN_ISSUES in {SCRIPT_PATH}'
            )

    if failures:
        print('GUI help-text problems (see the check-gui-docstrings hook in .pre-commit-config.yaml):')
        print('\n'.join(failures))
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
