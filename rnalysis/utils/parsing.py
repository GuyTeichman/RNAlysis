import difflib
import itertools
import platform
import re
import shlex
import warnings
from collections import Counter
from itertools import islice
from pathlib import Path
from typing import Any, Dict, Union, List, Tuple, Iterable, Literal

import mslex
import numpy as np
import pandas as pd
import polars as pl
import unicodedata

from rnalysis.utils import validation


def quote_path(pth: Union[Path, str]) -> str:
    if isinstance(pth, Path):
        pth = pth.as_posix()

    if platform.system() == 'Windows':
        return mslex.quote(pth)
    else:
        return shlex.quote(pth)


def r_make_names(names: List[str]):
    # Replace invalid characters with a dot
    valid_names = [re.sub(r'[^a-zA-Z0-9_]', '.', name) for name in names]
    # Ensure names start with a letter or a dot not followed by a number
    valid_names = [name if re.match(r'^[a-zA-Z].*|^\.[^0-9]', name) else '.' + name for name in valid_names]
    # Check for duplicates and append numbers if necessary
    counts = Counter(valid_names)
    name_map = {}
    for name, count in counts.items():
        if count > 1:  # Duplicate exists
            suffix_count = 1
            for i in range(len(valid_names)):
                if valid_names[i] == name:
                    new_name = name
                    # Append a number if this is not the first occurrence
                    if name in name_map:
                        new_name = f"{name}.{suffix_count}"
                        suffix_count += 1
                    valid_names[i] = new_name
                    name_map[name] = True  # Mark as processed
        else:
            name_map[name] = True  # Mark as unique
    return valid_names


def python_to_r_kwargs(kwargs: dict, delimiter: str = ',\n'):
    kwargs_str = ''
    for key, val in kwargs.items():
        val = _parse_r_arg(val)
        this_arg = f'{key} = {val}{delimiter}'
        kwargs_str += this_arg
    return kwargs_str[:-len(delimiter)]  # ignore the last delimiter


def _parse_r_arg(arg):
    if isinstance(arg, bool):
        return str(arg).upper()
    elif isinstance(arg, (int, float)):
        return str(arg)
    elif isinstance(arg, str):
        return f'"{arg}"'
    elif arg is None:
        return 'NULL'
    elif validation.isinstanceiter(arg, str):
        base = ', '.join([f'"{item}"' for item in arg])
        return f'c({base})'
    else:
        raise TypeError(f"Cannot parse argument: {arg} of type: {type(arg)}")


def remove_suffixes(path: Union[str, Path]) -> Union[str, Path]:
    edited_path = Path(path)
    while edited_path.suffix:
        edited_path = edited_path.with_suffix('')
    if isinstance(path, Path):
        return edited_path
    return edited_path.as_posix()


def from_string(msg: str = '', del_spaces: bool = False, delimiter: str = '\n'):
    """
    Takes a manual string input from the user, and then splits it using a comma delimiter into a list of values. \
    Called when an FeatureSet instance is created without input, \
    or when FeatureSet.enrich_randomization is called without input.

    :param msg: a promprt to be printed to the user
    :param del_spaces: if True, will delete all spaces in each delimited value.
    :param delimiter: the delimiter used to separate the values. Default is '\n'
    :return: A list of the comma-seperated values the user inserted.

    """
    string = input(msg)
    split = string.split(sep=delimiter)
    if del_spaces:
        for i in range(len(split)):
            split[i] = split[i].replace(' ', '')
    if split[-1] == '':
        split = split[:-1]
    return split


def uniprot_tab_to_dict(tab_input: str) -> Tuple[Dict[str, str], list]:
    split_list = tab_input.split()
    if len(split_list) == 0:
        return {}, []
    split_list = split_list[2:]

    parsed = {}
    duplicates = []
    for key, val in zip(islice(split_list, 0, len(split_list), 2), islice(split_list, 1, len(split_list), 2)):
        if key in parsed:
            parsed[key].append(val)
        else:
            parsed[key] = [val]

    for key in list(parsed.keys()):
        if len(parsed[key]) > 1:
            duplicates.extend(parsed.pop(key))
        else:
            parsed[key] = parsed[key][0]

    return parsed, duplicates


def uniprot_tab_with_score_to_dict(tab_input: str, reverse_key_value: bool = False) -> Dict[str, str]:
    split_list = re.split('[\t\n]', tab_input)
    if len(split_list) == 0:
        return {}
    split_list = split_list[3:]

    parsed: Dict[str, Tuple[List[str], List[int]]] = {}
    for gene_id, rank, key in zip(islice(split_list, 0, len(split_list), 3), islice(split_list, 1, len(split_list), 3),
                                  islice(split_list, 2, len(split_list), 3)):
        if reverse_key_value:
            gene_id, key = key, gene_id
        numeric_rank = int(rank[0])
        if key in parsed:
            parsed[key][0].append(gene_id)
            parsed[key][1].append(numeric_rank)
        else:
            parsed[key] = ([gene_id], [numeric_rank])

    key_to_id = {}
    duplicates = {}
    for key, (gene_ids, ranks) in zip(parsed.keys(), parsed.values()):
        if len(gene_ids) == 1:
            key_to_id[key] = gene_ids[0]
        else:
            best_id = gene_ids[max(range(len(ranks)), key=lambda i: ranks[i])]
            key_to_id[key] = best_id
            duplicates[key] = best_id

    if len(duplicates) > 0:
        warnings.warn(f"Duplicate mappings were found for {len(duplicates)} genes. "
                      f"The following mapping was chosen for them based on their annotation score: {duplicates}")
    return key_to_id


def data_to_list(data: Any, sort: bool = False) -> list:
    if isinstance(data, list):
        lst = data
    elif isinstance(data, (set, tuple, np.ndarray)):
        lst = list(data)
    elif isinstance(data, (dict, int, float, bool, str)):
        lst = [data]
    elif data is None:
        lst = [None]
    elif isinstance(data, pl.DataFrame):
        return data.to_series().to_list()
    elif isinstance(data, pl.Series):
        return data.to_list()
    elif callable(data):
        lst = [data]
    else:
        try:
            lst = list(data)
        except TypeError:
            lst = [data]

    if sort:
        lst.sort()
    return lst


def data_to_tuple(data: Any, sort: bool = False) -> tuple:
    if isinstance(data, tuple):
        tpl = data
    elif isinstance(data, (set, list, np.ndarray)):
        tpl = tuple(data)
    elif isinstance(data, (int, float, bool, str)):
        tpl = data,
    elif data is None:
        tpl = None,
    elif isinstance(data, pl.DataFrame):
        return tuple(data.to_series().to_list())
    elif isinstance(data, pl.Series):
        return tuple(data.to_list())
    elif callable(data):
        tpl = data,
    else:
        try:
            tpl = tuple(data)
        except TypeError:
            tpl = (data,)
    if sort:
        tpl = tuple(sorted(tpl))
    return tpl


def data_to_set(data: Any) -> set:
    if isinstance(data, set):
        return data
    elif isinstance(data, (list, tuple, np.ndarray)):
        return set(data)
    elif isinstance(data, (int, float, bool, str)):
        return {data}
    elif isinstance(data, pl.DataFrame):
        return set(data.to_series().to_list())
    elif isinstance(data, pl.Series):
        return set(data.to_list())
    elif data is None:
        return {None}
    elif callable(data):
        return {data}
    else:
        try:
            return set(data)
        except TypeError:
            return {data}


def partition_list(lst: Union[list, tuple], chunk_size: int) -> Union[List[tuple], List[list]]:
    assert isinstance(lst, (list, tuple)), f"'lst' must be a list or tuple; instead got type {type(lst)}."
    assert isinstance(chunk_size, int), f"'chunk_size' must be an integer; instead got type {type(chunk_size)}."
    assert chunk_size > 0, f"'chunk_size' must be >0; instead got {chunk_size}."
    if len(lst) == 0:
        return [type(lst)()]
    return [lst[i: i + chunk_size] for i in range(0, len(lst), chunk_size)]


def flatten(lst: list) -> list:
    """
    Flatten a list of arbitrary depth.

    :param lst: the list to be flattened
    :type lst: list
    :return: a flattened list
    :rtype: list
    """
    output = []
    for item in lst:
        if isinstance(item, list):
            output.extend(flatten(item))
        else:
            output.append(item)
    return output


def parse_docstring(docstring: str) -> Tuple[str, Dict[str, str]]:
    """
    Parse a given docstring (str) to retreive the description text, as well as a dictionary of parameter descriptions.

    :param docstring: the docstring to be parsed
    :type docstring: str
    :return: a string matching the description, and a dictionary containing the parameter descriptions
    """
    docstring = re.sub(' +', ' ', docstring)
    split = docstring.split('\n\n')
    desc = split[0]
    params_str = split[1]
    free_text_match = '["#\w\s\.\(\)\-\:%\*=@!\?\+\_,/`><\[\]' + "'" + ']'
    params_matches = list(
        re.finditer('^:param ([a-zA-Z_0-9]+):(' + free_text_match + '+?)(?=^\:.*\:)', params_str, re.MULTILINE))
    params = {match.group(1): match.group(2).replace('. ', '. \n') for match in params_matches}
    return desc, params


def generate_upset_series(objs: dict):
    """
    Receives a dictionary of sets from enrichment._fetch_sets(), \
    and reformats it as a pandas Series to be used by the python package 'upsetplot'.

    :param objs: the output of the enrichment._fetch_sets() function.
    :type objs: dict of sets
    :return: a pandas Series in the format requested by the 'upsetplot' package.

    """
    names = list(objs.keys())
    multi_ind = pd.MultiIndex.from_product([[True, False] for _ in range(len(names))], names=names)[:-1]
    srs = pd.Series(index=multi_ind, dtype='uint32')
    for ind in multi_ind:
        intersection_sets = list(itertools.compress(names, ind))
        difference_sets = list(itertools.compress(names, (not i for i in ind)))
        group = set.intersection(*[objs[s] for s in intersection_sets]).difference(*[objs[s] for s in difference_sets])
        group_size = len(group)
        srs.loc[ind] = group_size
    return srs


def parse_version(version: str):
    split = version.split('.')
    return [int(i) for i in split]


def parse_gtf_attributes(attr_str: str):
    attributes_dict = {}
    for this_attr in attr_str.split('; '):
        this_attr = this_attr.strip()
        key_end = this_attr.find(' ')
        if key_end == -1:
            continue
        key = this_attr[:key_end]
        val_start = this_attr.find('"', key_end)
        val_end = this_attr.find('"', val_start + 1)
        val = this_attr[val_start + 1:val_end]
        attributes_dict[key] = val
    return attributes_dict


def parse_gff3_attributes(attr_str: str):
    attributes_dict = {}
    for this_attr in attr_str.rstrip().split(';'):
        this_attr = this_attr.strip()
        if len(this_attr) == 0:
            continue
        key, val = this_attr.split('=')
        val = val.split(',')
        if len(val) == 1:
            val = val[0]
        attributes_dict[key] = val
    return attributes_dict


def format_dict_for_display(d: dict) -> str:
    """
    Formats a dictionary for display in a human-readable format.

    :param d: The dictionary to format.
    :type d: dict
    :return: A string representation of the dictionary with keys and values formatted for readability.
    :rtype: str
    """
    output = []
    for key, val in d.items():
        output.append(f"{key} = {repr(val)}")
    return ', \n'.join(output)


def replace_last_occurrence(regex, repl, item):
    """
    Replaces the last occurrence of a given regular expression in a string with a replacement string.

    :param regex: The regular expression to search for.
    :type regex: str
    :param repl: The replacement string.
    :type repl: str
    :param item: The input string.
    :type item: str
    :return: The modified string.
    :rtype: str
    """
    match = None
    for m in re.finditer(regex, item):
        match = m
    if match is None:
        # If no match was found, return the original string
        return item
    else:
        # Replace the last occurrence of the regex with the replacement string
        start, end = match.span()
        return item[:start] + re.sub(regex, repl, item[start:end]) + item[end:]


def items_to_html_table(items):
    """
    Converts a list of strings to an HTML table with a single column.

    :param items: The list of strings to convert.
    :type items: list of str
    :return: A string representation of the HTML table.
    :rtype: str
        """
    table = '<table border="1" class="dataframe">\n'
    for item in items:
        table += f'<tr><td style="border: 1px solid black; border-collapse: collapse;"><b>{item}</b></td></tr>\n'
    table += "</table>"
    return table


def df_to_html(df: pl.DataFrame, max_rows: int = 5, max_cols: int = 5):
    """
    Convert a Polars DataFrame to an HTML table.

    :param max_cols: maximum number of columns to display in the HTML table.
    :type max_cols: int (default=4)
    :param max_rows: maximum number of rows to display in the HTML table.
    :type max_rows: int (default=5)
    :param df: A Polars DataFrame to be converted.
    :type df: polars.DataFrame
    :return: A string representation of the HTML table.
    :rtype: str
    """
    if isinstance(df, pl.Series):
        df = df.to_frame()
    if df.shape[-1] == 0:
        return ''

    # Limit the number of rows and columns
    df_subset = df.head(max_rows).select(pl.col(df.columns[:max_cols])).to_pandas()
    df_subset = df_subset.set_index(df_subset.columns[0])
    df_subset.index.names = [None]

    styler = df_subset.style.format(precision=2)
    styler.set_table_styles(
        [{'selector': 'td', 'props': 'border: 1px solid grey; border-collapse: collapse;'},
         {'selector': 'th', 'props': 'border: 1px solid grey; border-collapse: collapse;'}], )
    html = styler.to_html(float_format=lambda x: f"{x:.2f}")
    if df.shape[0] > max_rows and df.shape[1] > max_cols:
        # remove a redundant '...' from the end of the table
        html = replace_last_occurrence(r'<td class="data col[\d]+ row_trim" >...<\/td>', '', html)
    return html


def parse_uniprot_id(input_string):
    """
    Parses a UniProt ID from a string.

    :param input_string: The string to parse.
    :type input_string: str
    """
    pattern = r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
    match = re.search(pattern, input_string)
    if match:
        return match.group()


def slugify(value, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize("NFKC", value)
    else:
        value = (
            unicodedata.normalize("NFKD", value)
            .encode("ascii", "ignore")
            .decode("ascii")
        )
    value = re.sub(r"[^\w\s-]", "", value)
    return re.sub(r"[-\s]+", "-", value).strip("-_")


def longest_common_substring(s1, s2):
    """Return the longest common substring between two strings using difflib."""
    matcher = difflib.SequenceMatcher(None, s1, s2)
    match = matcher.find_longest_match(0, len(s1), 0, len(s2))
    return s1[match.a: match.a + match.size]


def common_suffix(strings):
    """Find the common suffix among a list of strings."""
    if not strings:
        return ''
    if len(strings) == 1:
        return strings[0]
    s1 = min(strings)
    s2 = max(strings)
    for i, c in enumerate(reversed(s1)):
        if c != s2[-(i + 1)]:
            return s1[len(s1) - i:]
    return s1


def remove_suffix(s: str, suffix: Union[str, List[str]]):
    """Remove a suffix from a string."""
    suffixes = sorted(data_to_list(suffix), key=len, reverse=True)
    for this_suffix in suffixes:
        if s.endswith(this_suffix):
            return s[:-len(this_suffix)]
    return s


def generate_common_name(file_pairs):
    suffixes = ['_trimmed', '_sorted', '_']
    file_pairs = [(remove_suffix(p1, suffixes), remove_suffix(p2, suffixes)) for p1, p2 in file_pairs]
    lcs_list = [(pair[0], pair[1], longest_common_substring(pair[0], pair[1])) for pair in file_pairs]

    # Find LCS among all LCSs
    overall_lcs = lcs_list[0][2] if len(lcs_list) > 1 else ""
    for _, _, lcs in lcs_list[1:]:
        overall_lcs = longest_common_substring(overall_lcs, lcs)

    # Prepare output based on comparison
    initial_results = []
    for s1, s2, lcs in lcs_list:
        if len(lcs) > len(overall_lcs):
            initial_results.append(lcs)
        else:
            initial_results.append(s1 + s2)
    # Find and trim common suffix from LCSs
    lcs_only = [result for result in initial_results if len(result) <= len(s1 + s2)]
    suffix = common_suffix(lcs_only)
    if suffix and len(lcs_only) > 1:
        trimmed_results = [
            result[:-len(suffix)] if result.endswith(suffix) and len(result[:-len(suffix)]) > 0 else result for result
            in initial_results]
    else:
        trimmed_results = initial_results

    return trimmed_results


def make_group_names(sample_grouping: List[Iterable[str]], mode: Literal['auto', 'display'] = 'display'):
    assert mode in ['auto', 'display'], f"Invalid mode: {mode}. Must be 'auto' or 'display'."
    sep = '\n' if mode == 'display' else ';'
    group_names = []
    for i, group in enumerate(sample_grouping):
        if isinstance(group, str):
            group_names.append(group)
        elif validation.isiterable(group):
            group_names.append(sep.join(group))
    return group_names
