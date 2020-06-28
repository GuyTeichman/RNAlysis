import json
import os
import warnings
from typing import List, Set, Union

import requests

from rnalysis.utils import parsing, validation, __path__


def map_ids(ids: Union[str, Set[str], List[str]], map_from: str, map_to: str = 'UniProtKB AC'):
    url = 'https://www.uniprot.org/uploadlists/'
    id_dict = load_id_abbreviation_dict()
    validation.validate_uniprot_dataset_name(id_dict, map_to, map_from)

    params = {
        'from': id_dict[map_from],
        'to': id_dict[map_to],
        'format': 'tab',
        'query': format_ids(ids)
    }
    print(f"Mapping {len(params['query'])} entries from '{map_from}' to '{map_to}'...")
    req = requests.get(url, params=params)
    if req.status_code != 200:
        raise ConnectionError(f"Request failed with status code {req.status_code}.")
    output = parsing.uniprot_tab_to_dict(req.text)
    if len(output) < len(params['query']):
        warnings.warn(f"Failed to map {len(output) - len(params['query'])} entries from '{map_from}' to '{map_to}'. "
                      f"Returning the remaining {len(output)} entries.")
    return output


def format_ids(ids: Union[str, int, list, set]):
    if isinstance(ids, str):
        return ids
    elif isinstance(ids, int):
        return str(ids)
    return " ".join((str(item) for item in ids))


def load_id_abbreviation_dict(dict_path: str = os.path.join(__path__[0], 'uniprot_dataset_abbreviation_dict.json')):
    with open(dict_path) as f:
        return json.load(f)
