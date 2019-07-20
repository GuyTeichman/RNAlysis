import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from rnalysis.general import *
from rnalysis.feature_name_to_wbgene import *


def test_wbgene_translator_api():
    fname = "all_feature_96_new.csv"
    assert(GeneNameTranslator(fname))

