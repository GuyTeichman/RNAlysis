"""
This module can perform enrichment analyses on a given set of genomic features and visualize their intersections. \
These include gene ontology/tissue/phenotype enrichment, enrichment for user-defined attributes, \
set visualization ,etc. \
Results of enrichment analyses can be saved to .csv files.
"""
import itertools
import types
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib_venn as vn
import numpy as np
import pandas as pd
import upsetplot
from matplotlib.cm import ScalarMappable

from rnalysis.filtering import Filter
from rnalysis.utils import io, parsing, ref_tables, validation, enrichment_runner


class FeatureSet:
    """ Receives a filtered gene set and the set's name (optional) and preforms various enrichment analyses on them. """
    __slots__ = {'gene_set': 'set of feature names/indices', 'set_name': 'name of the FeatureSet'}

    def __init__(self, gene_set: Union[List[str], Set[str], 'Filter'] = None, set_name: str = ''):

        """
        :param gene_set: the set of genomic features to be used in downstream analyses
        :type gene_set: filtering.Filter object, set of strings or list of strings
        :param set_name: name of the FeatureSet
        :type set_name: str


        :Examples:
            >>> from rnalysis import enrichment, filtering
            >>> my_set = enrichment.FeatureSet({'gene1','gene2','gene2'}, 'name of my set')

            >>> filter_obj = filtering.CountFilter('tests/test_files/counted.csv')
            >>> my_other_set = enrichment.FeatureSet(filter_obj, 'name of my other set')

        """
        assert isinstance(set_name, str), f"'set_name' must be of type str, instead got {type(set_name)}."
        if gene_set is None:
            gene_set = parsing.data_to_set(parsing.from_string(
                "Please insert genomic features/indices separated by newline \n"
                "(example: \n'WBGene00000001\nWBGene00000002\nWBGene00000003')", delimiter='\n'))
        elif validation.isinstanceinh(gene_set, Filter):
            gene_set = gene_set.index_set
        else:
            gene_set = parsing.data_to_set(gene_set)
        self.gene_set = gene_set
        self.set_name = set_name

    def __repr__(self):
        return f"{self.__class__.__name__}: '{self.set_name}'"

    def __len__(self):
        return len(self.gene_set)

    def __contains__(self, item):
        return True if item in self.gene_set else False

    def change_set_name(self, new_name: str):
        """
        Change the 'set_name' of a FeatureSet to a new name.

        :param new_name: the new set name
        :type new_name: str

        """
        assert isinstance(new_name, str), f"New set name must be of type str. Instead, got {type(new_name)}"
        self.set_name = new_name

    def save_txt(self, fname: Union[str, Path]):

        """
        Save the list of features in the FeatureSet object under the specified filename and path.

        :type fname: str or pathlib.Path
        :param fname: full filename/path for the output file. Can include the '.txt' suffix but doesn't have to.

        """
        assert isinstance(fname, (str, Path)), "fname must be str or pathlib.Path!"
        if isinstance(fname, str):
            if not fname.endswith('.txt'):
                fname = fname + '.txt'
        elif isinstance(fname, Path):
            if not fname.suffix == '.txt':
                fname = Path(f"{str(fname.parent)}{fname.name}.txt")
        with open(fname, 'w') as f:
            for gene in self.gene_set:
                f.write(gene + '\n')

    def _set_ops(self, others: Union[set, 'FeatureSet', Tuple[Union[set, 'FeatureSet']]],
                 op: types.FunctionType) -> set:

        """
        Performs a given set operation on self and on another object (FeatureSet or set).
        :type others: FeatureSet or set
        :param others: Other object to perform set operation with.
        :type: op: Callable (set.union, set.intersection, set.difference or set.symmetric difference)
        :param op: The set operation to be performed.
        :return: A set resulting from the set operation.

        """
        others = list(others)
        for i, other in enumerate(others):
            if isinstance(other, set):
                pass
            elif isinstance(other, FeatureSet):
                others[i] = other.gene_set
            else:
                raise TypeError("'other' must be an FeatureSet object or a set!")
        try:
            return op(self.gene_set, *others)
        except TypeError as e:
            if op == set.symmetric_difference:
                raise TypeError(
                    f"Symmetric difference can only be calculated for two objects, {len(others) + 1} were given!")
            else:
                raise e

    def union(self, *others: Union[set, 'FeatureSet']) -> 'FeatureSet':

        """
         Calculates the set union of the indices from multiple FeatureSet objects \
        (the indices that exist in at least one of the FeatureSet objects).

        :type others: FeatureSet or set
        :param others: The objects against which the current object will be compared.
        :return: a new FeatureSet with elements from this FeatureSet and all other objects.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000004','WBGene00000005','WBGene00000006'}, 'set name')
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> en.union(s, en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000003', 'WBGene00000004', 'WBGene00000001', 'WBGene00000002', 'WBGene00000006', 'WBGene00000005'}

        """
        return FeatureSet(self._set_ops(others, set.union))

    def intersection(self, *others: Union[set, 'FeatureSet']) -> 'FeatureSet':

        """
        Calculates the set intersection of the indices from multiple FeatureSet objects \
        (the indices that exist in ALL of the FeatureSet objects).

        :type others: FeatureSet, RankedSet or set
        :param others: The objects against which the current object will be compared.
        :return: a new FeatureSet with elements common to this FeatureSet and all other objects.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000001','WBGene00000002','WBGene00000006'}, 'set name')
            >>> s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> en.intersection(s, en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000001'}

        """
        return FeatureSet(self._set_ops(others, set.intersection))

    def difference(self, *others: Union[set, 'FeatureSet']) -> 'FeatureSet':

        """
        Calculates the set difference of the indices from multiple FeatureSet objects \
        (the indices that appear in the first FeatureSet object but NOT in the other objects).

        :type others: FeatureSet, RankedSet or set
        :param others: The objects against which the current object will be compared.
        :return: a new FeatureSet with elements in this FeatureSet that are not in the other objects.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000001','WBGene00000002','WBGene00000006'}, 'set name')
            >>> s = {'WBGene00000001','WBGene00000002','WBGene00000003'}
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> en.difference(s, en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000006'}

        """
        return FeatureSet(self._set_ops(others, set.difference))

    def symmetric_difference(self, other: Union[set, 'FeatureSet']) -> 'FeatureSet':

        """
        Calculates the set symmetric difference of the indices from two FeatureSet objects \
        (the indices that appear in EXACTLY ONE of the FeatureSet objects, and not both/neither). \
        A-symmetric difference-B is equivalent to (A-difference-B)-union-(B-difference-A).

        :type other: FeatureSet, RankedSet or set
        :param other: A second object against which the current object will be compared.
        :return: a new FeatureSet with elements in either this FeatureSet or the other object, but not both.
        :rtype: FeatureSet

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000001','WBGene00000002','WBGene00000006'}, 'set name')
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> en.symmetric_difference(en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000002', 'WBGene00000006', 'WBGene00000004'}

        """
        return FeatureSet(self._set_ops((other,), set.symmetric_difference))

    @staticmethod
    def _enrichment_save_csv(df: pd.DataFrame, fname: str):

        """
        Internal method, used to save enrichment results to .csv files. Static class method.

        :param df: pandas DataFrame to be saved.
        :param fname: Name and full path under which the DataFrame will be saved

        """
        if fname is None:
            fname = input("Please insert the full name and path to save the file to")
        else:
            assert isinstance(fname, (str, Path))
        if isinstance(fname, Path):
            fname = str(Path)
        io.save_csv(df, filename=fname if fname.endswith('.csv') else fname + '.csv')

    def go_enrichment(self, organism: Union[str, int] = 'auto', gene_id_type: str = 'UniProtKB', alpha: float = 0.05,
                      statistical_test: str = 'fisher', biotype: str = 'protein_coding',
                      background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                      biotype_ref_path: str = 'predefined', propagate_annotations: str = 'elim',
                      aspects: Union[str, Iterable[str]] = 'any', evidence_types: Union[str, Iterable[str]] = 'any',
                      excluded_evidence_types: Union[str, Iterable[str]] = (),
                      databases: Union[str, Iterable[str]] = 'any',
                      excluded_databases: Union[str, Iterable[str]] = (),
                      qualifiers: Union[str, Iterable[str]] = 'any',
                      excluded_qualifiers: Union[str, Iterable[str]] = 'not', return_nonsignificant: bool = False,
                      save_csv: bool = False, fname=None, return_fig: bool = False, plot_horizontal: bool = True,
                      randomization_reps: int = 10000, random_seed: int = None,
                      parallel: bool = True) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        """
        Calculates enrichment and depletion of the FeatureSet for Gene Ontology (GO) terms against a background set \
        using the Hypergeometric Test. The GO terms and annotations are drawn via the GO Solr search engine GOlr, \
        using the search terms defined by the user. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. \
        P-values are calculated using a hypergeometric test: \
        Given M genes in the background set, n genes in the test set, \
        with N genes from the background set belonging to a specific attribute (or 'success') \
        and X genes from the test set belonging to that attribute. \
        If we were to randomly draw n genes from the background set (without replacement), \
        what is the probability of drawing X or more (in case of enrichment)/X or less (in case of depletion) \
        genes belonging to the given attribute? \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment score = 0) \
        appears with the smallest value in the scale.

        :param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the GOLR server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str (default='UniProtKB')
        :type alpha: float between 0 and 1
        :param alpha: Indicates the FDR threshold for significance.
        :param statistical_test: determines the statistical test to be used for enrichment analysis. \
        Note that some propagation methods support only some of the available statistical tests.
        :type statistical_test: 'fisher', 'hypergeometric' or 'randomization' (default 'fisher')
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :param propagate_annotations: determines the propagation method of GO annotations. \
        'no' does not propagate annotations at all; 'classic' propagates all annotations up to the DAG tree's root; \
        'elim' terminates propagation at nodes which show significant enrichment; 'weight' performs propagation in a \
        weighted manner based on the significance of children nodes relatively to their parents; and 'allm' uses a \
        combination of all proopagation methods. To read more about the propagation methods, \
        see Alexa et al: https://pubmed.ncbi.nlm.nih.gov/16606683/
        :type propagate_annotations: 'classic', 'elim', 'weight', 'all.m', or 'no' (default 'elim')
        :param aspects: only annotations from the specified GO aspects will be included in the analysis. \
        Legal aspects are 'biological_process' (P), 'molecular_function' (F), and 'cellular_component' (C).
        :type aspects: str, Iterable of str, 'biological_process', 'molecular_function', 'cellular_component', \
        or 'any' (default='any')
        :param evidence_types: only annotations with the specified evidence types will be included in the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', 'author', \
        'curator', 'electronic', or 'any' (default='any')
        :param excluded_evidence_types: annotations with the specified evidence types will be \
        excluded from the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type excluded_evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', \
        'author', 'curator', 'electronic', or None (default=None)
        :param databases: only annotations from the specified databases will be included in the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type databases: str, Iterable of str, or 'any' (default)
        :param excluded_databases: annotations from the specified databases will be excluded from the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type excluded_databases: str, Iterable of str, or None (default)
        :param qualifiers: only annotations with the speficied qualifiers will be included in the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type qualifiers: str, Iterable of str, or 'any' (default)
        :param excluded_qualifiers: annotations with the speficied qualifiers will be excluded from the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type excluded_qualifiers: str, Iterable of str, or None (default 'not')
        :param return_nonsignificant: if True, the results DataFrame will include all tested GO terms - \
        both significant and non-significant terms. If False (default), only significant GO terms will be returned.
        :type return_nonsignificant: bool (default False)
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. \
        Otherwise, results will be plotted with a vertical plot.
        :type random_seed: non-negative integer (default None)
        :type random_seed: if using a randomization test, determine the random seed used to initialize \
        the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs. If not using a randomization test, this parameter will not affect the analysis.
        :param randomization_reps: if using a randomization test, determine how many randomization repititions to run. \
        Otherwise, this parameter will not affect the analysis.
        :type randomization_reps: int larger than 0 (default 10000)
        :type parallel: bool (default False)
        :param parallel: if True, will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results_go.png
           :align:   center
           :scale: 60 %

           Example plot of go_enrichment()


        .. figure::  plot_enrichment_results_go_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of go_enrichment(plot_horizontal = False)
        """
        propagate_annotations = propagate_annotations.lower()
        if validation.isinstanceinh(background_genes, FeatureSet):
            background_genes = background_genes.gene_set
        plot_dag = False
        if statistical_test.lower() == 'randomization':
            kwargs = dict(reps=randomization_reps, random_seed=random_seed)
        else:
            kwargs = {}
        runner = enrichment_runner.GOEnrichmentRunner(self.gene_set, organism, gene_id_type, alpha,
                                                      propagate_annotations, aspects, evidence_types,
                                                      excluded_evidence_types, databases, excluded_databases,
                                                      qualifiers, excluded_qualifiers, return_nonsignificant, save_csv,
                                                      fname, return_fig, plot_horizontal, plot_dag, self.set_name,
                                                      parallel, statistical_test, biotype, background_genes,
                                                      biotype_ref_path, **kwargs)

        return runner.run()

    def enrich_randomization_parallel(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                                      fdr: float = 0.05, reps: int = 10000,
                                      biotype: str = 'protein_coding',
                                      background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                                      attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                                      save_csv: bool = False, fname=None, return_fig: bool = False,
                                      plot_horizontal: bool = True, random_seed: int = None
                                      ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:

        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set \
        using a randomization test running in parallel. \
        This function is depracated, and is replaced by enrich_randomization(parallel_processing=True). \
        It will be removed in future versions of RNAlysis. \
        The attributes are drawn from an Attribute Reference Table. \
        Parallel processing makes this function generally faster than FeatureSet.enrich_randomization. \
        Results should otherwise be the same between the two functions. \
        To use it you must first start a parallel session, using the function 'general.start_parallel_session()'. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. P-values are calculated using a \
        randomization test with the formula p = (successes + 1)/(repeats + 1). \
        This formula results in a positively-biased estimator of the real p-value \
        (a conservative estimate of p-value). When the number of reps approaches infinity, \
        the formula results in an unbiased estimator of the real p-value. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
        :type reps: int larger than 0
        :param reps: How many repetitions to run the randomization for. \
        10,000 is the default. Recommended 10,000 or higher.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        Cannot be specified together with 'biotype'.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :type random_seed: non-negative integer (default None)
        :type random_seed: The random seed used to initialize the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization_parallel()


        .. figure::  plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization_parallel(plot_horizontal = False)

       """
        warnings.warn("FeatureSet.enrich_randomization_parallel() is deprecated and will be removed "
                      "in a future release. Please use Featureset.enrich_randomization() "
                      "with the parameter 'parallel_processing=True' instead.", DeprecationWarning)
        return self.enrich_randomization(attributes, fdr, reps, biotype, background_genes, attr_ref_path,
                                         biotype_ref_path, save_csv, fname, return_fig, plot_horizontal, random_seed,
                                         parallel=True)

    def enrich_randomization(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                             alpha: float = 0.05, reps: int = 10000, biotype: str = 'protein_coding',
                             background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                             attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                             save_csv: bool = False, fname=None, return_fig: bool = False, plot_horizontal: bool = True,
                             random_seed: int = None, parallel: bool = False
                             ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:

        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set \
        using a randomization test. The attributes are drawn from an Attribute Reference Table. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. P-values are calculated using a \
        randomization test with the formula p = (successes + 1)/(repeats + 1). \
        This formula results in a positively-biased estimator of the real p-value \
        (a conservative estimate of p-value). When the number of reps approaches infinity, \
        the formula results in an unbiased estimator of the real p-value. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type alpha: float between 0 and 1 (default 0.05)
        :param alpha: Indicates the FDR threshold for significance.
        :type reps: int larger than 0 (default 10000)
        :param reps: How many repetitions to run the randomization for. \
        10,000 is the default. Recommended 10,000 or higher.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        (default 'protein_coding')
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object \
        (default None)
        :param background_genes: a set of specific feature indices to be used as background genes. \
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        Cannot be specified together with 'biotype'.
        :type save_csv: bool (default False)
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path (default None)
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :type random_seed: non-negative integer (default None)
        :type random_seed: The random seed used to initialize the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs.
        :type parallel: bool (default False)
        :param parallel: if True, will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization()


        .. figure::  plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization(plot_horizontal = False)

        """
        if validation.isinstanceinh(background_genes, FeatureSet):
            background_genes = background_genes.gene_set
        runner = enrichment_runner.EnrichmentRunner(self.gene_set, attributes, alpha, attr_ref_path, save_csv, fname,
                                                    return_fig, plot_horizontal, self.set_name, parallel,
                                                    'randomization', biotype, background_genes, biotype_ref_path,
                                                    single_list=False, random_seed=random_seed, reps=reps)
        return runner.run()

    def enrich_hypergeometric(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                              alpha: float = 0.05, biotype: str = 'protein_coding',
                              background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                              attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                              save_csv: bool = False, fname=None, return_fig: bool = False,
                              plot_horizontal: bool = True, parallel: bool = True
                              ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:

        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set \
        using the Hypergeometric Test. The attributes are drawn from an Attribute Reference Table. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. \
        P-values are calculated using a hypergeometric test: \
        Given M genes in the background set, n genes in the test set, \
        with N genes from the background set belonging to a specific attribute (or 'success') \
        and X genes from the test set belonging to that attribute. \
        If we were to randomly draw n genes from the background set (without replacement), \
        what is the probability of drawing X or more (in case of enrichment)/X or less (in case of depletion) \
        genes belonging to the given attribute? \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment score = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type alpha: float between 0 and 1 (default 0.05)
        :param alpha: Indicates the FDR threshold for significance.
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all' \
        (default 'protein_coding')
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object \
        (default None)
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :type save_csv: bool (default False)
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path (default None)
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :type parallel: bool (default False)
        :param parallel: if True, will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_hypergeometric()


        .. figure::  plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_hypergeometric(plot_horizontal = False)

        """
        if validation.isinstanceinh(background_genes, FeatureSet):
            background_genes = background_genes.gene_set
        runner = enrichment_runner.EnrichmentRunner(self.gene_set, attributes, alpha, attr_ref_path, save_csv, fname,
                                                    return_fig, plot_horizontal, self.set_name, parallel,
                                                    'hypergeometric', biotype, background_genes, biotype_ref_path,
                                                    single_list=False)
        return runner.run()

    @staticmethod
    def plot_enrichment_results(df: pd.DataFrame, fdr=0.05, en_score_col: str = 'log2_fold_enrichment',
                                name_col: str = None, title: str = '', center_bars: bool = True,
                                plot_horizontal: bool = True, ylabel: str = r"$\log_2$(Fold Enrichment)") -> plt.Figure:

        """
        Receives a DataFrame output from an enrichment function and plots it in a bar plot. \
        For the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :param fdr:
        :type fdr:
        :param name_col:
        :type name_col:
        :param df: a pandas DataFrame created by FeatureSet.enrich_randomization.
        :type df: pd.DataFrame
        :param en_score_col: name of the DataFrame column that contains the enrichment scores.
        :type en_score_col: str (default 'log2_fold_enrichment')
        :param title: plot title.
        :type title: str
        :param ylabel: plot ylabel.
        :type ylabel: str
        :param plot_horizontal:
        :type plot_horizontal: bool (default True)
        :param center_bars: if True, centers the bars around Y=0. Otherwise, ylim is determined by min/max values.
        :type center_bars: bool (default True)
        :return: Figure object containing the bar plot
        :rtype: matplotlib.figure.Figure instance
        """
        plt.style.use('seaborn-white')
        # choose functions and parameters according to the graph's orientation (horizontal vs vertical)
        if plot_horizontal:
            figsize = [14, 0.4 * (6.4 + df.shape[0])]
            bar_func = plt.Axes.barh
            line_func = plt.Axes.axvline
            cbar_kwargs = dict(location='bottom')
            tick_func = plt.Axes.set_yticks
            ticklabels_func = plt.Axes.set_yticklabels
            ticklabels_kwargs = dict(fontsize=13, rotation=0)
            df = df[::-1]
        else:
            figsize = [0.5 * (6.4 + df.shape[0]), 5.6]
            bar_func = plt.Axes.bar
            line_func = plt.Axes.axhline
            cbar_kwargs = dict(location='left')
            tick_func = plt.Axes.set_xticks
            ticklabels_func = plt.Axes.set_xticklabels
            ticklabels_kwargs = dict(fontsize=13, rotation=45)

        # pull names/scores/pvals out to avoid accidentally changing the results DataFrame in-place
        enrichment_names = df.index.values.tolist() if name_col is None else df[name_col].values.tolist()
        enrichment_scores = df[en_score_col].values.tolist()
        enrichment_pvalue = df['padj'].values.tolist()

        # set enrichment scores which are 'inf' or '-inf' to be the second highest/lowest enrichment score in the list
        scores_no_inf = [i for i in enrichment_scores if i != np.inf and i != -np.inf and i < 0]
        if len(scores_no_inf) == 0:
            scores_no_inf.append(-1)
        for i in range(len(enrichment_scores)):
            if enrichment_scores[i] == -np.inf:
                enrichment_scores[i] = min(scores_no_inf)
        max_score = max(np.max(np.abs(enrichment_scores)), 2)

        # get color values for bars
        data_color_norm = [0.5 * (1 + i / (np.floor(max_score) + 1)) * 255 for i in enrichment_scores]
        data_color_norm_8bit = [int(i) if i != np.inf and i != -np.inf else np.sign(i) * max(np.abs(scores_no_inf)) for
                                i in data_color_norm]
        my_cmap = plt.cm.get_cmap('coolwarm')
        colors = my_cmap(data_color_norm_8bit)

        # generate bar plot
        fig, ax = plt.subplots(constrained_layout=True, figsize=figsize)
        bar = bar_func(ax, range(len(enrichment_names)), enrichment_scores, color=colors, edgecolor='black',
                       linewidth=1, zorder=2)
        bar.tick_labels = enrichment_names
        # determine bound, and enlarge the bound by a small margin (0.2%) so nothing gets cut out of the figure
        bounds = np.array([np.ceil(-max_score) - 1, (np.floor(max_score) + 1)]) * 1.002
        # add black line at y=0 and grey lines at every round positive/negative integer in range
        for ind in range(int(bounds[0]) + 1, int(bounds[1]) + 1):
            color = 'black' if ind == 0 else 'grey'
            linewidth = 1 if ind == 0 else 0.5
            linestyle = '-' if ind == 0 else '-.'
            line_func(ax, ind, color=color, linewidth=linewidth, linestyle=linestyle, zorder=0)
        # add colorbar
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(*bounds))
        sm.set_array(np.array([]))
        cbar_label_kwargs = dict(label=ylabel, fontsize=16, labelpad=15)
        cbar = fig.colorbar(sm, ticks=range(int(bounds[0]), int(bounds[1]) + 1), **cbar_kwargs)
        cbar.set_label(**cbar_label_kwargs)
        cbar.ax.tick_params(labelsize=14, pad=6)
        # apply xticks
        tick_func(ax, range(len(enrichment_names)))
        ticklabels_func(ax, enrichment_names, **ticklabels_kwargs)
        # title
        ax.set_title(title, fontsize=18)
        # add significance asterisks
        for col, sig in zip(bar, enrichment_pvalue):
            asterisks, fontweight = enrichment_runner.EnrichmentRunner._get_pval_asterisk(sig, fdr)
            if plot_horizontal:
                x = col._width
                y = col.xy[1] + 0.5 * col._height
                valign = 'center'
                halign = 'left' if np.sign(col._width) == 1 else 'right'
                rotation = 270 if np.sign(col._width) == 1 else 90
            else:
                x = col.xy[0] + 0.5 * col._width
                y = col._height
                valign = 'bottom' if np.sign(col._height) == 1 else 'top'
                halign = 'center'
                rotation = 0

            ax.text(x=x, y=y, s=asterisks, fontname='DejaVu Sans', fontweight=fontweight, rotation=rotation,
                    fontsize=9, horizontalalignment=halign, verticalalignment=valign, zorder=1)
        # despine
        _ = [ax.spines[side].set_visible(False) for side in ['top', 'right']]
        # center bars
        if center_bars:
            if plot_horizontal:
                ax.set_xbound(bounds)
                plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            else:
                ax.set_ybound(bounds)
                plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        plt.show()
        return fig

    def enrich_non_categorical(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                               alpha: float = 0.05, parametric_test: bool = False, biotype: str = 'protein_coding',
                               background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                               attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                               plot_log_scale: bool = True, plot_style: str = 'overlap', n_bins: int = 50,
                               save_csv: bool = False, fname=None, return_fig: bool = False
                               ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, List[plt.Figure]]]:

        if validation.isinstanceinh(background_genes, FeatureSet):
            background_genes = background_genes.gene_set
        runner = enrichment_runner.NonCategoricalEnrichmentRunner(self.gene_set, attributes, alpha, biotype,
                                                                  background_genes, attr_ref_path,
                                                                  biotype_ref_path, save_csv, fname,
                                                                  return_fig, plot_log_scale, plot_style,
                                                                  n_bins, self.set_name, parallel=False,
                                                                  parametric_test=parametric_test)
        return runner.run()

    @staticmethod
    def plot_enrichment_hist(attribute: str, gene_set: set, set_name: str, attr_srs: pd.Series,
                             plot_log_scale: bool, plot_style: str, n_bins: int, parametric_test: bool,
                             pval: float) -> plt.Figure:
        assert isinstance(n_bins,
                          int) and n_bins > 0, f"'n_bins' must be a positive integer. Instead got {type(n_bins)}."
        # generate observed and expected Series, either linear or in log10 scale
        exp = attr_srs
        obs = exp.loc[gene_set]
        if plot_log_scale:
            xlabel = r"$\log_{10}$" + f"({attribute})"
            exp = np.log10(exp)
            obs = np.log10(obs)
        else:
            xlabel = f"{attribute}"

        # determine bins according to value range and 'n_bins'
        bins = np.linspace(np.min(exp), np.max(exp), n_bins).squeeze()

        # generate histogram according to plot style
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        kwargs = dict(bins=bins, density=True, alpha=0.5, edgecolor='black', linewidth=1)
        colors = ['C0', 'C1']
        if plot_style.lower() == 'interleaved':
            y, x, _ = ax.hist([exp.values, obs.values], **kwargs, color=colors, label=['Expected', 'Observed'])
            max_y_val = np.max(y)
        elif plot_style.lower() == 'overlap':
            y, _, _ = ax.hist(exp.values, **kwargs, color=colors[0], label='Expected')
            y2, _, _ = ax.hist(obs.values, **kwargs, color=colors[1], label='Observed')
            max_y_val = np.max([np.max(y), np.max(y2)])
        else:
            raise ValueError(f"Invalid value for 'plot_style': '{plot_style}'")

        # set either mean or median as the measure of centrality
        if parametric_test:
            x_exp, x_obs = exp.mean(), obs.mean()
            label_exp, label_obs = 'Expected mean', 'Observed mean'
        else:
            x_exp, x_obs = exp.median(), obs.median()
            label_exp, label_obs = 'Expected median', 'Observed median'

        # add lines for mean/median of observed and expected distributions
        ax.vlines(x_exp, ymin=0, ymax=max_y_val * 1.1, color='blue', linestyle='dashed', linewidth=2,
                  label=label_exp)
        ax.vlines(x_obs, ymin=0, ymax=max_y_val * 1.1, color='red', linestyle='dashed', linewidth=2,
                  label=label_obs)

        # add significance notation
        asterisks, fontweight = FeatureSet._get_pval_asterisk(pval)
        ax.vlines([x_exp, x_obs], ymin=max_y_val * 1.12, ymax=max_y_val * 1.16, color='k', linewidth=1)
        ax.hlines(max_y_val * 1.16, xmin=min(x_exp, x_obs), xmax=max(x_exp, x_obs), color='k', linewidth=1)
        ax.text(np.mean([x_exp, x_obs]), max_y_val * 1.17, asterisks, horizontalalignment='center',
                fontweight=fontweight)

        # legend and titles
        ax.legend()
        obs_name = 'Observed' if set_name == '' else set_name
        ax.set_title(f"Histogram of {attribute} - {obs_name} vs Expected", fontsize=17)
        ax.set_ylabel("Probability density", fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        plt.show()

        return fig

    def biotypes(self, ref: str = 'predefined'):

        """
        Returns a DataFrame of the biotypes in the gene set and their count.

        :type ref: str or pathlib.Path (default 'predefined')
        :param ref: Path of the reference file used to determine biotype. \
        Default is the path predefined in the settings file.

        :Examples:
            >>> from rnalysis import enrichment, filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> en = enrichment.FeatureSet(d)
            >>> en.biotypes(ref='tests/biotype_ref_table_for_tests.csv')
                            gene
            biotype
            protein_coding    26
            pseudogene         1
            unknown            1

        """

        ref = ref_tables.get_biotype_ref_path(ref)
        ref_df = io.load_csv(ref)
        validation.validate_biotype_table(ref_df)
        ref_df.columns = ref_df.columns.str.lower()
        not_in_ref = pd.Index(self.gene_set).difference(set(ref_df['gene']))
        if len(not_in_ref) > 0:
            warnings.warn(
                f'{len(not_in_ref)} of the features in the Filter object do not appear in the Biotype Reference Table. ')
            ref_df = ref_df.append(pd.DataFrame({'gene': not_in_ref, 'biotype': 'not_in_biotype_reference'}))
        return ref_df.set_index('gene', drop=False).loc[self.gene_set].groupby('biotype').count()


class RankedSet(FeatureSet):
    """
    Receives a ranked gene set, sorted by any biologically-meaningful value (expression level, fold-change, etc)\
     and preforms various enrichment analyses on them. \
     ALl functions that can be applied to FeatureSet objects can also be applied to RankedSet objects.
    """
    __slots__ = {'ranked_genes': 'a vector of feature names/indices ordered by rank'}

    def __init__(self, ranked_genes: Union[Filter, List[str], Tuple[str], np.ndarray], set_name: str = ''):

        if validation.isinstanceinh(ranked_genes, Filter):
            self.ranked_genes = ranked_genes.df.index.values.astype('str', copy=True)
        elif isinstance(ranked_genes, (list, tuple)):
            self.ranked_genes = np.array(ranked_genes, dtype='str')
        elif isinstance(ranked_genes, np.ndarray):
            self.ranked_genes = ranked_genes.astype('str', copy=True)
        elif isinstance(ranked_genes, set):
            raise TypeError("'ranked_genes' must be an array, list, tuple or Filter object, sorted by rank. "
                            "Python sets are not a valid type ofr 'ranked_genes'.")
        elif isinstance(ranked_genes, dict):
            raise TypeError("Generating a RankedSet from a dictionary is not implemented yet.")
        else:
            raise TypeError(f"'ranked_genes' must be an array, list, tuple or Filter object, sorted by rank. "
                            f"Instead got {type(ranked_genes)}.")

        super().__init__(ranked_genes, set_name)
        assert len(self.ranked_genes) == len(self.gene_set), f"'ranked_genes' must have no repeating elements!"

    def _set_ops(self, others: Union[set, 'FeatureSet'], op: types.FunctionType):
        warnings.warn("Warning: when performing set operations with RankedSet objects, "
                      "the return type will always be FeatureSet and not RankedSet.")
        return super()._set_ops(others, op)

    def go_enrichment_single_list(self, organism: Union[str, int] = 'auto', gene_id_type: str = 'UniProtKB',
                                  alpha: float = 0.05, propagate_annotations: str = 'elim',
                                  aspects: Union[str, Iterable[str]] = 'any',
                                  evidence_types: Union[str, Iterable[str]] = 'any',
                                  excluded_evidence_types: Union[str, Iterable[str]] = (),
                                  databases: Union[str, Iterable[str]] = 'any',
                                  excluded_databases: Union[str, Iterable[str]] = (),
                                  qualifiers: Union[str, Iterable[str]] = 'any',
                                  excluded_qualifiers: Union[str, Iterable[str]] = 'not',
                                  return_nonsignificant: bool = False,
                                  save_csv: bool = False, fname=None,
                                  return_fig: bool = False, plot_horizontal: bool = True, parallel: bool = False
                                  ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        """
        Calculates enrichment and depletion of the sorted RankedSet for Gene Ontology (GO) terms \
        WITHOUT a background set, using the generalized Minimum Hypergeometric Test (XL-mHG, developed by  \
        `Prof. Zohar Yakhini and colleagues <https://dx.doi.org/10.1371/journal.pcbi.0030039/>`_ \
        and generalized by \
        `Florian Wagner <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143196/>`_). \
        The GO terms and annotations are drawn via the GO Solr search engine GOlr, \
        using the search terms defined by the user. \
        P-values are calculated using using the generalized Minimum Hypergeometric Test. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the GOLR server do not match your gene_id_type, \
        RNAlysis will attempt to map  the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str (default='UniProtKB')
        :type alpha: float between 0 and 1
        :param alpha: Indicates the FDR threshold for significance.
        :param propagate_annotations: determines the propagation method of GO annotations. \
        'no' does not propagate annotations at all; 'classic' propagates all annotations up to the DAG tree's root; \
        'elim' terminates propagation at nodes which show significant enrichment; 'weight' performs propagation in a \
        weighted manner based on the significance of children nodes relatively to their parents; and 'allm' uses a \
        combination of all proopagation methods. To read more about the propagation methods, \
        see Alexa et al: https://pubmed.ncbi.nlm.nih.gov/16606683/
        :type propagate_annotations: 'classic', 'elim', 'weight', 'all.m', or 'no' (default 'elim')
        :param aspects: only annotations from the specified GO aspects will be included in the analysis. \
        Legal aspects are 'biological_process' (P), 'molecular_function' (F), and 'cellular_component' (C).
        :type aspects: str, Iterable of str, 'biological_process', 'molecular_function', 'cellular_component', \
        or 'any' (default='any')
        :param evidence_types: only annotations with the specified evidence types \
        will be included in the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', 'author', \
        'curator', 'electronic', or 'any' (default='any')
        :param excluded_evidence_types: annotations with the specified evidence types will be \
        excluded from the analysis. \
        For a full list of legal evidence codes and evidence code categories see the GO Consortium website: \
        http://geneontology.org/docs/guide-go-evidence-codes/
        :type excluded_evidence_types: str, Iterable of str, 'experimental', 'phylogenetic' ,'computational', \
        'author', 'curator', 'electronic', or None (default=None)
        :param databases: only annotations from the specified databases will be included in the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type databases: str, Iterable of str, or 'any' (default)
        :param excluded_databases: annotations from the specified databases \
        will be excluded from the analysis. \
        For a full list of legal databases see the GO Consortium website:
        http://amigo.geneontology.org/xrefs
        :type excluded_databases: str, Iterable of str, or None (default)
        :param qualifiers: only annotations with the speficied qualifiers will be included in the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type qualifiers: str, Iterable of str, or 'any' (default)
        :param excluded_qualifiers: annotations with the speficied qualifiers \
        will be excluded from the analysis. \
        Legal qualifiers are 'not', 'contributes_to', and/or 'colocalizes_with'.
        :type excluded_qualifiers: str, Iterable of str, or None (default 'not')
        :param return_nonsignificant: if True, the results DataFrame will include all tested GO terms - \
        both significant and non-significant terms. If False (default), \
        only significant GO terms will be returned.
        :type return_nonsignificant: bool (default False)
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), \
        fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. \
        Otherwise, results will be plotted with a vertical plot.
        :type parallel: bool (default False)
        :param parallel: if True, will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results_go_singlelist.png
           :align:   center
           :scale: 60 %

           Example plot of go_enrichment_single_list()


        .. figure::  plot_enrichment_results_go_singlelist_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of go_enrichment_single_list(plot_horizontal = False)
        """
        plot_dag = False
        runner = enrichment_runner.GOEnrichmentRunner(self.ranked_genes, organism, gene_id_type, alpha,
                                                      propagate_annotations, aspects, evidence_types,
                                                      excluded_evidence_types, databases, excluded_databases,
                                                      qualifiers, excluded_qualifiers, return_nonsignificant, save_csv,
                                                      fname, return_fig, plot_horizontal, plot_dag, self.set_name,
                                                      parallel=parallel, enrichment_func_name='xlmhg', single_list=True)

        return runner.run()

    def enrich_single_list(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, alpha: float = 0.05,
                           attr_ref_path: str = 'predefined', save_csv: bool = False, fname=None,
                           return_fig: bool = False, plot_horizontal: bool = True, parallel: bool = False):
        """
        Calculates enrichment and depletion of the sorted RankedSet for user-defined attributes \
        WITHOUT a background set, using the generalized Minimum Hypergeometric Test (XL-mHG, developed by  \
        `Prof. Zohar Yakhini and colleagues <https://dx.doi.org/10.1371/journal.pcbi.0030039/>`_ \
        and generalized by \
        `Florian Wagner <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143196/>`_). \
        The attributes are drawn from an Attribute Reference Table. \
        P-values are calculated using using the generalized Minimum Hypergeometric Test. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type alpha: float between 0 and 1
        :param alpha: Indicates the FDR threshold for significance.
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :type parallel: bool (default False)
        :param parallel: if True, will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  plot_enrichment_results_single_list.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_single_list()


        .. figure::  plot_enrichment_results_single_list_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_single_list(plot_horizontal = False)

        """
        runner = enrichment_runner.EnrichmentRunner(self.ranked_genes, attributes, alpha, attr_ref_path, save_csv,
                                                    fname, return_fig, plot_horizontal, self.set_name,
                                                    parallel=parallel, enrichment_func_name='xlmhg', single_list=True)
        return runner.run()


def _fetch_sets(objs: dict, ref: str = 'predefined'):
    """
    Receives the 'objs' input from enrichment.upset_plot() and enrichment.venn_diagram(), and turns the values in it \
    into python sets.

    :param objs: the 'objs' input given to the function enrichment.upset_plot() or enrichment.venn_diagram().
    :type objs: a dictionary, where the keys are names of sets, and the values are either\
     python sets, FeatureSets or names of columns in the Attribute Reference Table.
    :param ref: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type ref: str or pathlib.Path (default 'predefined')
    :return: a dictionary, where the keys are names of sets and the values are python sets of feature indices.
    """
    assert isinstance(objs, dict), f"objs must be a dictionary. Instaed got {type(objs)}"
    for obj in objs:
        if isinstance(objs[obj], set):
            pass
        elif validation.isinstanceinh(objs[obj], Filter):
            objs[obj] = objs[obj].index_set
        elif isinstance(objs[obj], FeatureSet):
            objs[obj] = objs[obj].gene_set
        elif isinstance(objs[obj], str):
            if 'attr_table' not in locals():
                pth = ref_tables.get_attr_ref_path(ref)
                attr_table = io.load_csv(pth, 0)
            attr = objs[obj]
            myset = set(attr_table[attr].loc[attr_table[attr].notna()].index)
            objs[obj] = myset
        else:
            raise TypeError
    return objs


def upset_plot(objs: Dict[str, Union[str, FeatureSet, Set[str]]], title: str = '', ref: str = 'predefined'):
    """
    Generate an UpSet plot of 2 or more sets, FeatureSets or attributes from the Attribute Reference Table.


    :param objs: the FeatureSets, python sets or user-defined attributes to plot.
    :type objs: a dictionary with 2 or more entries, where the keys are the names of the sets, and the values are either \
    a FeatureSet, a python set of feature indices, or a name of a column in the Attribute Reference Table. \
    For example: \
    {'first set':{'gene1','gene2','gene3'}, 'second set':'name_of_attribute_from_reference_table'}
    :param title: determines the title of the plot.
    :type title: str
    :param ref: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type ref: str or pathlib.Path (default 'predefined')
    :returns: a dictionary of matplotlib axes, where the keys are 'matrix', 'intersections', 'totals', 'shading'.


        .. figure::  upsetplot.png
           :align:   center
           :scale: 70 %

           Example plot of upset_plot()
    """

    upset_df = _generate_upset_srs(_fetch_sets(objs=objs, ref=ref))
    upset = upsetplot.plot(upset_df)
    plt.title(title)
    return upset


def venn_diagram(objs: Dict[str, Union[str, FeatureSet, Set[str]]], title: str = 'default', ref: str = 'predefined',
                 set_colors: tuple = ('r', 'g', 'b'),
                 alpha: float = 0.4, weighted: bool = True, lines: bool = True, linecolor: str = 'black',
                 linestyle='solid', linewidth=2.0,
                 normalize_to: float = 1.0):
    """
    Generate a Venn diagram of 2 to 3 sets, FeatureSets or attributes from the Attribute Reference Table.

    :param objs: the FeatureSets, python sets or user-defined attributes to plot.
    :type objs: a dictionary with 2-3 entries, where the keys are the names of the sets, and the values are either \
    a FeatureSet, a python set of feature indices, or a name of a column in the Attribute Reference Table. \
    For example: \
    {'first set':{'gene1','gene2','gene3'}, 'second set':'name_of_attribute_from_reference_table'}
    :type title: str
    :param set_colors: determines the colors of the circles in the diagram.
    :param ref: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type ref: str or pathlib.Path (default 'predefined')
    :param title: determines the title of the plot.
    :type set_colors: tuple of matplotlib-format colors, the same size as 'objs'
    :param alpha: determines the opacity of the circles.
    :type alpha: a float between 0 and 1
    :param weighted: if True, the plot will be area-weighted.
    :type weighted: bool (default True)
    :param lines: if True, adds an outline to the circles.
    :type lines: bool (default True)
    :param linecolor: Determines the color of the circles' outline.
    :type linecolor: matplotlib-format color (default 'black')
    :param linestyle: the style of the circles' outline.
    :type linestyle: 'solid' or 'dashed' (default 'solid')
    :param linewidth: the widdth of the circles' outlines.
    :type linewidth: float (default 2.0)
    :param normalize_to:
    :type normalize_to: float (default 1.0)
    :return: a tuple of a VennDiagram object; and a list of 2-3 Circle patches.


        .. figure::  venn.png
           :align:   center
           :scale: 70 %

           Example plot of venn_diagram()
    """
    if len(objs) > 3 or len(objs) < 2:
        raise ValueError(f'Venn can only accept between 2 and 3 sets. Instead got {len(objs)}')
    assert isinstance(title, str), f'Title must be a string. Instead got {type(title)}'
    objs = _fetch_sets(objs=objs, ref=ref)
    if len(objs) == 2:
        func = vn.venn2 if weighted else vn.venn2_unweighted
        func_circles = vn.venn2_circles
        set_colors = set_colors[0:2]
    else:
        func = vn.venn3 if weighted else vn.venn3_unweighted
        func_circles = vn.venn3_circles
    _ = plt.figure()
    plot_obj = func(tuple(objs.values()), tuple(objs.keys()), set_colors=set_colors, alpha=alpha,
                    normalize_to=normalize_to)
    if lines and weighted:
        circle_obj = func_circles(tuple(objs.values()), color=linecolor, linestyle=linestyle, linewidth=linewidth,
                                  normalize_to=normalize_to)
    elif lines and not weighted:
        warnings.warn('Cannot draw lines on an unweighted venn diagram. ')
        circle_obj = None
    else:
        circle_obj = None
    if title == 'default':
        title = 'Venn diagram of ' + ''.join([name + ' ' for name in objs.keys()])
    plt.title(title)
    return plot_obj, circle_obj


def _generate_upset_srs(objs: dict):
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
