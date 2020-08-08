"""
This module can perform enrichment analyses on a given set of genomic features and visualize their intersections. \
These include gene ontology/tissue/phenotype enrichment, enrichment for user-defined attributes, \
set visualization ,etc. \
Results of enrichment analyses can be saved to .csv files.
"""
import warnings
from itertools import compress, repeat, chain
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Set, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib_venn as vn
from numba import jit
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as multitest
import upsetplot as upset
import collections
from ipyparallel import Client
from matplotlib.cm import ScalarMappable
from scipy.stats import hypergeom, ttest_1samp
from statsmodels.stats.descriptivestats import sign_test
from xlmhg import get_xlmhg_test_result as xlmhg_test

from rnalysis.utils import io, parsing, validation, ref_tables
from rnalysis.filtering import Filter


class FeatureSet:
    """ Receives a filtered gene set and the set's name (optional) and preforms various enrichment analyses on them. """
    __slots__ = {'gene_set': 'set of feature names/indices', 'set_name': 'name of the FeatureSet'}
    __goa_queries__ = {}
    __go_basic__ = None

    # TODO: add a half-way memoization of goa_queries (before the table stage)

    def __init__(self, gene_set: Union[List[str], Set[str], Filter] = None, set_name: str = ''):

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
            self.gene_set = parsing.from_string(
                "Please insert genomic features/indices separated by newline \n"
                "(example: \n'WBGene00000001\nWBGene00000002\nWBGene00000003')", delimiter='\n')
        elif isinstance(gene_set, set):
            pass
        elif isinstance(gene_set, (list, tuple)):
            gene_set = set(gene_set)
        elif validation.isinstanceinh(gene_set, Filter):
            gene_set = gene_set.index_set
        else:
            raise TypeError(f"Error: 'gene_set' must be a set, list or tuple! Is a {type(gene_set)} instead. ")
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

    def _set_ops(self, others, op: Callable):

        """
        Performs a given set operation on self and on another object (FeatureSet or set).
        :type others: FeatureSet, set or str
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

    def union(self, *others):

        """
         Calculates the set union of the indices from multiple FeatureSet objects \
        (the indices that exist in at least one of the FeatureSet objects).

        :type others: FeatureSet, RankedSet or set
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

    def intersection(self, *others):

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

    def difference(self, *others):

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

    def symmetric_difference(self, other):

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
        return FeatureSet(self._set_ops([other], set.symmetric_difference))

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
        io.save_csv(df, filename=fname + '.csv')

    def _go_enrichment_fetch_annotations(self, organism: Union[str, int], gene_id_type: str,
                                         aspects: Union[str, Iterable[str]],
                                         evidence_types: Union[str, Iterable[str]],
                                         excluded_evidence_types: Union[str, Iterable[str]],
                                         databases: Union[str, Iterable[str]],
                                         excluded_databases: Union[str, Iterable[str]],
                                         qualifiers: Union[str, Iterable[str]],
                                         excluded_qualifiers: Union[str, Iterable[str]],
                                         propagate_annotations: str):
        goa_query_key = (organism, gene_id_type, parsing.data_to_tuple(evidence_types), parsing.data_to_tuple(aspects),
                         parsing.data_to_tuple(databases), parsing.data_to_tuple(qualifiers),
                         parsing.data_to_tuple(excluded_qualifiers), propagate_annotations.lower() != 'no')

        if goa_query_key in self.__goa_queries__:
            goa_df, go_id_to_term_dict = self.__goa_queries__[goa_query_key]
        else:
            taxon_id, organism_name = io.map_taxon_id(organism)
            if propagate_annotations != 'no':
                out = f"Fetching and propagating GO annotations for organism '{organism_name}' (taxon ID:{taxon_id})."
            else:
                out = f"Fetching GO annotations for organism '{organism_name}' (taxon ID:{taxon_id})."
            print(out)
            annotations_iter = io.golr_annotations_iterator(taxon_id, aspects, evidence_types, excluded_evidence_types,
                                                            databases, excluded_databases, qualifiers,
                                                            excluded_qualifiers)
            sparse_annotation_dict = {}
            go_id_to_term_dict = {}
            source_to_gene_id_dict = {}
            for annotation in annotations_iter:
                # add annotation to annotation dictionary
                if annotation['bioentity_internal_id'] not in sparse_annotation_dict:
                    sparse_annotation_dict[annotation['bioentity_internal_id']] = (set())
                sparse_annotation_dict[annotation['bioentity_internal_id']].add(annotation['annotation_class'])
                if propagate_annotations != 'no':
                    sparse_annotation_dict[annotation['bioentity_internal_id']].update(
                        annotation['isa_partof_closure_map'].keys())

                # add gene id and source to source dict
                if annotation['source'] not in source_to_gene_id_dict:
                    source_to_gene_id_dict[annotation['source']] = set()
                source_to_gene_id_dict[annotation['source']].add(annotation['bioentity_internal_id'])
                # add go term to term dictionary
                if annotation['annotation_class'] not in go_id_to_term_dict:
                    try:
                        go_id_to_term_dict[annotation['annotation_class']] = annotation['annotation_class_label']
                    except KeyError:
                        go_id_to_term_dict[annotation['annotation_class']] = 'None'
                    if propagate_annotations != 'no':
                        go_id_to_term_dict.update(annotation['isa_partof_closure_map'])
            print(f"Found annotations for {len(sparse_annotation_dict)} genes.")
            # translate gene IDs
            translated_sparse_annotation_dict = {}
            for source in source_to_gene_id_dict:
                translator = io.map_gene_ids(source_to_gene_id_dict[source], source, gene_id_type)
                for gene_id in sparse_annotation_dict.copy():
                    if gene_id in translator:
                        translated_sparse_annotation_dict[translator[gene_id]] = sparse_annotation_dict.pop(gene_id)

            # get boolean DataFrame for enrichment
            print("Generating Gene Ontology Reference Table...")
            goa_df = parsing.sparse_dict_to_bool_df(translated_sparse_annotation_dict)

            self.__goa_queries__[goa_query_key] = (goa_df, go_id_to_term_dict)

        return goa_df, go_id_to_term_dict

    def go_enrichment(self, organism: Union[str, int], gene_id_type: str = 'UniProtKB', fdr: float = 0.05,
                      biotype: str = 'protein_coding', background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                      biotype_ref_path: str = 'predefined',
                      propagate_annotations: str = 'elim',
                      aspects: Union[str, Iterable[str]] = 'any',
                      evidence_types: Union[str, Iterable[str]] = 'any',
                      excluded_evidence_types: Union[str, Iterable[str]] = None,
                      databases: Union[str, Iterable[str]] = 'any',
                      excluded_databases: Union[str, Iterable[str]] = None,
                      qualifiers: Union[str, Iterable[str]] = 'any',
                      excluded_qualifiers: Union[str, Iterable[str]] = 'not',
                      return_nonsignificant: bool = False,
                      save_csv: bool = False, fname=None, return_fig: bool = False, plot_horizontal: bool = True):
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
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
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
        :param propagate_annotations: determines the propagation method of annotations. If 'elim' (default), #TODO
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
        goa_df, go_id_to_term_dict = self._go_enrichment_fetch_annotations(organism, gene_id_type, aspects,
                                                                           evidence_types, excluded_evidence_types,
                                                                           databases, excluded_databases, qualifiers,
                                                                           excluded_qualifiers, propagate_annotations)

        biotype_ref_path = ref_tables.get_biotype_ref_path(biotype_ref_path)
        goa_df, gene_set = self._enrichment_get_reference(biotype, background_genes, goa_df, biotype_ref_path,
                                                          reindex=True)
        goa_df.fillna(False)

        print(f"Calculating enrichment using the '{propagate_annotations}' method...")
        args = [gene_set, goa_df, go_id_to_term_dict]
        if propagate_annotations in {'classic', 'no'}:
            res_dict = self._go_classic_pvalues(*args)
        elif propagate_annotations in {'elim', 'weight', 'all.m'}:
            dag_tree = io.fetch_go_basic() if self.__go_basic__ is None else self.__go_basic__
            if propagate_annotations == 'elim':
                args.extend([dag_tree, fdr])
                res_dict = self._go_elim_pvalues(*args)
            elif propagate_annotations == 'weight':
                args.append(dag_tree)
                res_dict = self._go_weight_pvalues(*args)
            else:
                args.extend([dag_tree, fdr])
                res_dict = self._go_allm_pvalues(*args)
        else:
            raise ValueError(f"Invalid value for 'propagate_annotations': '{propagate_annotations}'.")
        return self._go_enrichment_output(res_dict, fdr, return_nonsignificant, save_csv, fname, True, return_fig,
                                          plot_horizontal)

    @staticmethod
    def _calc_go_stats(go_id: str, bg_size: int, go_size: int, de_size: int, go_de_size: int, term: str,
                       res_dict: dict) -> None:
        obs, exp = go_de_size, de_size * (go_size / bg_size)
        log2_fold_enrichment = np.log2(obs / exp) if obs > 0 else -np.inf
        pvalue = FeatureSet._calc_hypergeometric_pval(bg_size, go_size, de_size, go_de_size)
        res_dict[go_id] = (term, de_size, obs, exp, log2_fold_enrichment, pvalue)

    @staticmethod
    def _go_classic_pvalues(gene_set: set, goa_df: pd.DataFrame, go_id_to_term_dict: Dict[str, str]) -> dict:
        res_dict = {}
        go_sizes = goa_df.sum(axis=0)
        bg_size = goa_df.shape[0]
        de_size = len(gene_set)
        go_de_sizes = goa_df.loc[gene_set].sum(axis=0)

        for go_id in goa_df.columns:
            go_size = go_sizes[go_id]
            go_de_size = go_de_sizes[go_id]
            FeatureSet._calc_go_stats(go_id, bg_size, go_size, de_size, go_de_size, go_id_to_term_dict[go_id], res_dict)

        return res_dict

    @staticmethod
    def _go_elim_pvalues(gene_set: set, goa_df: pd.DataFrame, go_id_to_term_dict: Dict[str, str],
                         dag_tree: parsing.DAGTreeParser, fdr: float, inplace: bool = False) -> dict:
        if not inplace:
            goa_df = goa_df.copy(deep=True)

        threshold = fdr / sum([len(level) for level in dag_tree.levels])
        res_dict = {}
        marked_nodes = {}
        bg_size = goa_df.shape[0]
        de_size = len(gene_set)
        for go_id in dag_tree.level_iterator():
            if go_id not in goa_df.columns:  # skip any GO ID that has no annotations whatsoever (direct or inherited)
                continue
            if go_id in marked_nodes:  # if this node was marked, remove from it all marked genes
                goa_df[go_id].loc[marked_nodes[go_id]] = 0
            go_size = goa_df[go_id].sum()
            go_de_size = goa_df[go_id].loc[gene_set].sum()
            FeatureSet._calc_go_stats(go_id, bg_size, go_size, de_size, go_de_size, go_id_to_term_dict[go_id], res_dict)

            if res_dict[go_id][-1] <= threshold:  # if current GO ID is significant, mark its ancestors
                new_marked_genes = set(goa_df[go_id][goa_df[go_id] == 1].index)
                for ancestor in dag_tree.upper_induced_graph_iterator(go_id):
                    if ancestor not in marked_nodes:
                        marked_nodes[ancestor] = set()
                    marked_nodes[ancestor] = marked_nodes[ancestor].union(new_marked_genes)
        return res_dict

    @staticmethod
    def _go_weight_pvalues(gene_set: set, goa_df: pd.DataFrame, go_id_to_term_dict: Dict[str, str],
                           dag_tree: parsing.DAGTreeParser, inplace: bool = False) -> dict:
        if not inplace:
            goa_df = goa_df.copy(deep=True)

        res_dict = {}
        weights = collections.defaultdict(lambda: 1)  # default weight for all nodes is 1
        bg_size = goa_df.shape[0]
        de_size = len(gene_set)
        for go_id in dag_tree.level_iterator():
            if go_id not in goa_df.columns:  # skip any GO ID that has no annotations whatsoever (direct or inherited)
                continue
            children = {child for child in dag_tree.go_terms[go_id].get_children() if child in goa_df.columns}
            FeatureSet._compute_term_sig(go_id, children, dag_tree, gene_set, goa_df, go_id_to_term_dict, weights,
                                         res_dict, bg_size, de_size)

        return res_dict

    @staticmethod
    def _compute_term_sig(go_id: str, children: set, dag_tree: parsing.DAGTreeParser, gene_set: set,
                          goa_df: pd.DataFrame, go_id_to_term_dict: Dict[str, str], weights: dict, res_dict: dict,
                          bg_size: int, de_size: int) -> None:
        # calculate stats for go_id
        go_size = np.ceil(goa_df[go_id].sum())
        go_de_size = np.ceil(goa_df[go_id].loc[gene_set].sum())
        FeatureSet._calc_go_stats(go_id, bg_size, go_size, de_size, go_de_size, go_id_to_term_dict[go_id], res_dict)

        if len(children) == 0:
            return

        sig_children = set()
        for child in children:
            weights[child] = res_dict[child][-1] / res_dict[go_id][-1]  # TODO: test potential inf issues
            if weights[child] >= 1:
                sig_children.add(child)

        # CASE 1: if go_id is more significant than all children, re-weigh the children and recompute their stats
        if len(sig_children) == 0:
            for child in children:
                goa_df[child] *= weights[child]
                go_size = np.ceil(goa_df[child].sum())
                go_de_size = np.ceil(goa_df[child].loc[gene_set].sum())
                FeatureSet._calc_go_stats(child, bg_size, go_size, de_size, go_de_size, go_id_to_term_dict[child],
                                          res_dict)
            return

        # CASE 2: if some children are more significant than parent, re-weigh ancesctors (including 'go_id'),
        # and then recompute stats for go_id
        for sig_child in sig_children:
            for inclusive_ancestor in chain([go_id], dag_tree.upper_induced_graph_iterator(go_id)):
                goa_df[inclusive_ancestor] *= 1 / weights[sig_child]
        # re-run compute_term_sig, only with the children which were not more significant than their parents
        FeatureSet._compute_term_sig(go_id, children.difference(sig_children), dag_tree, gene_set, goa_df,
                                     go_id_to_term_dict, weights, res_dict, bg_size, de_size)

    @staticmethod
    def _go_allm_pvalues(gene_set: set, goa_df: pd.DataFrame, go_id_to_term_dict: Dict[str, str],
                         dag_tree: parsing.DAGTreeParser, fdr: float) -> dict:
        classic = FeatureSet._go_classic_pvalues(gene_set, goa_df, go_id_to_term_dict)
        elim = FeatureSet._go_elim_pvalues(gene_set, goa_df, go_id_to_term_dict, dag_tree, fdr, inplace=False)
        weight = FeatureSet._go_weight_pvalues(gene_set, goa_df, go_id_to_term_dict, dag_tree, inplace=False)
        # TODO: make sure no in-place shenanigans occur with goa_df
        res_dict = {}
        for go_id in classic.keys():
            pvalue = np.exp(np.mean(np.log([classic[go_id[-1]], elim[go_id[-1]], weight[go_id[-1]]])))
            res_dict[go_id] = (*classic[go_id][:-1], pvalue)
        return res_dict

    @staticmethod
    def _single_enrichment(gene_set, attribute, attr_ref_df: pd.DataFrame, reps: int):
        assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
        srs = attr_ref_df[attribute]
        bg_array = np.isnan(srs.values)
        obs_array = np.isnan(srs.loc[gene_set].values)
        n = len(gene_set)
        expected_fraction = (bg_array.shape[0] - np.sum(bg_array)) / bg_array.shape[0]
        observed_fraction = (n - np.sum(obs_array)) / n
        log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
        pval = FeatureSet._calc_randomization_pval(n, log2_fold_enrichment, bg_array, reps, observed_fraction)

        return [attribute, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval]

    @staticmethod
    @jit(nopython=True)
    def _calc_randomization_pval(n: int, log2fc: float, bg_array: np.ndarray, reps: int, obs_frac: float):
        ind_range = np.arange(bg_array.shape[0])
        success = 0
        if log2fc >= 0:
            for _ in range(reps):
                success += (n - np.sum(bg_array[np.random.choice(ind_range, n, replace=False)])) / n >= obs_frac

        else:
            for _ in range(reps):
                success += (n - np.sum(bg_array[np.random.choice(ind_range, n, replace=False)])) / n <= obs_frac
        pval = (success + 1) / (reps + 1)
        return pval

    @staticmethod
    def _enrichment_get_attrs(attributes, attr_ref_path):
        if attributes is None:
            attributes = parsing.from_string(
                "Please insert attributes separated by newline "
                "(for example: \n'epigenetic_related_genes\nnrde-3 targets\nALG-3/4 class small RNAs')")
        elif isinstance(attributes, (str, int)):
            attributes = [attributes]
        else:
            assert isinstance(attributes, (list, tuple, set)), "'attributes' must be a list, tuple or set!"
            for attr in attributes:
                if isinstance(attr, int):
                    assert attr >= 0, f"Error in attribute number {attr}: index must be non-negative!"
                else:
                    assert isinstance(attr, str), f"Invalid type of attribute {attr}: {type(attr)}"

        try:
            with open(attr_ref_path) as f:
                all_attrs = f.readline().split(',')[1::]
        except FileNotFoundError:
            raise FileNotFoundError(f"Invalid or nonexistent Attribute Reference Table path! path:'{attr_ref_path}'")
        if all_attrs[-1].endswith('\n'):
            all_attrs[-1] = all_attrs[-1][:-1]

        if attributes == ['all']:
            attributes = all_attrs
        elif np.all([True if isinstance(i, int) else False for i in attributes]):
            return [all_attrs[ind] for ind in attributes]
        return attributes

    def _enrichment_get_reference(self, biotype: Union[str, List[str], Set[str], Tuple[str]], background_genes,
                                  attr_ref_df: pd.DataFrame, biotype_ref_path: str, reindex: bool = False):
        gene_set = self.gene_set

        assert (isinstance(biotype, (str, list, set, tuple)))

        if background_genes is not None:
            if isinstance(background_genes, FeatureSet):
                bg = background_genes.gene_set
            elif validation.isinstanceinh(background_genes, Filter):
                bg = background_genes.index_set
            elif isinstance(background_genes, set):
                bg = background_genes
            else:
                raise TypeError(f"background_genes must be a set, enrichment.FeatureSet or filtering.Filter; "
                                f"instead got {type(background_genes)}")

            if biotype != 'all':
                warnings.warn("both 'biotype' and 'background_genes' were specified. Therefore 'biotype' is ignored.")

        else:
            if biotype == 'all':
                bg = attr_ref_df.index
            else:
                biotype_ref_df = io.load_csv(biotype_ref_path)
                validation.validate_biotype_table(biotype_ref_df)
                biotype_ref_df.set_index('gene', inplace=True)
                biotype = parsing.data_to_list(biotype)
                mask = pd.Series(np.zeros_like(biotype_ref_df['biotype'].values, dtype=bool),
                                 biotype_ref_df['biotype'].index,
                                 name='biotype')
                for bio in biotype:
                    mask = mask | (biotype_ref_df['biotype'] == bio)
                bg = biotype_ref_df[mask].index

        if reindex:
            attr_ref_df = attr_ref_df.reindex(attr_ref_df.index.union(bg))
        else:
            attr_ref_df = attr_ref_df.loc[attr_ref_df.index.intersection(bg)]
        if len(attr_ref_df.index) < len(bg):
            warnings.warn(
                f"{len(bg) - len(attr_ref_df.index)} indices from the requested "
                f"background genes do not appear in the Attribute Reference Table, and are therefore ignored. \n"
                f"This leaves a total of {len(attr_ref_df.index)} background genes. ")
        print(f"{len(attr_ref_df.index)} background genes are used. ")

        not_in_bg = gene_set.difference(set(attr_ref_df.index))
        if len(not_in_bg) > 0:
            gene_set = gene_set.difference(not_in_bg)
            warnings.warn(f"{len(not_in_bg)} genes in the enrichment set do not appear in the "
                          f"Attribute Reference Table and/or the background genes. \n"
                          f"Enrichment will be run on the remaining {len(gene_set)} genes.")
        attr_ref_df.sort_index(inplace=True)
        return attr_ref_df, gene_set

    def _enrichment_setup(self, biotype: str, background_genes: Union[Iterable[str], str, Iterable[int], int, None],
                          attr_ref_path: str, biotype_ref_path: str, attributes: Union[str, List[str], List[int]]):
        """
        Perform setup for enrichment functions. This function receives most of the input variables \
        from enrichment functions, generates a full list of attributes for enrichment, gets the relevant full path to \
        Attribute/Biotype reference tables, loads the Reference tables into a DataFrame, \
        filters the Attribute DataFrame to include only relevant Attributes, \
        generates a background gene set according to the user's specifications, \
        filters the tested gene set based on the background gene set and/or biotype, \
        and standardizes the data scale input.

        """
        attr_ref_path = ref_tables.get_attr_ref_path(attr_ref_path)
        biotype_ref_path = ref_tables.get_biotype_ref_path(biotype_ref_path)

        attr_ref_df = io.load_csv(attr_ref_path)
        validation.validate_attr_table(attr_ref_df)
        attr_ref_df.set_index('gene', inplace=True)

        attr_ref_df, gene_set = self._enrichment_get_reference(biotype=biotype, background_genes=background_genes,
                                                               attr_ref_df=attr_ref_df,
                                                               biotype_ref_path=biotype_ref_path)

        attributes = self._enrichment_get_attrs(attributes=attributes, attr_ref_path=attr_ref_path)
        return attr_ref_df, gene_set, attributes

    def _go_enrichment_output(self, res_dict: dict, fdr: float, return_nonsignificant: bool, save_csv: bool,
                              fname: str, plot: bool, return_fig: bool, plot_horizontal: bool,
                              single_list: bool = False):
        """
        Formats the enrich list into a results Dataframe, saves the DataFrame to csv if requested, \
        plots the enrichment results, and returns either the Dataframe alone or the Dataframe and the Figure object.
        Called at the end of every enrichment function \
        (enrich_randomization, enrich_randomization_parallel, enrich_statistic...).

        """
        n_plot = min(10, len(res_dict))
        if single_list:
            en_score_col = 'log2_enrichment_score'
            columns = ['go_id', 'term', 'samples', en_score_col, 'pval']
            title = f"Single-list enrichment for {self.set_name}\ntop {n_plot} most significant GO terms"
            ylabel = r"$\log_2$(XL-mHG enrichment score)"
        else:
            en_score_col = 'log2_fold_enrichment'
            columns = ['term', 'samples', 'obs', 'exp', en_score_col, 'pval']
            title = f"Enrichment for {self.set_name}\ntop {n_plot} most significant GO terms"
            ylabel = r"$\log_2$(Fold Enrichment)"
        res_df = pd.DataFrame.from_dict(res_dict, orient='index', columns=columns)
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.rename_axis('go_id')
        res_df.sort_values('padj', inplace=True)
        if not return_nonsignificant:
            res_df = res_df[res_df['significant']]

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        if plot:
            n_plot = min(10, res_df.shape[0])
            fig = self.plot_enrichment_results(res_df.iloc[:n_plot], en_score_col, res_df.columns[0],
                                               title=title, ylabel=ylabel, plot_horizontal=plot_horizontal)
            if return_fig:
                return res_df, fig
        return res_df

    def _enrichment_output(self, enriched_list: list, fdr: float, save_csv: bool, fname: str, plot: bool,
                           return_fig: bool, plot_horizontal: bool = True):
        """
        Formats the enrich list into a results Dataframe, saves the DataFrame to csv if requested, \
        plots the enrichment results, and returns either the Dataframe alone or the Dataframe and the Figure object.
        Called at the end of every enrichment function \
        (enrich_randomization, enrich_randomization_parallel, enrich_statistic...).

        """
        res_df = pd.DataFrame(enriched_list,
                              columns=['name', 'samples', 'obs', 'exp', 'log2_fold_enrichment', 'pval'])
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.set_index('name', inplace=True)
        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        if plot:
            fig = self.plot_enrichment_results(res_df, title=f"Enrichment for {self.set_name}",
                                               plot_horizontal=plot_horizontal)
            if return_fig:
                return res_df, fig
        return res_df

    def enrich_randomization_parallel(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                                      fdr: float = 0.05, reps: int = 10000,
                                      biotype: str = 'protein_coding',
                                      background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                                      attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                                      save_csv: bool = False, fname=None, return_fig: bool = False,
                                      plot_horizontal: bool = True, random_seed: int = None):

        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set \
        using a randomization test running in parallel. The attributes are drawn from an Attribute Reference Table. \
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
        attr_ref_df, gene_set, attributes, = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)

        client = Client()
        dview = client[:]
        dview.execute("""import numpy as np
              import pandas as pd""")
        if random_seed is not None:
            assert isinstance(random_seed, int) and random_seed >= 0, f"random_seed must be a non-negative integer. " \
                                                                      f"Value {random_seed} invalid."
            dview.execute(f"np.random.seed({random_seed})")
        k = len(attributes)
        gene_set_rep = list(repeat(gene_set, k))
        attr_ref_df_rep = list(repeat(attr_ref_df, k))
        reps_rep = list(repeat(reps, k))

        res = dview.map(FeatureSet._single_enrichment, gene_set_rep, attributes, attr_ref_df_rep, reps_rep)
        enriched_list = res.result()
        return self._enrichment_output(enriched_list, fdr, save_csv, fname, True, return_fig, plot_horizontal)

    def enrich_randomization(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                             reps: int = 10000, biotype: str = 'protein_coding',
                             background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                             attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                             save_csv: bool = False, fname=None, return_fig: bool = False, plot_horizontal: bool = True,
                             random_seed: int = None):

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

           Example plot of enrich_randomization()


        .. figure::  plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of enrich_randomization(plot_horizontal = False)

        """
        attr_ref_df, gene_set, attributes = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)
        enriched_list = []
        if random_seed is not None:
            assert isinstance(random_seed, int) and random_seed >= 0, f"random_seed must be a non-negative integer. " \
                                                                      f"Value {random_seed} invalid."
            np.random.seed(random_seed)

        for k, attribute in enumerate(attributes):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            enriched_list.append(self._single_enrichment(gene_set, attribute, attr_ref_df, reps))
            print(f"Finished {k + 1} attributes out of {len(attributes)}")

        return self._enrichment_output(enriched_list, fdr, save_csv, fname, True, return_fig, plot_horizontal)

    def enrich_hypergeometric(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                              biotype: str = 'protein_coding',
                              background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                              attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                              save_csv: bool = False, fname=None, return_fig: bool = False,
                              plot_horizontal: bool = True):

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
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
        :type attr_ref_path: str or pathlib.Path (default 'predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default 'predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
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
        attr_ref_df, gene_set, attributes = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)
        enriched_list = []
        bg_size = attr_ref_df.shape[0]
        de_size = len(gene_set)
        for k, attribute in enumerate(attributes):
            srs = attr_ref_df[attribute]
            go_size = srs.notna().sum()
            go_de_size = srs.loc[gene_set].notna().sum()

            expected_fraction = go_size / bg_size
            observed_fraction = go_de_size / de_size
            log2_fold_enrichment = np.log2(
                observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
            pval = self._calc_hypergeometric_pval(bg_size=bg_size, go_size=go_size,
                                                  de_size=de_size, go_de_size=go_de_size)
            obs, exp = int(de_size * observed_fraction), de_size * expected_fraction

            enriched_list.append(
                (attribute, de_size, obs, exp, log2_fold_enrichment, pval))

        return self._enrichment_output(enriched_list, fdr, save_csv, fname, True, return_fig, plot_horizontal)

    @staticmethod
    def _calc_hypergeometric_pval(bg_size: int, go_size: int, de_size: int, go_de_size: int):

        """
        Performs a hypergeometric test on the given enrichment set. \
        Given M genes in the background set, n genes in the test set, \
        with N genes from the background set belonging to a specific attribute (or 'success') \
        and X genes from the test set belonging to that attribute. \
        If we were to randomly draw n genes from the background set (without replacement), \
        what is the probability of drawing X or more (in case of enrichment)/X or less (in case of depletion) \
        genes belonging to the given attribute?

        :param bg_size: size of the background set. Usually denoted as 'M'.
        :type bg_size: positive int
        :param go_size: number of features in the background set corresponding to the attribute, \
        or number of successes in the population. Usually denoted as 'n'.
        :type go_size: positive int
        :param de_size: size of the differentially-expressed set, or size of test set. usually denoted as 'N'.
        :type de_size: positive int
        :param go_de_size: or number of successes in the test set. Usually denoted as 'x' or 'k'. s
        :type go_de_size: non-negative int
        :return: p-value of the hypergeometric test.
        :rtype: float between 0 and 1

        """
        if go_de_size / de_size < go_size / bg_size:
            return hypergeom.cdf(go_de_size, bg_size, go_size, de_size)
        return hypergeom.sf(go_de_size - 1, bg_size, go_size, de_size)

    @staticmethod
    def plot_enrichment_results(df: pd.DataFrame, en_score_col: str = 'log2_fold_enrichment', name_col: str = None,
                                title: str = '', ylabel: str = r"$\log_2$(Fold Enrichment)",
                                plot_horizontal: bool = True, center_bars: bool = True):

        """
        Receives a DataFrame output from an enrichment function and plots it in a bar plot. \
        For the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

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
            figsize = [5.6, 0.4 * (6.4 + df.shape[0])]
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
            asterisks, fontweight = FeatureSet._get_pval_asterisk(sig)
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

    @staticmethod
    def _get_pval_asterisk(pval: float):
        fontweight = 'bold'
        if pval < 0.0001:
            asterisks = u'\u2217' * 4
        elif pval < 0.001:
            asterisks = u'\u2217' * 3
        elif pval < 0.01:
            asterisks = u'\u2217' * 2
        elif pval < 0.05:
            asterisks = u'\u2217'
        else:
            asterisks = 'ns'
            fontweight = 'normal'
        return asterisks, fontweight

    def enrich_non_categorical(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                               fdr: float = 0.05, parametric_test: bool = False, biotype: str = 'protein_coding',
                               background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                               attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                               plot_log_scale: bool = True, plot_style: str = 'overlap', n_bins: int = 50,
                               save_csv: bool = False, fname=None, return_fig: bool = False):
        attr_ref_df, gene_set, attributes = \
            self._enrichment_setup(biotype, background_genes, attr_ref_path, biotype_ref_path, attributes)
        enriched_list = []
        for k, attribute in enumerate(attributes):
            srs = attr_ref_df[attribute]
            if not parametric_test:
                exp, obs = srs.median(), srs[gene_set].median()
                _, pval = sign_test(srs[gene_set].values, exp)
            else:
                exp, obs = srs.mean(), srs[gene_set].mean()
                _, pval = ttest_1samp(srs[gene_set], popmean=exp, nan_policy='propagate')

            enriched_list.append(
                (attribute, len(gene_set), obs, exp, np.nan, pval))

        results_df = self._enrichment_output(enriched_list, fdr, save_csv, fname, return_fig, False)
        results_df.dropna(axis=1, inplace=True)
        for attribute, pval in zip(attributes, results_df['padj']):
            self._enrichment_plot_histogram(attribute, gene_set, self.set_name, attr_ref_df[attribute], plot_log_scale,
                                            plot_style, n_bins, parametric_test, pval)
        return results_df

    @staticmethod
    def _enrichment_plot_histogram(attribute: str, gene_set: set, set_name: str, attr_srs: pd.Series,
                                   plot_log_scale: bool, plot_style: str, n_bins: int, parametric_test: bool,
                                   pval: float):
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

    def _set_ops(self, others, op: Callable):
        warnings.warn("Warning: when performing set operations with RankedSet objects, "
                      "the return type will always be FeatureSet and not RankedSet.")
        super()._set_ops(others, op)

    def go_enrichment_single_list(self, organism: Union[str, int], gene_id_type: str = 'UniProtKB', fdr: float = 0.05,
                                  propagate_annotations: str = 'elim',
                                  aspects: Union[str, Iterable[str]] = 'any',
                                  evidence_types: Union[str, Iterable[str]] = 'any',
                                  excluded_evidence_types: Union[str, Iterable[str]] = None,
                                  databases: Union[str, Iterable[str]] = 'any',
                                  excluded_databases: Union[str, Iterable[str]] = None,
                                  qualifiers: Union[str, Iterable[str]] = 'any',
                                  excluded_qualifiers: Union[str, Iterable[str]] = 'not',
                                  return_nonsignificant: bool = False,
                                  save_csv: bool = False, fname=None,
                                  return_fig: bool = False, plot_horizontal: bool = True):
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
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
        :param propagate_annotations: determines the propagation method of annotations. If 'elim' (default), #TODO
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
        goa_df, go_id_to_term_dict = self._go_enrichment_fetch_annotations(organism, gene_id_type, aspects,
                                                                           evidence_types, excluded_evidence_types,
                                                                           databases, excluded_databases, qualifiers,
                                                                           excluded_qualifiers, propagate_annotations)

        print("Calculating enrichment...")
        goa_df, gene_set = self._enrichment_get_reference('all', self.gene_set, goa_df, '', reindex=True)
        goa_df.fillna(False)
        # TODO: implement propagate_annotations in go_enrichment_single_list()
        res_dict = {}
        for k, go_id in enumerate(goa_df.columns):
            pval, en_score = self._calc_xlmhg_stats(
                self._xlmhg_index_vector(self.ranked_genes, go_id, goa_df, notna=False), len(self.ranked_genes))
            log2_en_score = np.log2(en_score) if en_score > 0 else -np.inf
            res_dict[go_id] = [go_id_to_term_dict[go_id], len(self.ranked_genes), log2_en_score, pval]

        return self._go_enrichment_output(res_dict, fdr, return_nonsignificant, save_csv, fname, True, return_fig,
                                          plot_horizontal)

    def enrich_single_list(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                           attr_ref_path: str = 'predefined', save_csv: bool = False, fname=None,
                           return_fig: bool = False, plot_horizontal: bool = True):
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
        :type fdr: float between 0 and 1
        :param fdr: Indicates the FDR threshold for significance.
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
        # TODO: prevent printout of 'XYZ background genes are used' when doing single-list enrichment!!
        # TODO: look into floating point precision for XL-mHG p-values
        attr_ref_df, gene_set, attributes = self._enrichment_setup(biotype='all', background_genes=None,
                                                                   attr_ref_path=attr_ref_path, biotype_ref_path='',
                                                                   attributes=attributes)
        if len(gene_set) == len(self.ranked_genes):
            ranked_genes = self.ranked_genes
        else:
            ranked_genes = np.empty((len(gene_set),), dtype=self.ranked_genes.dtype)
            i = 0
            for elem in self.ranked_genes:
                if elem in gene_set:
                    ranked_genes[i] = elem
                    i += 1

        enriched_list = []
        for attribute in attributes:
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            pval, en_score = self._calc_xlmhg_stats(self._xlmhg_index_vector(ranked_genes, attribute, attr_ref_df),
                                                    len(ranked_genes))
            log2_en_score = np.log2(en_score) if en_score > 0 else -np.inf
            enriched_list.append([attribute, len(ranked_genes), log2_en_score, pval])

        return self._xlmhg_output(enriched_list, fdr, save_csv, fname, return_fig, plot_horizontal)

    def _xlmhg_output(self, enriched_list: list, fdr: float, save_csv: bool, fname: str, return_fig: bool,
                      plot_horizontal: bool):
        """
        Formats the enrich list into a results Dataframe, saves the DataFrame to csv if requested, \
        plots the enrichment results, and returns either the Dataframe alone or the Dataframe and the Figure object.
        Called at the end of enrich_single_list().

        """
        en_score_col = 'log2_enrichment_score'
        res_df = pd.DataFrame(enriched_list, columns=['name', 'samples', en_score_col, 'pval'])
        res_df.replace(-np.inf, -np.max(np.abs(res_df['log2_enrichment_score'].values)))
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.set_index('name', inplace=True)

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        fig = self.plot_enrichment_results(res_df, en_score_col=en_score_col,
                                           title=f"Single-list enrichment for {self.set_name}",
                                           ylabel=r"$\log_2$(XL-mHG enrichment score)", plot_horizontal=plot_horizontal)
        if return_fig:
            return res_df, fig
        return res_df

    @staticmethod
    def _calc_xlmhg_stats(index_vec: np.ndarray, ranked_genes_len: int):
        index_vec = np.uint16(index_vec)
        rev_index_vec = np.uint16([ranked_genes_len - 1 - index_vec[i - 1] for i in range(len(index_vec), 0, -1)])
        # X = the minimal amount of 'positive' elements above the hypergeometric cutoffs out of all of the positive
        # elements in the ranked set. Determined to be the minimum between x_min and ceil(x_frac * k),
        # where 'k' is the number of 'positive' elements in the ranked set.
        x_frac = 0.5
        x_min = 10
        # L = the lowest possible cutoff (n) to be tested out of the entire list.
        # Determined to be floor(l_frac * N), where 'N' is total number of elements in the ranked set (ranked_genes_len).
        l_frac = 0.1
        # pre-allocate empty array to speed up computation
        table = np.empty((len(index_vec) + 1, ranked_genes_len - len(index_vec) + 1), dtype=np.longdouble)
        res_obj_fwd = xlmhg_test(N=ranked_genes_len, indices=index_vec, L=int(np.floor(l_frac * ranked_genes_len)),
                                 X=min(x_min, int(np.ceil(x_frac * len(index_vec)))), table=table)
        res_obj_rev = xlmhg_test(N=ranked_genes_len, indices=rev_index_vec, L=int(np.floor(l_frac * ranked_genes_len)),
                                 X=min(x_min, int(np.ceil(x_frac * len(index_vec)))), table=table)

        if res_obj_fwd.pval <= res_obj_rev.pval:
            pval, en_score = res_obj_fwd.pval, res_obj_fwd.escore
        else:
            pval, en_score = res_obj_rev.pval, 1 / res_obj_rev.escore
        pval = pval if not np.isnan(pval) else 1
        en_score = en_score if not np.isnan(en_score) else 1
        return pval, en_score

    @staticmethod
    def _xlmhg_index_vector(ranked_genes, attribute, attr_ref_df, notna: bool = True):
        ranked_srs = attr_ref_df.loc[ranked_genes, attribute]
        assert ranked_srs.shape[0] == len(ranked_genes)
        if notna:
            return np.uint16(np.nonzero(ranked_srs.notna().values)[0])
        return np.uint16(np.nonzero(ranked_srs.values)[0])


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
    upsetplot = upset.plot(upset_df)
    plt.title(title)
    return upsetplot


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
    srs = pd.Series(index=multi_ind)
    for ind in multi_ind:
        sets = list(compress(names, ind))
        group_size = len(set.intersection(*[objs[s] for s in sets]))
        srs.loc[ind] = group_size
    return srs
