"""
This module can perform enrichment analyses on a given set of genomic features and visualize their intersections. \
These include gene ontology/tissue/phenotype enrichment, enrichment for user-defined attributes, \
set visualization ,etc. \
Results of enrichment analyses can be saved to .csv files.
"""
import functools
import itertools
import types
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple, Union, Sequence, Literal

import matplotlib
import matplotlib.pyplot as plt
import matplotlib_venn as vn
import numpy as np
import pandas as pd
import upsetplot

from rnalysis.filtering import Filter, readable_name
from rnalysis.utils import io, parsing, settings, validation, enrichment_runner, generic, param_typing, ontology
from rnalysis.utils.param_typing import GO_ASPECTS, GO_EVIDENCE_TYPES, GO_QUALIFIERS, DEFAULT_ORGANISMS, \
    PARALLEL_BACKENDS, get_gene_id_types, PositiveInt


class FeatureSet(set):
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
        # elif isinstance(gene_set,type(self)):
        #     gene_set = gene_set.gene_set
        else:
            gene_set = parsing.data_to_set(gene_set)
        self.gene_set = gene_set
        self.set_name = set_name

        super().__init__(self.gene_set)

    def __copy__(self):
        obj = type(self)(self.gene_set.copy(), self.set_name)
        return obj

    def __repr__(self):
        return f"{self.__class__.__name__}: '{self.set_name}'"

    def __eq__(self, other):
        if type(self) != type(other):
            return False

        if self.set_name != other.set_name:
            return False

        if self.gene_set != other.gene_set:
            return False

        return True

    def _update(self, **kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self, key, val)
            except AttributeError:
                raise AttributeError(f"Cannot update attribute {key} for {type(self)} object: attribute does not exist")

    def _inplace(self, func, func_kwargs, inplace: bool, **update_kwargs):
        """
        Executes the user's choice whether to filter in-place or create a new instance of the FeatureSet object.
        """
        filter_obj = self._convert_to_filter_obj()
        applied = func(filter_obj, **func_kwargs, inplace=inplace)
        if inplace:
            applied = filter_obj

        new_set = applied.index_set
        suffix = '_' + applied.fname.stem.split('_')[-1]
        new_name = self.set_name + suffix

        # if inplace, modify self, name and other properties of self
        if inplace:
            self.intersection_update(new_set)
            self._update(gene_set=new_set, set_name=new_name, **update_kwargs)
        # if not inplace, copy self, modify the self, name, and other properties of the copy, and return it
        else:
            new_obj = type(self)(new_set, new_name)
            new_obj._update(**update_kwargs)
            return new_obj

    def _convert_to_filter_obj(self) -> Filter:
        return Filter.from_dataframe(pd.DataFrame(index=parsing.data_to_list(self.gene_set)), self.set_name)

    @readable_name('Translate gene IDs')
    def translate_gene_ids(self, translate_to: Union[str, Literal[get_gene_id_types()]],
                           translate_from: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                           remove_unmapped_genes: bool = False, inplace: bool = True) -> 'FeatureSet':
        """
                Translates gene names/IDs from one type to another. \
                Mapping is done using the UniProtKB Gene ID Mapping service. \
                You can choose to optionally drop from the table all rows that failed to be translated.

                :param translate_to: the gene ID type to translate gene names/IDs to. \
                For example: UniProtKB, Ensembl, Wormbase.
                :type translate_to: str
                :param translate_from: the gene ID type to translate gene names/IDs from. \
                For example: UniProtKB, Ensembl, Wormbase. If translate_from='auto', \
                *RNAlysis* will attempt to automatically determine the gene ID type of the features in the table.
                :type translate_from: str or 'auto' (default='auto')
                :param remove_unmapped_genes: if True, rows with gene names/IDs that could not be translated \
                will be dropped from the table. \
                Otherwise, they will remain in the table with their original gene name/ID.
                :type remove_unmapped_genes: bool (default=False)
                :return: returns a new and translated FeatureSet.
    """
        kwargs = dict(translate_to=translate_to,
                      translate_from=translate_from,
                      remove_unmapped_genes=remove_unmapped_genes)
        return self._inplace(Filter.translate_gene_ids, kwargs, inplace)

    @readable_name('Filter by KEGG Pathways annotations')
    def filter_by_kegg_annotations(self, kegg_ids: Union[str, List[str]],
                                   mode: Literal['union', 'intersection'] = 'union',
                                   organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                                   gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                   opposite: bool = False, inplace: bool = True):
        """
        Filters genes according to KEGG pathways, keeping only genes that belong to specific KEGG pathway. \
        When multiple KEGG IDs are given, filtering can be done in 'union' mode \
        (where genes that belong to at least one pathway are not filtered out), or in 'intersection' mode \
        (where only genes that belong to ALL pathways are not filtered out).

        :param kegg_ids: the KEGG pathway IDs according to which the table will be filtered. \
        An example for a legal KEGG pathway ID would be 'path:cel04020' for the C. elegans calcium signaling pathway.
        :type kegg_ids: str or list of str
        :type mode: 'union' or 'intersection'.
        :param mode: If 'union', filters out every genomic feature that does not belong to one or more \
        of the indicated attributes. If 'intersection', \
        filters out every genomic feature that does not belong to ALL of the indicated attributes.
        param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :type opposite: bool (default=False)
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current FeatureSet object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new, filtered instance of the FeatureSet object.
        """

        kwargs = dict(kegg_ids=kegg_ids, mode=mode, organism=organism, gene_id_type=gene_id_type, opposite=opposite)
        return self._inplace(Filter.filter_by_kegg_annotations, kwargs, inplace)

    @readable_name('Filter by Gene Ontology (GO) annotation')
    def filter_by_go_annotations(self, go_ids: Union[str, List[str]], mode: Literal['union', 'intersection'] = 'union',
                                 organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                                 gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                 propagate_annotations: bool = True,
                                 evidence_types: Union[Literal[('any',) + GO_EVIDENCE_TYPES], Iterable[
                                     Literal[GO_EVIDENCE_TYPES]]] = 'any',
                                 excluded_evidence_types: Union[Literal[GO_EVIDENCE_TYPES], Iterable[Literal[
                                     GO_EVIDENCE_TYPES]]] = (),
                                 databases: Union[str, Iterable[str]] = 'any',
                                 excluded_databases: Union[str, Iterable[str]] = (),
                                 qualifiers: Union[
                                     Literal[('any',) + GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'any',
                                 excluded_qualifiers: Union[
                                     Literal[GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'not',
                                 opposite: bool = False, inplace: bool = True):
        """
        Filters genes according to GO annotations, keeping only genes that are annotated with a specific GO term. \
        When multiple GO terms are given, filtering can be done in 'union' mode \
        (where genes that belong to at least one GO term are not filtered out), or in 'intersection' mode \
        (where only genes that belong to ALL GO terms are not filtered out).

        :param go_ids:
        :type go_ids: str or list of str
        :type mode: 'union' or 'intersection'.
        :param mode: If 'union', filters out every genomic feature that does not belong to one or more \
        of the indicated attributes. If 'intersection', \
        filters out every genomic feature that does not belong to ALL of the indicated attributes.
        param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :param propagate_annotations: determines the propagation method of GO annotations. \
        'no' does not propagate annotations at all; 'classic' propagates all annotations up to the DAG tree's root; \
        'elim' terminates propagation at nodes which show significant enrichment; 'weight' performs propagation in a \
        weighted manner based on the significance of children nodes relatively to their parents; and 'allm' uses a \
        combination of all proopagation methods. To read more about the propagation methods, \
        see Alexa et al: https://pubmed.ncbi.nlm.nih.gov/16606683/
        :type propagate_annotations: 'classic', 'elim', 'weight', 'all.m', or 'no' (default='elim')
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
        :type excluded_qualifiers: str, Iterable of str, or None (default='not')
        :type opposite: bool (default=False)
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current FeatureSet object. If False, \
        the function will return a new FeatureSet instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new, filtered instance of the FeatureSet object.
        """
        kwargs = dict(go_ids=go_ids, mode=mode, organism=organism, gene_id_type=gene_id_type,
                      propagate_annotations=propagate_annotations, evidence_types=evidence_types,
                      excluded_evidence_types=excluded_evidence_types, databases=databases,
                      excluded_databases=excluded_databases, qualifiers=qualifiers,
                      excluded_qualifiers=excluded_qualifiers, opposite=opposite)
        return self._inplace(Filter.filter_by_go_annotations, kwargs, inplace)

    @readable_name('Filter by user-defined attribute')
    def filter_by_attribute(self, attributes: Union[str, List[str]] = None,
                            mode: Literal['union', 'intersection'] = 'union',
                            ref: Union[str, Path, Literal['predefined']] = 'predefined',
                            opposite: bool = False, inplace: bool = True):

        """
        Filters features according to user-defined attributes from an Attribute Reference Table. \
        When multiple attributes are given, filtering can be done in 'union' mode \
        (where features that belong to at least one attribute are not filtered out), or in 'intersection' mode \
        (where only features that belong to ALL attributes are not filtered out). \
        To learn more about user-defined attributes and Attribute Reference Tables, read the user guide.

        :type attributes: string or list of strings, \
        which are column titles in the user-defined Attribute Reference Table.
        :param attributes: attributes to filter by.
        :type mode: 'union' or 'intersection'.
        :param mode: If 'union', filters out every genomic feature that does not belong to one or more \
        of the indicated attributes. If 'intersection', \
        filters out every genomic feature that does not belong to ALL of the indicated attributes.
        :type ref: str or pathlib.Path (default='predefined')
        :param ref: filename/path of the attribute reference table to be used as reference.
        :type opposite: bool (default=False)
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new and filtered instance of the Filter object.
        """
        kwargs = dict(attributes=attributes, mode=mode, opposite=opposite, ref=ref)
        return self._inplace(Filter.filter_by_attribute, kwargs, inplace)

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
            for i, gene in enumerate(self.gene_set):
                line = gene + '\n' if (i + 1) < len(self.gene_set) else gene
                f.write(line)

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
                raise TypeError(f"'other' must be an FeatureSet object or a set, instead got {type(other)}")
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

    @readable_name("GO Enrichment")
    def go_enrichment(self, organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                      gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                      alpha: param_typing.Fraction = 0.05,
                      statistical_test: Literal['fisher', 'hypergeometric', 'randomization'] = 'fisher',
                      biotype: Union[str, List[str], Literal['all']] = 'all',
                      background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                      biotype_ref_path: Union[str, Path, Literal['predefined']] = 'predefined',
                      propagate_annotations: Literal['classic', 'elim', 'weight', 'all.m', 'no'] = 'elim',
                      aspects: Union[Literal[('any',) + GO_ASPECTS], Iterable[Literal[GO_ASPECTS]]] = 'any',
                      evidence_types: Union[
                          Literal[('any',) + GO_EVIDENCE_TYPES], Iterable[Literal[GO_EVIDENCE_TYPES]]] = 'any',
                      excluded_evidence_types: Union[Literal[GO_EVIDENCE_TYPES], Iterable[Literal[
                          GO_EVIDENCE_TYPES]]] = (),
                      databases: Union[str, Iterable[str], Literal['any']] = 'any',
                      excluded_databases: Union[str, Iterable[str]] = (),
                      qualifiers: Union[Literal[('any',) + GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'any',
                      excluded_qualifiers: Union[Literal[GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'not',
                      exclude_unannotated_genes: bool = True, return_nonsignificant: bool = False,
                      save_csv: bool = False, fname=None, return_fig: bool = False, plot_horizontal: bool = True,
                      show_expected: bool = False, plot_style: Literal['bar', 'lollipop'] = 'bar',
                      plot_ontology_graph: bool = True,
                      ontology_graph_format: Literal[param_typing.GRAPHVIZ_FORMATS] = 'none',
                      randomization_reps: PositiveInt = 10000, random_seed: Union[int, None] = None,
                      parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                      gui_mode: bool = False) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        """
        Calculates enrichment and depletion of the FeatureSet for Gene Ontology (GO) terms against a background set. \
        The GO terms and annotations are drawn via the GO Solr search engine GOlr, \
        using the search terms defined by the user. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. \
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
        :type gene_id_type: str or 'auto' (default='auto')
        :type alpha: float between 0 and 1 (default=0.05)
        :param alpha: Indicates the FDR threshold for significance.
        :param statistical_test: determines the statistical test to be used for enrichment analysis. \
        Note that some propagation methods support only some of the available statistical tests.
        :type statistical_test: 'fisher', 'hypergeometric' or 'randomization' (default='fisher')
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :type biotype_ref_path: str or pathlib.Path (default='predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :param propagate_annotations: determines the propagation method of GO annotations. \
        'no' does not propagate annotations at all; 'classic' propagates all annotations up to the DAG tree's root; \
        'elim' terminates propagation at nodes which show significant enrichment; 'weight' performs propagation in a \
        weighted manner based on the significance of children nodes relatively to their parents; and 'allm' uses a \
        combination of all proopagation methods. To read more about the propagation methods, \
        see Alexa et al: https://pubmed.ncbi.nlm.nih.gov/16606683/
        :type propagate_annotations: 'classic', 'elim', 'weight', 'all.m', or 'no' (default='elim')
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
        :type excluded_qualifiers: str, Iterable of str, or None (default='not')
        :param exclude_unannotated_genes: if True, genes that have no annotation associated with them will be \
        excluded from the enrichment analysis. This is the recommended practice for enrichment analysis, since \
        keeping unannotated genes in the analysis increases the chance of discovering spurious enrichment results.
        :type exclude_unannotated_genes: bool (deafult=True)
        :param return_nonsignificant: if True, the results DataFrame will include all tested GO terms - \
        both significant and non-significant terms. If False (default), only significant GO terms will be returned.
        :type return_nonsignificant: bool (default=False)
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default=False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_ontology_graph: bool (default=True)
        :param plot_ontology_graph: if True, will generate an ontology graph depicting the \
        significant GO terms and their parent nodes.
        :type ontology_graph_format: 'pdf', 'png', 'svg', or 'none' (default='none')
        :param ontology_graph_format: if ontology_graph_format is not 'none', the ontology graph will additonally be \
        generated in the specified file format.
        :type plot_horizontal: bool (default=True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. \
        Otherwise, results will be plotted with a vertical plot.
        :param show_expected: if True, the observed/expected values will be shown on the plot.
        :type show_expected: bool (default=False)
        :param plot_style: style for the plot. Either 'bar' for a bar plot or 'lollipop' for a lollipop plot \
        in which the lollipop size indicates the size of the observed gene set.
        :type plot_style: 'bar' or 'lollipop' (default='bar')
        :type random_seed: non-negative integer (default=None)
        :type random_seed: if using a randomization test, determine the random seed used to initialize \
        the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs. If not using a randomization test, this parameter will not affect the analysis.
        :param randomization_reps: if using a randomization test, determine how many randomization repititions to run. \
        Otherwise, this parameter will not affect the analysis.
        :type randomization_reps: int larger than 0 (default=10000)
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with GO terms as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure:: /figures/ontology_graph.png
           :align:   center
           :scale: 60 %

           Example plot of go_enrichment(plot_ontology_graph=True)


        .. figure:: /figures/plot_enrichment_results_go.png
           :align:   center
           :scale: 60 %

           Example plot of go_enrichment()


        .. figure:: /figures/plot_enrichment_results_go_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of go_enrichment(plot_horizontal = False)
        """
        propagate_annotations = propagate_annotations.lower()
        try:
            background_genes = background_genes.gene_set
        except AttributeError:
            pass
        if statistical_test.lower() == 'randomization':
            kwargs = dict(reps=randomization_reps, random_seed=random_seed)
        else:
            kwargs = {}
        runner = enrichment_runner.GOEnrichmentRunner(self.gene_set, organism, gene_id_type, alpha,
                                                      propagate_annotations, aspects, evidence_types,
                                                      excluded_evidence_types, databases, excluded_databases,
                                                      qualifiers, excluded_qualifiers, return_nonsignificant, save_csv,
                                                      fname, return_fig, plot_horizontal, plot_ontology_graph,
                                                      self.set_name, parallel_backend, statistical_test, biotype,
                                                      background_genes, biotype_ref_path, exclude_unannotated_genes,
                                                      ontology_graph_format=ontology_graph_format,
                                                      plot_style=plot_style, show_expected=show_expected, **kwargs)

        if gui_mode:
            return runner.run(plot=False), runner
        return runner.run()

    @readable_name("Enrichment for KEGG Pathways")
    def kegg_enrichment(self, organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                        gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                        alpha: param_typing.Fraction = 0.05,
                        statistical_test: Literal['fisher', 'hypergeometric', 'randomization'] = 'fisher',
                        biotype: Union[str, List[str], Literal['all']] = 'all',
                        background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                        biotype_ref_path: Union[str, Path, Literal['predefined']] = 'predefined',
                        exclude_unannotated_genes: bool = True, return_nonsignificant: bool = False,
                        save_csv: bool = False, fname=None, return_fig: bool = False, plot_horizontal: bool = True,
                        show_expected: bool = False, plot_style: Literal['bar', 'lollipop'] = 'bar',
                        plot_pathway_graphs: bool = True,
                        pathway_graphs_format: Literal[param_typing.GRAPHVIZ_FORMATS] = 'none',
                        randomization_reps: PositiveInt = 10000, random_seed: Union[int, None] = None,
                        parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                        gui_mode: bool = False
                        ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        """
        Calculates enrichment and depletion of the FeatureSet for Kyoto Encyclopedia of Genes and Genomes (KEGG) \
        curated pathways against a background set. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment score = 0) \
        appears with the smallest value in the scale.

        :param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, RNAlysis will attempt to map \
        the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :type alpha: float between 0 and 1 (default=0.05)
        :param alpha: Indicates the FDR threshold for significance.
        :param statistical_test: determines the statistical test to be used for enrichment analysis. \
        Note that some propagation methods support only some of the available statistical tests.
        :type statistical_test: 'fisher', 'hypergeometric' or 'randomization' (default='fisher')
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :type biotype_ref_path: str or pathlib.Path (default='predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :param exclude_unannotated_genes: if True, genes that have no annotation associated with them will be \
        excluded from the enrichment analysis. This is the recommended practice for enrichment analysis, since \
        keeping unannotated genes in the analysis increases the chance of discovering spurious enrichment results.
        :type exclude_unannotated_genes: bool (deafult=True)
        :param return_nonsignificant: if True, the results DataFrame will include all tested pathways - \
        both significant and non-significant ones. If False (default), only significant pathways will be returned.
        :type return_nonsignificant: bool (default=False)
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default=False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default=True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. \
        Otherwise, results will be plotted with a vertical plot.
        :param show_expected: if True, the observed/expected values will be shown on the plot.
        :type show_expected: bool (default=False)
        :param plot_style: style for the plot. Either 'bar' for a bar plot or 'lollipop' for a lollipop plot \
        in which the lollipop size indicates the size of the observed gene set.
        :type plot_style: 'bar' or 'lollipop' (default='bar')
        :type plot_pathway_graphs: bool (default=True)
        :param plot_pathway_graphs: if True, will generate pathway graphs depicting the significant KEGG pathways.
        :type pathway_graphs_format: 'pdf', 'png', 'svg', or None (default=None)
        :param pathway_graphs_format: if pathway_graphs_format is not 'none', the pathway graphs will additonally be \
        generated in the specified file format.
        :type random_seed: non-negative integer (default=None)
        :type random_seed: if using a randomization test, determine the random seed used to initialize \
        the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs. If not using a randomization test, this parameter will not affect the analysis.
        :param randomization_reps: if using a randomization test, determine how many randomization repititions to run. \
        Otherwise, this parameter will not affect the analysis.
        :type randomization_reps: int larger than 0 (default=10000)
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated pathway names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.


        .. figure:: /figures/plot_enrichment_results_kegg.png
           :align:   center
           :scale: 60 %

           Example plot of kegg_enrichment()


        .. figure:: /figures/plot_enrichment_results_kegg_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of kegg_enrichment(plot_horizontal = False)
        """
        try:
            background_genes = background_genes.gene_set
        except AttributeError:
            pass
        if statistical_test.lower() == 'randomization':
            kwargs = dict(reps=randomization_reps, random_seed=random_seed)
        else:
            kwargs = {}
        runner = enrichment_runner.KEGGEnrichmentRunner(self.gene_set, organism, gene_id_type, alpha,
                                                        return_nonsignificant, save_csv, fname, return_fig,
                                                        plot_horizontal, plot_pathway_graphs, self.set_name,
                                                        parallel_backend, statistical_test, biotype, background_genes,
                                                        biotype_ref_path, exclude_unannotated_genes,
                                                        pathway_graphs_format=pathway_graphs_format,
                                                        plot_style=plot_style, show_expected=show_expected, **kwargs)

        if gui_mode:
            return runner.run(plot=False), runner
        return runner.run()

    @readable_name("Enrichment for user-defined attributes")
    def user_defined_enrichment(self, attributes: Union[List[str], str, List[int], int, Literal['all']],
                                statistical_test: Literal['fisher', 'hypergeometric', 'randomization'] = 'fisher',
                                alpha: param_typing.Fraction = 0.05,
                                biotype: Union[str, List[str], Literal['all']] = 'all',
                                background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                                attr_ref_path: Union[str, Path, Literal['predefined']] = 'predefined',
                                biotype_ref_path: Union[str, Path, Literal['predefined']] = 'predefined',
                                exclude_unannotated_genes: bool = True, return_nonsignificant: bool = True,
                                save_csv: bool = False, fname=None, return_fig: bool = False,
                                plot_horizontal: bool = True,
                                show_expected: bool = False, plot_style: Literal['bar', 'lollipop'] = 'bar',
                                randomization_reps: PositiveInt = 10000,
                                random_seed: Union[int, None] = None,
                                parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                                gui_mode: bool = False) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        """
        Calculates enrichment and depletion of the FeatureSet for user-defined attributes against a background set.\
        The attributes are drawn from an Attribute Reference Table. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment score = 0) \
        appears with the smallest value in the scale.

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :param statistical_test: determines the statistical test to be used for enrichment analysis. \
        Note that some propagation methods support only some of the available statistical tests.
        :type statistical_test: 'fisher', 'hypergeometric' or 'randomization' (default='fisher')
        :type alpha: float between 0 and 1 (default=0.05)
        :param alpha: Indicates the FDR threshold for significance.
        :type attr_ref_path: str or pathlib.Path (default='predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default='predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all' \
        (default='protein_coding')
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object \
        (default=None)
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :param exclude_unannotated_genes: if True, genes that have no annotation associated with them will be \
        excluded from the enrichment analysis. This is the recommended practice for enrichment analysis, since \
        keeping unannotated genes in the analysis increases the chance of discovering spurious enrichment results.
        :type exclude_unannotated_genes: bool (deafult=True)
        :param return_nonsignificant: if True (default), the results DataFrame will include all tested attributes - \
        both significant and non-significant ones. If False, only significant attributes will be returned.
        :type return_nonsignificant: bool (default=True)
        :type save_csv: bool (default=False)
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path (default=None)
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default=False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default=True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :param show_expected: if True, the observed/expected values will be shown on the plot.
        :type show_expected: bool (default=False)
        :param plot_style: style for the plot. Either 'bar' for a bar plot or 'lollipop' for a lollipop plot \
        in which the lollipop size indicates the size of the observed gene set.
        :type plot_style: 'bar' or 'lollipop' (default='bar')
        :type random_seed: non-negative integer (default=None)
        :type random_seed: if using a randomization test, determine the random seed used to initialize \
        the pseudorandom generator for the randomization test. \
        By default it is picked at random, but you can set it to a particular integer to get consistents results \
        over multiple runs. If not using a randomization test, this parameter will not affect the analysis.
        :param randomization_reps: if using a randomization test, determine how many randomization repititions to run. \
        Otherwise, this parameter will not affect the analysis.
        :type randomization_reps: int larger than 0 (default=10000)
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure:: /figures/plot_enrichment_results.png
           :align:   center
           :scale: 60 %

           Example plot of user_defined_enrichment()


        .. figure:: /figures/plot_enrichment_results_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of user_defined_enrichment(plot_horizontal = False)

        """
        try:
            background_genes = background_genes.gene_set
        except AttributeError:
            pass

        if statistical_test == 'randomization':
            kwargs = dict(reps=randomization_reps)
        else:
            kwargs = dict()
        runner = enrichment_runner.EnrichmentRunner(self.gene_set, attributes, alpha, attr_ref_path,
                                                    return_nonsignificant, save_csv, fname, return_fig, plot_horizontal,
                                                    self.set_name, parallel_backend, statistical_test, biotype,
                                                    background_genes, biotype_ref_path, exclude_unannotated_genes,
                                                    single_set=False, random_seed=random_seed, plot_style=plot_style,
                                                    show_expected=show_expected, **kwargs)
        if gui_mode:
            return runner.run(plot=False), runner
        return runner.run()

    @readable_name("Enrichment for user-defined non-categorical attributes")
    def non_categorical_enrichment(self, attributes: Union[List[str], str, List[int], int, Literal['all']] = None,
                                   alpha: param_typing.Fraction = 0.05, parametric_test: bool = False,
                                   biotype: Union[str, List[str], Literal['all']] = 'all',
                                   background_genes: Union[Set[str], Filter, 'FeatureSet'] = None,
                                   attr_ref_path: Union[str, Path, Literal['predefined']] = 'predefined',
                                   biotype_ref_path: Union[str, Path, Literal['predefined']] = 'predefined',
                                   plot_log_scale: bool = True,
                                   plot_style: Literal['interleaved', 'overlap'] = 'overlap', n_bins: PositiveInt = 50,
                                   save_csv: bool = False, fname=None, return_fig: bool = False, gui_mode: bool = False
                                   ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, List[plt.Figure]]]:
        """
        Calculates enrichment and depletion of the FeatureSet for user-defined non-categorical attributes \
        against a background set using either a one-sample T-test or Sign test. \
        The attributes are drawn from an Attribute Reference Table. \
        The background set is determined by either the input variable ‘background_genes’, \
        or by the input variable ‘biotype’ and a Biotype Reference Table. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method).

        :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'.
        :param attributes: An iterable of attribute names or attribute numbers \
        (according to their order in the Attribute Reference Table). \
        If 'all', all of the attributes in the Attribute Reference Table will be used. \
        If None, a manual input prompt will be raised.
        :type alpha: float between 0 and 1 (default=0.05)
        :param alpha: Indicates the FDR threshold for significance.
        :param parametric_test: if True, performs a parametric statistical test (one-sample t-test). \
        If False (default), performs a non-parametric statistical test (sign test).
        :type parametric_test: bool (default=False)
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all' \
        (default='protein_coding')
        :param biotype: determines the background genes by their biotype. Requires specifying a Biotype Reference Table. \
        'all' will include all genomic features in the reference table, \
        'protein_coding' will include only protein-coding genes from the reference table, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object \
        (default=None)
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :type attr_ref_path: str or pathlib.Path (default='predefined')
        :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :type biotype_ref_path: str or pathlib.Path (default='predefined')
        :param biotype_ref_path: the path of the Biotype Reference Table. \
        Will be used to generate background set if 'biotype' is specified.
        :param plot_log_scale: if True (default), the Y-axis of the enrichment plot will be logarithmic. \
        Otherwise, the Y-axis of the enrichment plot will be linear.
        :type plot_log_scale: bool (default=True)
        :param plot_style: 'interleaved' will plot an interleaved histogram. \
        'overlap' will plot a semi-transparent histogram where the obsreved and expected are overlapping.
        :type plot_style: 'overlap' or 'interleaved' (default='overlap')
        :param n_bins: the number of bins to display in the enrichment plot histograms
        :type n_bins: int larger than 0 (default=50)
        :type save_csv: bool (default=False)
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path (default=None)
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default=False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.


        .. figure:: /figures/hist_overlap.png
           :align:   center
           :scale: 60 %

           Example plot of non_categorical_enrichment(plot_style`='overlap')

        .. figure:: /figures/hist_interleaved.png
           :align:   center
           :scale: 60 %

           Example plot of non_categorical_enrichment(plot_style='interleaved')
        """

        try:
            background_genes = background_genes.gene_set
        except AttributeError:
            pass
        runner = enrichment_runner.NonCategoricalEnrichmentRunner(self.gene_set, attributes, alpha, biotype,
                                                                  background_genes, attr_ref_path,
                                                                  biotype_ref_path, save_csv, fname,
                                                                  return_fig, plot_log_scale, plot_style,
                                                                  n_bins, self.set_name, parallel_backend='sequential',
                                                                  parametric_test=parametric_test)
        if gui_mode:
            return runner.run(plot=False), runner
        return runner.run()

    @readable_name('Summarize feature biotypes (based on a reference table)')
    def biotypes_from_ref_table(self, ref: Union[str, Path, Literal['predefined']] = 'predefined'):
        """
        Returns a DataFrame of the biotypes in the gene set and their count.

        :type ref: str or pathlib.Path (default='predefined')
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

        ref = settings.get_biotype_ref_path(ref)
        ref_df = io.load_table(ref)
        validation.validate_biotype_table(ref_df)
        ref_df.columns = ref_df.columns.str.lower()
        not_in_ref = pd.Index(self.gene_set).difference(set(ref_df['gene']))
        if len(not_in_ref) > 0:
            warnings.warn(
                f'{len(not_in_ref)} of the features in the Filter object do not appear in the Biotype Reference Table. ')
            ref_df = pd.concat(
                [ref_df, pd.DataFrame({'gene': not_in_ref, 'biotype': '_missing_from_biotype_reference'})])
        return ref_df.set_index('gene', drop=False).loc[parsing.data_to_list(self.gene_set)].groupby('biotype').count()


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
        assert len(self.ranked_genes) == len(self.gene_set), "'ranked_genes' must have no repeating elements!"

    def __copy__(self):
        obj = type(self)(self.ranked_genes.copy(), self.set_name)
        return obj

    def __iter__(self):
        return self.ranked_genes.__iter__()

    def __eq__(self, other):
        if type(self) != type(other):
            return False

        if self.set_name != other.set_name:
            return False

        if len(self.ranked_genes) != len(other.ranked_genes) or not np.all(self.ranked_genes == other.ranked_genes):
            return False

        return True

    def _convert_to_filter_obj(self) -> Filter:
        return Filter.from_dataframe(pd.DataFrame(index=self.ranked_genes), self.set_name)

    def _inplace(self, func, func_kwargs, inplace: bool, **update_kwargs):
        """
        Executes the user's choice whether to filter in-place or create a new instance of the FeatureSet object.
        """
        filter_obj = self._convert_to_filter_obj()
        applied = func(filter_obj, **func_kwargs, inplace=inplace)
        if inplace:
            applied = filter_obj

        new_set = applied.index_set
        suffix = '_' + applied.fname.stem.split('_')[-1]
        new_name = self.set_name + suffix
        new_ranked = []
        if len(new_set.difference(self.ranked_genes)) > 0:
            new_ranked = applied.df.index
        else:
            for item in self.ranked_genes:
                if item in new_set:
                    new_ranked.append(item)

        new_ranked = np.array(new_ranked, dtype='str')
        # if inplace, modify self, name and other properties of self
        if inplace:
            self.intersection_update(new_set)
            self._update(ranked_genes=new_ranked, gene_set=new_set, set_name=new_name, **update_kwargs)
        # if not inplace, copy self, modify the self, name, and other properties of the copy, and return it
        else:
            new_obj = type(self)(new_ranked, new_name)
            new_obj._update(**update_kwargs)
            return new_obj

    def _set_ops(self, others: Union[set, 'FeatureSet'], op: types.FunctionType):
        warnings.warn("Warning: when performing set operations with RankedSet objects, "
                      "the return type will always be FeatureSet and not RankedSet.")
        return super()._set_ops(others, op)

    @readable_name("Single-Set GO Enrichment")
    def single_set_go_enrichment(self, organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                                 gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                 alpha: param_typing.Fraction = 0.05,
                                 propagate_annotations: Literal['classic', 'elim', 'weight', 'all.m', 'no'] = 'elim',
                                 aspects: Union[Literal[('any',) + GO_ASPECTS], Iterable[Literal[GO_ASPECTS]]] = 'any',
                                 evidence_types: Union[
                                     Literal[('any',) + GO_EVIDENCE_TYPES], Iterable[Literal[
                                         GO_EVIDENCE_TYPES]]] = 'any',
                                 excluded_evidence_types: Union[
                                     Literal[GO_EVIDENCE_TYPES], Iterable[Literal[GO_EVIDENCE_TYPES]]] = (),
                                 databases: Union[str, Iterable[str], Literal['any']] = 'any',
                                 excluded_databases: Union[str, Iterable[str]] = (),
                                 qualifiers: Union[
                                     Literal[('any',) + GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'any',
                                 excluded_qualifiers: Union[
                                     Literal[GO_QUALIFIERS], Iterable[Literal[GO_QUALIFIERS]]] = 'not',
                                 exclude_unannotated_genes: bool = True, return_nonsignificant: bool = False,
                                 save_csv: bool = False, fname=None,
                                 return_fig: bool = False, plot_horizontal: bool = True,
                                 show_expected: bool = False, plot_style: Literal['bar', 'lollipop'] = 'bar',
                                 plot_ontology_graph: bool = True,
                                 ontology_graph_format: Literal[param_typing.GRAPHVIZ_FORMATS] = 'none',
                                 parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                                 gui_mode: bool = False
                                 ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        """
        Calculates enrichment and depletion of the sorted RankedSet for Gene Ontology (GO) terms \
        WITHOUT a background set, using the generalized Minimum Hypergeometric Test (XL-mHG, developed by  \
        `Prof. Zohar Yakhini and colleagues <https://dx.doi.org/10.1371/journal.pcbi.0030039/>`_ \
        and generalized by \
        `Florian Wagner <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143196/>`_). \
        The GO terms and annotations are drawn via the GO Solr search engine GOlr, \
        using the search terms defined by the user. \
        P-values are calculated using the generalized Minimum Hypergeometric Test. \
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
        :type gene_id_type: str or 'auto' (default='auto')
        :type alpha: float between 0 and 1
        :param alpha: Indicates the FDR threshold for significance.
        :param propagate_annotations: determines the propagation method of GO annotations. \
        'no' does not propagate annotations at all; 'classic' propagates all annotations up to the DAG tree's root; \
        'elim' terminates propagation at nodes which show significant enrichment; 'weight' performs propagation in a \
        weighted manner based on the significance of children nodes relatively to their parents; and 'allm' uses a \
        combination of all proopagation methods. To read more about the propagation methods, \
        see Alexa et al: https://pubmed.ncbi.nlm.nih.gov/16606683/
        :type propagate_annotations: 'classic', 'elim', 'weight', 'all.m', or 'no' (default='elim')
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
        :type excluded_qualifiers: str, Iterable of str, or None (default='not')
        :param exclude_unannotated_genes: if True, genes that have no annotation associated with them will be \
        excluded from the enrichment analysis. This is the recommended practice for enrichment analysis, since \
        keeping unannotated genes in the analysis increases the chance of discovering spurious enrichment results.
        :type exclude_unannotated_genes: bool (deafult=True)
        :param return_nonsignificant: if True, the results DataFrame will include all tested GO terms - \
        both significant and non-significant terms. If False (default), \
        only significant GO terms will be returned.
        :type return_nonsignificant: bool (default=False)
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), \
        fname will be requested in a manual prompt.
        :type return_fig: bool (default=False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default=True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. \
        Otherwise, results will be plotted with a vertical plot.
        :param show_expected: if True, the observed/expected values will be shown on the plot.
        :type show_expected: bool (default=False)
        :param plot_style: style for the plot. Either 'bar' for a bar plot or 'lollipop' for a lollipop plot \
        in which the lollipop size indicates the size of the observed gene set.
        :type plot_style: 'bar' or 'lollipop' (default='bar')
        :type plot_ontology_graph: bool (default=True)
        :param plot_ontology_graph: if True, will generate an ontology graph depicting the \
        significant GO terms and their parent nodes.
        :type ontology_graph_format: 'pdf', 'png', 'svg', or None (default=None)
        :param ontology_graph_format: if ontology_graph_format is not 'none', the ontology graph will additonally be \
        generated in the specified file format.
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure:: /figures/ontology_graph_singlelist.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_go_enrichment(plot_ontology_graph=True)


        .. figure:: /figures/plot_enrichment_results_go_singlelist.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_go_enrichment()


        .. figure:: /figures/plot_enrichment_results_go_singlelist_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_go_enrichment(plot_horizontal = False)
        """
        runner = enrichment_runner.GOEnrichmentRunner(self.ranked_genes, organism, gene_id_type, alpha,
                                                      propagate_annotations, aspects, evidence_types,
                                                      excluded_evidence_types, databases, excluded_databases,
                                                      qualifiers, excluded_qualifiers, return_nonsignificant, save_csv,
                                                      fname, return_fig, plot_horizontal, plot_ontology_graph,
                                                      self.set_name,
                                                      parallel_backend=parallel_backend, enrichment_func_name='xlmhg',
                                                      exclude_unannotated_genes=exclude_unannotated_genes,
                                                      single_set=True, ontology_graph_format=ontology_graph_format,
                                                      plot_style=plot_style, show_expected=show_expected)

        if gui_mode:
            return runner.run(plot=False), runner
        return runner.run()

    @readable_name("Single-Set Enrichment for KEGG Pathways")
    def single_set_kegg_enrichment(self,
                                   organism: Union[str, int, Literal['auto'], Literal[DEFAULT_ORGANISMS]] = 'auto',
                                   gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                                   alpha: param_typing.Fraction = 0.05, exclude_unannotated_genes: bool = True,
                                   return_nonsignificant: bool = False, save_csv: bool = False,
                                   fname=None, return_fig: bool = False, plot_horizontal: bool = True,
                                   show_expected: bool = False, plot_style: Literal['bar', 'lollipop'] = 'bar',
                                   plot_pathway_graphs: bool = True,
                                   pathway_graphs_format: Literal[param_typing.GRAPHVIZ_FORMATS] = 'none',
                                   parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky', gui_mode: bool = False
                                   ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, plt.Figure]]:
        """
        Calculates enrichment and depletion of the sorted RankedSet for Kyoto Encyclopedia of Genes and Genomes (KEGG) \
        curated pathways WITHOUT a background set, using the generalized Minimum Hypergeometric Test \
        (XL-mHG, developed by `Prof. Zohar Yakhini and colleagues <https://dx.doi.org/10.1371/journal.pcbi.0030039/>`_ \
        and generalized by \
        `Florian Wagner <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143196/>`_). \
        P-values are calculated using the generalized Minimum Hypergeometric Test. \
        P-values are corrected for multiple comparisons using the Benjamini–Hochberg step-up procedure \
        (original FDR method). In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :param organism: organism name or NCBI taxon ID for which the function will fetch GO annotations.
        :type organism: str or int
        :param gene_id_type: the identifier type of the genes/features in the FeatureSet object \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
        If the annotations fetched from the KEGG server do not match your gene_id_type, \
        RNAlysis will attempt to map  the annotations' gene IDs to your identifier type. \
        For a full list of legal 'gene_id_type' names, see the UniProt website: \
        https://www.uniprot.org/help/api_idmapping
        :type gene_id_type: str or 'auto' (default='auto')
        :type alpha: float between 0 and 1
        :param alpha: Indicates the FDR threshold for significance.
        :param exclude_unannotated_genes: if True, genes that have no annotation associated with them will be \
        excluded from the enrichment analysis. This is the recommended practice for enrichment analysis, since \
        keeping unannotated genes in the analysis increases the chance of discovering spurious enrichment results.
        :type exclude_unannotated_genes: bool (deafult=True)
        :param return_nonsignificant: if True, the results DataFrame will include all tested GO terms - \
        both significant and non-significant terms. If False (default), \
        only significant KEGG pathways will be returned.
        :type return_nonsignificant: bool (default=False)
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), \
        fname will be requested in a manual prompt.
        :type return_fig: bool (default=False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default=True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. \
        Otherwise, results will be plotted with a vertical plot.
        :param show_expected: if True, the observed/expected values will be shown on the plot.
        :type show_expected: bool (default=False)
        :param plot_style: style for the plot. Either 'bar' for a bar plot or 'lollipop' for a lollipop plot \
        in which the lollipop size indicates the size of the observed gene set.
        :type plot_style: 'bar' or 'lollipop' (default='bar')
        :type plot_pathway_graphs: bool (default=True)
        :param plot_pathway_graphs: if True, will generate pathway graphs depicting the significant KEGG pathways.
        :type pathway_graphs_format: 'pdf', 'png', 'svg', or None (default=None)
        :param pathway_graphs_format: if pathway_graphs_format is not 'none', the pathway graphs will additonally be \
        generated in the specified file format.
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure:: /figures/pathway_graph_singlelist.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_kegg_enrichment(plot_pathway_graphs=True)


        .. figure:: /figures/plot_enrichment_results_kegg_single_set.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_kegg_enrichment()


        .. figure:: /figures/plot_enrichment_results_kegg_single_set_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_kegg_enrichment(plot_horizontal = False)
        """
        runner = enrichment_runner.KEGGEnrichmentRunner(self.ranked_genes, organism, gene_id_type, alpha,
                                                        return_nonsignificant, save_csv, fname, return_fig,
                                                        plot_horizontal, plot_pathway_graphs, self.set_name,
                                                        parallel_backend=parallel_backend, enrichment_func_name='xlmhg',
                                                        exclude_unannotated_genes=exclude_unannotated_genes,
                                                        single_set=True, pathway_graphs_format=pathway_graphs_format,
                                                        plot_style=plot_style, show_expected=show_expected)

        if gui_mode:
            return runner.run(plot=False), runner
        return runner.run()

    @readable_name("Single-Set Enrichment for user-defined attributes")
    def single_set_enrichment(self, attributes: Union[List[str], str, List[int], int, Literal['all']],
                              alpha: param_typing.Fraction = 0.05,
                              attr_ref_path: Union[str, Path, Literal['predefined']] = 'predefined',
                              exclude_unannotated_genes: bool = True, return_nonsignificant: bool = True,
                              save_csv: bool = False, fname=None, return_fig: bool = False,
                              plot_horizontal: bool = True,
                              show_expected: bool = False, plot_style: Literal['bar', 'lollipop'] = 'bar',
                              parallel_backend: Literal[PARALLEL_BACKENDS] = 'loky',
                              gui_mode: bool = False):
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
        :type attr_ref_path: str or pathlib.Path (default='predefined')
        :param attr_ref_path: path of the Attribute Reference Table from which user-defined attributes will be drawn.
        :param exclude_unannotated_genes: if True, genes that have no annotation associated with them will be \
        excluded from the enrichment analysis. This is the recommended practice for enrichment analysis, since \
        keeping unannotated genes in the analysis increases the chance of discovering spurious enrichment results.
        :type exclude_unannotated_genes: bool (default=True)
        :param return_nonsignificant: if True (default), the results DataFrame will include all tested attributes - \
        both significant and non-significant ones. If False, only significant attributes will be returned.
        :type return_nonsignificant: bool (default=True)
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        'C:/dir/file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :type return_fig: bool (default=False)
        :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :type plot_horizontal: bool (default=True)
        :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. Otherwise, results \
        will be plotted with a vertical plot.
        :type parallel_backend: Literal[PARALLEL_BACKENDS] (default='loky')
        :param parallel_backend: Determines the babckend used to run the analysis. \
        if parallel_backend not 'sequential', will calculate the statistical tests using parallel processing. \
        In most cases parallel processing will lead to shorter computation time, but does not affect the results of \
        the analysis otherwise.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index; \
        and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure:: /figures/plot_enrichment_results_single_set.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_enrichment()


        .. figure:: /figures/plot_enrichment_results_single_set_vertical.png
           :align:   center
           :scale: 60 %

           Example plot of single_set_enrichment(plot_horizontal = False)

        """
        runner = enrichment_runner.EnrichmentRunner(self.ranked_genes, attributes, alpha, attr_ref_path,
                                                    return_nonsignificant, save_csv, fname, return_fig, plot_horizontal,
                                                    self.set_name, parallel_backend=parallel_backend,
                                                    enrichment_func_name='xlmhg',
                                                    exclude_unannotated_genes=exclude_unannotated_genes,
                                                    single_set=True, plot_style=plot_style, show_expected=show_expected)
        if gui_mode:
            return runner.run(plot=False), runner
        return runner.run()


def enrichment_bar_plot(results_table_path: Union[str, Path], alpha: param_typing.Fraction = 0.05,
                        enrichment_score_column: Union[
                            str, Literal['log2_fold_enrichment', 'log2_enrichment_score']] = 'log2_fold_enrichment',
                        n_bars: Union[param_typing.PositiveInt, Literal['all']] = 'all',
                        title: str = 'Enrichment results', center_bars: bool = True,
                        plot_horizontal: bool = True, ylabel: Union[str, Literal[
        r"$\log_2$(Fold Enrichment)", r"$\log_2$(Enrichment Score)"]] = r"$\log_2$(Fold Enrichment)",
                        ylim: Union[float, Literal['auto']] = 'auto',
                        plot_style: Literal['bar', 'lollipop'] = 'bar', show_expected: bool = False) -> plt.Figure:
    """
    Generate an enrichment bar-plot based on an enrichment results table. \
    For the clarity of display, complete depletion (linear enrichment = 0) \
    appears with the smallest value in the scale.

    :param results_table_path: Path to the results table returned by enrichment functions.
    :type results_table_path: str or Path
    :param alpha: the threshold for statistical significance. Used to draw significance asterisks on the graph.
    :type alpha: float between 0 and 1 (default=0.05)
    :param enrichment_score_column: name of the table column containing enrichment scores.
    :type enrichment_score_column: str (default='log2_fold_enrichment')
    :param n_bars: number of bars to display in the bar plot. If n_bars='all', \
     all the results will be displayed on the graph. Otherwise, only the top n results will be displayed on the graph.
    :type n_bars: int > 1 or 'all' (default='all')
    :param title: plot title.
    :type title: str
    :param plot_horizontal: if True, results will be plotted with a horizontal bar plot. \
        Otherwise, results will be plotted with a vertical plot.
    :type plot_horizontal: bool (default=True)
    :param ylabel: plot y-axis label.
    :type ylabel: str (default=r"$\log_2$(Fold Enrichment)")
    :param center_bars: if True, center the bars around Y=0. Otherwise, ylim is determined by min/max values.
    :type center_bars: bool (default=True)
    :param ylim: set the Y-axis limits. If `ylim`='auto', determines the axis limits automatically based on the data. \
    If `ylim` is a number, set the Y-axis limits to [-ylim, ylim].
    :type ylim: float or 'auto' (default='auto')
    :param plot_style: style for the plot. Either 'bar' for a bar plot or 'lollipop' for a lollipop plot \
    in which the lollipop size indicates the size of the observed gene set.
    :type plot_style: 'bar' or 'lollipop' (default='bar')
    :param show_expected: if True, the observed/expected values will be shown on the plot.
    :type show_expected: bool (default=False)
    :return: Figure object containing the bar plot
    :rtype: matplotlib.figure.Figure instance
    """
    results_table = io.load_table(results_table_path, 0)
    runner = enrichment_runner.EnrichmentRunner(set(), results_table.index, alpha, '', True, False, '', True,
                                                plot_horizontal, '', False, 'hypergeometric', 'all',
                                                plot_style=plot_style, show_expected=show_expected)
    runner.en_score_col = enrichment_score_column
    runner.results = results_table
    return runner.enrichment_bar_plot(n_bars, center_bars, ylabel=ylabel, title=title, ylim=ylim)


def _fetch_sets(objs: dict, ref: Union[str, Path, Literal['predefined']] = 'predefined'):
    """
    Receives the 'objs' input from enrichment.upset_plot() and enrichment.venn_diagram(), and turns the values in it \
    into python sets.

    :param objs: the 'objs' input given to the function enrichment.upset_plot() or enrichment.venn_diagram().
    :type objs: a dictionary, where the keys are names of sets, and the values are either\
     python sets, FeatureSets or names of columns in the Attribute Reference Table.
    :param ref: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type ref: str or pathlib.Path (default='predefined')
    :return: a dictionary, where the keys are names of sets and the values are python sets of feature indices.
    """
    assert isinstance(objs, dict), f"objs must be a dictionary. Instaed got {type(objs)}"
    fetched_sets = dict()

    if validation.isinstanceiter_any(objs.values(), str):
        attr_ref_table = io.load_table(settings.get_attr_ref_path(ref))
        validation.validate_attr_table(attr_ref_table)
        attr_ref_table.set_index('gene', inplace=True)

    for set_name, set_obj in zip(objs.keys(), objs.values()):
        if validation.isinstanceinh(set_obj, FeatureSet):
            fetched_sets[set_name] = set_obj.gene_set
        elif isinstance(objs[set_name], (set, list, tuple)):
            fetched_sets[set_name] = parsing.data_to_set(set_obj)
        elif validation.isinstanceinh(set_obj, Filter):
            fetched_sets[set_name] = set_obj.index_set
        elif isinstance(set_obj, str):
            fetched_sets[set_name] = set(attr_ref_table[set_obj].loc[attr_ref_table[set_obj].notna()].index)
        else:
            raise TypeError(f"Invalid type for the set '{set_name}': {set_obj}.")

    return fetched_sets


def upset_plot(objs: Dict[str, Union[str, FeatureSet, Set[str]]], set_colors: param_typing.ColorList = ('black',),
               title: str = 'UpSet Plot', title_fontsize: float = 20, show_percentages: bool = True,
               attr_ref_table_path: Union[str, Path, Literal['predefined']] = 'predefined', fig: plt.Figure = None
               ) -> plt.Figure:
    """
    Generate an UpSet plot of 2 or more sets, FeatureSets or attributes from the Attribute Reference Table.

    :param objs: the FeatureSets, python sets or user-defined attributes to plot.
    :type objs: a dictionary with 2 or more entries, where the keys are the names of the sets, and the values are either \
    a FeatureSet, a python set of feature indices, or a name of a column in the Attribute Reference Table. \
    For example: \
    {'first set':{'gene1','gene2','gene3'}, 'second set':'name_of_attribute_from_reference_table'}
    :param set_colors: I\if one color is supplied, this will determine the color of all sets on the plot. \
    If multiple colors are supplied, this will determine the color of each set on the plot, and the subset colors \
    will be determined by mixing.
    :type set_colors: Iterable of colors (default=('black',)
    :param title: determines the title of the plot.
    :type title: str
    :param title_fontsize: font size for the plot's title
    :type title_fontsize: float (default=20)
    :param show_percentages: if True, shows the percentage that each set or subset takes out of the entire dataset.
    :type show_percentages: bool (default=True)
    :param attr_ref_table_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type attr_ref_table_path: str or pathlib.Path (default='predefined')
    :param fig: optionally, supply your own Figure to generate the plot onto.
    :type fig: matplotlib.Figure
    :returns: plt.Figure


        .. figure:: /figures/upsetplot.png
           :align:   center
           :scale: 70 %

           Example plot of upset_plot()
    """

    upset_df = parsing.generate_upset_series(_fetch_sets(objs=objs, ref=attr_ref_table_path))
    with pd.option_context("mode.copy_on_write", False):
        upset_obj = upsetplot.UpSet(upset_df, sort_by='degree', sort_categories_by=None,
                                    show_percentages=show_percentages)
        axes = upset_obj.plot(fig=fig)
        if fig is None:
            fig = plt.gcf()

        set_colors = [matplotlib.colors.to_rgb(color) for color in parsing.data_to_list(set_colors)]
        if len(set_colors) == 1:
            set_colors = [set_colors[0]] * len(objs)
        elif len(set_colors) > len(objs):
            set_colors = set_colors[0:len(objs)]
        elif len(set_colors) < len(objs):
            set_colors = set_colors + [(0, 0, 0)] * (len(objs) - len(set_colors))

        for main_set in range(len(objs)):
            color = set_colors[main_set]
            axes['totals'].patches[main_set].set_facecolor(color)

        patcn_ids = _get_tuple_patch_ids(len(objs))
        for subset in range(len(upset_obj.subset_styles)):
            id = patcn_ids[subset]
            colors_to_mix = [set_colors[ind] for ind, is_present in enumerate(id) if is_present]
            color = generic.mix_colors(*colors_to_mix)
            upset_obj.subset_styles[subset]['facecolor'] = color
            axes['intersections'].patches[subset].set_facecolor(color)
        matrix_ax = axes['matrix']
        matrix_ax.clear()
        upset_obj.plot_matrix(matrix_ax)
        fig.suptitle(title, fontsize=title_fontsize)

        plt.show()
        return fig


def _compare_ids(id1: Tuple[int, ...], id2: Tuple[int, ...]):
    if sum(id1) > sum(id2):
        return 1
    elif sum(id1) < sum(id2):
        return -1
    else:
        for i in range(len(id1), 0, -1):
            ind = i - 1
            if id1[ind] > id2[ind]:
                return 1
            elif id1[ind] < id2[ind]:
                return -1
        return 0


def _get_tuple_patch_ids(n_sets: int) -> List[Tuple[int, ...]]:
    unsorted_ids = list(itertools.product([0, 1], repeat=n_sets))[1:]
    sorted_ids = sorted(unsorted_ids, key=functools.cmp_to_key(_compare_ids))
    return sorted_ids


def venn_diagram(objs: Dict[str, Union[str, FeatureSet, Set[str]]], title: Union[str, Literal['default']] = 'default',
                 attr_ref_table_path: Union[str, Path, Literal['predefined']] = 'predefined',
                 set_colors: param_typing.ColorList = ('r', 'g', 'b'),
                 transparency: param_typing.Fraction = 0.4, weighted: bool = True, add_outline: bool = True,
                 linecolor: param_typing.Color = 'black', linestyle: Literal['solid', 'dashed'] = 'solid',
                 linewidth: float = 2.0, title_fontsize: float = 14, set_fontsize: float = 12,
                 subset_fontsize: float = 10, normalize_to: float = 1.0, fig: plt.Figure = None) -> plt.Figure:
    """
    Generate a Venn diagram of 2 to 3 sets, FeatureSets or attributes from the Attribute Reference Table.

    :param objs: the FeatureSets, python sets or user-defined attributes to plot.
    :type objs: a dictionary with 2-3 entries, where the keys are the names of the sets, and the values are either \
    a FeatureSet, a python set of feature indices, or a name of a column in the Attribute Reference Table. \
    For example: \
    {'first set':{'gene1','gene2','gene3'}, 'second set':'name_of_attribute_from_reference_table'}
    :type title: str or 'default' (default='default')
    :param set_colors: determines the colors of the circles in the diagram.
    :param attr_ref_table_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn, \
    if such attributes are included in 'objs'.
    :type attr_ref_table_path: str or pathlib.Path (default='predefined')
    :param title: determines the title of the plot.
    :type set_colors: tuple of matplotlib-format colors, the same size as 'objs'
    :param transparency: determines the opacity of the circles. \
    Opacity of 0 is completely transparent, while opacity of 1 is completely opaque.
    :type transparency: a float between 0 and 1
    :param weighted: if True, the plot will be area-weighted.
    :type weighted: bool (default=True)
    :param add_outline: if True, adds an outline to the circles.
    :type add_outline: bool (default=True)
    :param linecolor: Determines the color of the circles' outline.
    :type linecolor: matplotlib-format color (default='black')
    :param linestyle: the style of the circles' outline.
    :type linestyle: 'solid' or 'dashed' (default='solid')
    :param linewidth: the widdth of the circles' outlines.
    :type linewidth: float (default=2.0)
    :param title_fontsize: font size for the plot's title.
    :type title_fontsize: float (default=14)
    :param set_fontsize: font size for the set labels.
    :type set_fontsize: float (default=12)
    :param subset_fontsize: font size for the subset labels.
    :type subset_fontsize: float (default=10)
    :param fig: optionally, supply your own Figure to generate the plot onto.
    :type fig: matplotlib.Figure

    :param normalize_to: the total (on-axes) area of the circles to be drawn. Sometimes tuning it (together
    with the overall fiture size) may be useful to fit the text labels better.
    :type normalize_to: float (default=1.0)
    :return: a tuple of a VennDiagram object; and a list of 2-3 Circle patches.


        .. figure:: /figures/venn.png
           :align:   center
           :scale: 70 %

           Example plot of venn_diagram()
    """
    if len(objs) > 3 or len(objs) < 2:
        raise ValueError(f'Venn can only accept between 2 and 3 sets. Instead got {len(objs)}')
    assert isinstance(title, str), f'Title must be a string. Instead got {type(title)}'
    objs = _fetch_sets(objs=objs, ref=attr_ref_table_path)
    set_colors = parsing.data_to_tuple(set_colors)
    if len(set_colors) == 1:
        set_colors *= 3

    if len(objs) == 2:
        func = vn.venn2 if weighted else vn.venn2_unweighted
        func_circles = vn.venn2_circles
        set_colors = set_colors[0:2]
    else:
        func = vn.venn3 if weighted else vn.venn3_unweighted
        func_circles = vn.venn3_circles
        set_colors = set_colors[0:3]
    if fig is None:
        fig = plt.figure()
    ax = fig.add_subplot()
    plot_obj = func(tuple(objs.values()), tuple(objs.keys()), set_colors=set_colors, alpha=transparency,
                    normalize_to=normalize_to, ax=ax)
    if add_outline and weighted:
        circle_obj = func_circles(tuple(objs.values()), color=linecolor, linestyle=linestyle, linewidth=linewidth,
                                  normalize_to=normalize_to, ax=ax)
    elif add_outline and not weighted:
        circle_obj = func(tuple(objs.values()), tuple(objs.keys()), alpha=1, normalize_to=normalize_to, ax=ax)
        for patch in circle_obj.patches:
            patch.set_edgecolor(linecolor)
            patch.set_linewidth(linewidth)
            patch.set_linestyle(linestyle)
            patch.set_fill(False)
    else:
        circle_obj = None

    for label in plot_obj.set_labels:
        label.set_fontsize(set_fontsize)
    for sublabel in plot_obj.subset_labels:
        if sublabel is not None:
            sublabel.set_fontsize(subset_fontsize)

    if title == 'default':
        title = 'Venn diagram of ' + ''.join([name + ' ' for name in objs.keys()])
    ax.set_title(title, fontsize=title_fontsize)

    plt.show()
    return fig


def gene_ontology_graph(aspect: Literal[param_typing.GO_ASPECTS], results_table_path: Union[str, Path],
                        enrichment_score_column: Union[
                            str, Literal['log2_enrichment_score', 'log2_fold_enrichment']] = 'log2_fold_enrichment',
                        title: Union[str, Literal['auto']] = 'auto', ylabel: str = r"$\log_2$(Fold Enrichment)",
                        graph_format: Literal[param_typing.GRAPHVIZ_FORMATS] = 'none', dpi: PositiveInt = 300
                        ) -> Union[plt.Figure, None]:
    """
    Generate a GO enrichment ontology graph based on an enrichment results table.

    :param aspect: The GO aspect to generate an ontology graph for.
    :type aspect: 'biological_process', 'molecular_function', or 'cellular_component'
    :param results_table_path: Path to the results table returned by enrichment functions.
    :type results_table_path: str or Path
    :param enrichment_score_column: name of the table column that contains the enrichment scores.
    :type enrichment_score_column: str (default='log2_fold_enrichment')
    :param title: plot title.
    :type title: str or 'auto' (default='auto')
    :param ylabel: plot y-axis label.
    :type ylabel: str (default=r"$\log_2$(Fold Enrichment)")
    :param graph_format: if `graph_format` is not 'none', the ontology graph will additonally be \
    generated in the specified file format.
    :type graph_format: 'pdf', 'png', 'svg', or 'none' (default='none')
    :param dpi: resolution of the ontology graph in DPI (dots per inch).
    :type dpi: int (default=300)
    """
    results_df = io.load_table(results_table_path, index_col=0)
    assert enrichment_score_column in results_df, f"Invalid enrichment_score_col '{enrichment_score_column}'"
    dag_tree = ontology.fetch_go_basic()
    return dag_tree.plot_ontology(aspect, results_df, enrichment_score_column, title, ylabel, graph_format, dpi)


def kegg_pathway_graph(pathway_id: str, marked_genes: Union[Sequence[str], None],
                       gene_id_type: Union[str, Literal['auto'], Literal[get_gene_id_types()]] = 'auto',
                       title: Union[str, Literal['auto']] = 'auto', ylabel: str = '',
                       graph_format: Literal[param_typing.GRAPHVIZ_FORMATS] = 'none', dpi: PositiveInt = 300
                       ) -> Union[plt.Figure, None]:
    """
    Generate a KEGG Pathway graph.

    :param pathway_id: KEGG ID of the pathway to be plotted.
    :type pathway_id: str
    :param marked_genes: a set of genes/genomic features to be highlighted on the pathway graph. \
    The gene ID type of those genes should match the parameter `gene_id_type`.
    :type marked_genes: sequence of str or None
    :param gene_id_type: the identifier type you want to use when displaying genes in the graph \
        (for example: 'UniProtKB', 'WormBase', 'RNACentral', 'Entrez Gene ID'). \
    :type gene_id_type: str or 'auto' (default='auto')
    :param title: plot title.
    :type title: str or 'auto' (default='auto')
    :param ylabel: plot y-axis label.
    :type ylabel: str (default=r"$\log_2$(Fold Enrichment)")
    :param graph_format: if `graph_format` is not 'none', the ontology graph will additonally be \
    generated in the specified file format.
    :type graph_format: 'pdf', 'png', 'svg', or 'none' (default='none')
    :param dpi: resolution of the ontology graph in DPI (dots per inch).
    :type dpi: int (default=300)
    """
    if marked_genes is None:
        translator = None
    else:
        if gene_id_type.lower() == 'auto':
            translator, _, _ = io.find_best_gene_mapping(parsing.data_to_tuple(marked_genes), ('KEGG',), None)
        else:
            translator = io.map_gene_ids(parsing.data_to_tuple(marked_genes), 'KEGG', gene_id_type)

    pathway = ontology.fetch_kegg_pathway(pathway_id, translator)
    return pathway.plot_pathway(marked_genes, title, ylabel, graph_format, dpi)
