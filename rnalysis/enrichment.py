"""
This module can perform enrichment analyses on a given set of genomic features and visualize their intersections. \
These include gene ontology/tissue/phenotype enrichment, enrichment for user-defined attributes, \
set visualization ,etc. \
Results of enrichment analyses can be saved to .csv files.
"""
import random
import numpy as np
import pandas as pd
from rnalysis import general, filtering
import tissue_enrichment_analysis as tea
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from pathlib import Path
import statsmodels.stats.multitest as multitest
from ipyparallel import Client
from itertools import repeat, compress
import upsetplot as upset
import matplotlib_venn as vn
import warnings
from typing import Union, List, Set, Dict, Tuple, Iterable, Type, Callable


class FeatureSet:
    """ receives a filtered gene set and preforms various enrichment analyses"""
    __slots__ = {'gene_set': 'set of feature names/indices', 'set_name': 'name of the FeatureSet'}
    _go_dicts = {}

    def __init__(self, gene_set: Union[List[str], Set[
        str], filtering.Filter, filtering.CountFilter, filtering.DESeqFilter, filtering.FoldChangeFilter] = None,
                 set_name: str = ''):

        """
        :param gene_set: the set of genomic features to be used in downstream analyses
        :type gene_set: filtering.Filter object, set of strings or list of strings
        :param set_name: name of the FeatureSet
        :type set_name: str


        :Examples:
            >>> from rnalysis import enrichment, filtering
            >>> my_set = enrichment.FeatureSet({'gene1','gene2','gene2'}, 'name of my set')

            >>> filter_obj = filtering.CountFilter('tests/counted.csv')
            >>> my_other_set = enrichment.FeatureSet(filter_obj, 'name of my other set')

        """
        if gene_set is None:
            self.gene_set = general.parse_wbgene_string(input(
                "Please insert genomic features/indices separated by newline "
                "(example: \n'WBGene00000001\nWBGene00000002\nWBGene00000003')"))
        elif isinstance(gene_set, set):
            pass
        elif isinstance(gene_set, list):
            gene_set = set(gene_set)
        elif issubclass(gene_set.__class__, filtering.Filter):
            gene_set = gene_set.index_set
        else:
            raise TypeError(f"Error: 'gene_set' must be a set, list or tuple! Is a {type(gene_set)} instead. ")
        self.gene_set = gene_set
        self.set_name = set_name

    def __repr__(self):
        return f"FeatureSet: {self.set_name}\n" + self.gene_set.__str__()

    @staticmethod
    def _from_string(msg: str = '', del_spaces: bool = False, delimiter: str = '\n'):

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

    def _inplace(self, gene_set: set, inplace: bool):

        """
        Executes the user's choice whether to perform set operations in-place \
        or create a new instance of the FeatureSet object.

        :param gene_set: The set of features resulting from the set operations
        :param inplace: bool. If True, gene_set will be saved to the current FeatureSet object. \
        If False, gene_set will be used to created a new instance of FeatureSet.
        :return: If inplace is True, returns a new instance of FeatureSet.

        """
        if inplace:
            self.gene_set = gene_set
        else:
            return FeatureSet(gene_set)

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
        :type other: FeatureSet, set or str
        :param other: Other object to perform set operation with.
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
            elif isinstance(other, str):
                others[i] = general.parse_wbgene_string(other)
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

    def union(self, *others, inplace: bool = True):

        """
         Calculates the set union of the indices from multipple FeatureSet objects \
        (the indices that exist in at least one of the FeatureSet objects).

        :type others: FeatureSet, set or str
        :param others: The objects against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of FeatureSet. \
        If False, returns a new instance of FeatureSet.
        :return: if inplace is False, returns a new instance of FeatureSet.

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
        return self._inplace(self._set_ops(others, set.union), inplace)

    def intersection(self, *others, inplace: bool = True):

        """
        Calculates the set intersection of the indices from multiple FeatureSet objects \
        (the indices that exist in ALL of the FeatureSet objects).

        :type others: FeatureSet, set or str
        :param others: The objects against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of FeatureSet. \
        If False, returns a new instance of FeatureSet.
        :return: if inplace is False, returns a new instance of FeatureSet.

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
        return self._inplace(self._set_ops(others, set.intersection), inplace)

    def difference(self, *others, inplace: bool = True):

        """
        Calculates the set difference of the indices from multiple FeatureSet objects \
        (the indices that appear in the first FeatureSet object but NOT in the other objects).

        :type others: FeatureSet, set or str
        :param others: The objects against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of FeatureSet. \
        If False, returns a new instance of FeatureSet.
        :return: if inplace is False, returns a new instance of FeatureSet.

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
        return self._inplace(self._set_ops(others, set.difference), inplace)

    def symmetric_difference(self, other, inplace: bool = True):

        """
        Calculates the set symmetric difference of the indices from two FeatureSet objects \
        (the indices that appear in EXACTLY ONE of the FeatureSet objects, and not both/neither). \
        A-symmetric difference-B is equivalent to (A-difference-B)-union-(B-difference-A).

        :type other: FeatureSet, set or str
        :param other: A second  object against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of FeatureSet. \
        If False, returns a new instance of FeatureSet.
        :return: if inplace is False, returns a new instance of FeatureSet.

        :Examples:
            >>> from rnalysis import enrichment
            >>> en = enrichment.FeatureSet({'WBGene00000001','WBGene00000002','WBGene00000006'}, 'set name')
            >>> en2 = enrichment.FeatureSet({'WBGene00000004','WBGene00000001'})
            >>> en.symmetric_difference(en2)
            >>> print(en)
            FeatureSet: set name
            {'WBGene00000002', 'WBGene00000006', 'WBGene00000004'}

        """
        return self._inplace(self._set_ops([other], set.symmetric_difference), inplace)

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
        general.save_to_csv(df, filename=fname + '.csv')

    def go_enrichment(self, mode: str = 'all', alpha: float = 0.05, save_csv: bool = False, fname: str = None):

        """
        Analyzes GO, Tissue and/or Phenotype enrichment for the given group of genomic features. \
        Uses the the Anatomy, Phenotype and Gene Ontology annotations for C. elegans. \
        Corrected p-values are calculated using hypergeometric statistics. \
        For more details see GitHub page of the developers: https://github.com/dangeles/TissueEnrichmentAnalysis

        :type mode: 'go', 'tissue', 'phenotype' or 'all' (default 'all')
        :param mode: the enrichment you wish to perform. 'go' for gene ontology enrichment, \
        'tissue' for tissue enrichment, 'phenotype' for phenotype enrichment, or 'all' for all three.
        :type alpha: float between 0 and 1 (default 0.05)
        :param alpha: Significance threshold.
        :type save_csv: bool (default False)
        :param save_csv: If True, save the result to a csv.
        :type fname: str or pathlib.Path
        :param fname: Name and path in which to save the results. Must be specified if save_csv is True.
        :return: a DataFrame which contains the significant enrichmenet terms

        .. figure::  go_en.png
           :align:   center
           :scale: 40 %

           Example plot of GO enrichment

        .. figure::  tissue_en.png
           :align:   center
           :scale: 40 %

           Example plot of Tissue enrichment

        """
        assert isinstance(alpha, float), "alpha must be a float!"
        assert isinstance(mode, str), "'mode' must be a string!"
        plt.style.use('seaborn-white')
        if mode == 'all':
            d = []
            df_comb = pd.DataFrame()
            for k, arg in enumerate(('go', 'tissue', 'phenotype')):
                print(f'Calculating... {100 * k / 3 :.2f}% done')
                if arg in FeatureSet._go_dicts:
                    d.append(FeatureSet._go_dicts[arg])
                else:
                    d.append(tea.fetch_dictionary(arg))
                    FeatureSet._go_dicts[arg] = d[-1]
                df = tea.enrichment_analysis(self.gene_set, d[-1], alpha=alpha)
                if not df.empty:
                    df_comb = df_comb.append(df)
                    plt.figure()
                    tea.plot_enrichment_results(df, title=f'{arg.capitalize()} Enrichment Analysis', analysis=arg)
                    plt.title(f'{arg.capitalize()} Enrichment Analysis for sample {self.set_name}', fontsize=20)

        else:
            assert (mode == 'go' or mode == 'tissue' or mode == 'phenotype'), "Invalid mode!"
            d = tea.fetch_dictionary(mode)
            df_comb = tea.enrichment_analysis(self.gene_set, d, show=True)
            if not df_comb.empty:
                tea.plot_enrichment_results(df_comb, title=f'{mode.capitalize()} Enrichment Analysis', analysis=mode)
                plt.title(f'{mode.capitalize()} Enrichment Analysis', fontsize=20)

        if save_csv:
            self._enrichment_save_csv(df_comb, fname)
        plt.show()
        return df_comb

    @staticmethod
    def _single_enrichment(gene_set, attributes, attr_ref_df: pd.DataFrame, fraction: Callable, reps: int):
        attributes = [attributes] if not isinstance(attributes, list) else attributes
        for attribute in attributes:
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            df = attr_ref_df[[attribute, 'int_index']]
            srs = df[attribute]
            srs_int = (df.set_index('int_index', inplace=False))[attribute]
            obs_srs = srs.loc[gene_set]
            n = obs_srs.shape[0]
            expected_fraction = fraction(srs)
            observed_fraction = fraction(obs_srs)
            log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
            ind = srs_int.index
            if log2_fold_enrichment >= 0:
                success = sum(
                    (fraction(srs_int.loc[np.random.choice(ind, n, replace=False)]) >= observed_fraction
                     for _ in repeat(None, reps)))
            else:
                success = sum(
                    (fraction(srs_int.loc[np.random.choice(ind, n, replace=False)]) <= observed_fraction
                     for _ in repeat(None, reps)))
            pval = (success + 1) / (reps + 1)

            return [attribute, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval]

    @staticmethod
    def _enrichment_get_attrs(attributes, attr_ref_path):
        if attributes is None:
            attributes = FeatureSet._from_string(
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
        except:
            raise ValueError(f"Invalid or nonexistent Attribute Reference Table path! path:'{attr_ref_path}'")
        if all_attrs[-1].endswith('\n'):
            all_attrs[-1] = all_attrs[-1][:-1]

        if attributes == ['all']:
            attributes = all_attrs
        elif np.all([True if isinstance(i, int) else False for i in attributes]):
            select_attributes = []
            for i in attributes:
                select_attributes.append(all_attrs[i])
            return select_attributes
        return attributes

    def _enrichment_get_reference(self, biotype, background_genes, attr_ref_path, biotype_ref_path):
        gene_set = self.gene_set

        attr_ref_df = general.load_csv(attr_ref_path)
        general._attr_table_assertions(attr_ref_df)
        attr_ref_df.set_index('gene', inplace=True)

        assert (isinstance(biotype, (str, list, set, tuple)))

        if background_genes is None:
            pass
        else:
            assert isinstance(background_genes,
                              (set, FeatureSet)) or issubclass(background_genes.__class__,
                                                               filtering.Filter), f"background_genes must be a set, " \
                                                                                  f"enrichment.FeatureSet or filtering.Filter;" \
                                                                                  f" instead is {type(background_genes)}"
            if isinstance(background_genes, FeatureSet):
                background_genes = background_genes.gene_set
            elif issubclass(background_genes.__class__, filtering.Filter):
                background_genes = background_genes.index_set
            if biotype != 'all':
                warnings.warn(
                    "both 'biotype' and 'background_genes' were specified. Therefore 'biotype' is ignored. ")
                biotype = 'all'

            attr_ref_df = attr_ref_df.loc[background_genes.intersection(set(attr_ref_df.index))]
            if len(attr_ref_df.index) < len(background_genes):
                warnings.warn(
                    f"{len(background_genes) - len(attr_ref_df.index)} indices from the requested "
                    f"background genes do not appear in the Attribute Reference Table, and are therefore ignored. \n"
                    f"This leaves a total of {len(attr_ref_df.index)} background genes. ")
        if biotype == 'all':
            pass
        else:
            biotype_ref_df = general.load_csv(biotype_ref_path)
            general._biotype_table_assertions(biotype_ref_df)
            biotype_ref_df.set_index('gene', inplace=True)
            biotype_ref_df.columns = biotype_ref_df.columns.str.lower()
            if isinstance(biotype, (list, tuple, set)):
                mask = pd.Series(np.zeros_like(biotype_ref_df['biotype'].values, dtype=bool),
                                 biotype_ref_df['biotype'].index,
                                 name='biotype')
                for bio in biotype:
                    mask = mask | (biotype_ref_df['biotype'] == bio)
            else:
                biotype_ref_df = biotype_ref_df.loc[biotype_ref_df.index.intersection(attr_ref_df.index)]
                mask = biotype_ref_df['biotype'] == biotype
            attr_ref_df = attr_ref_df.loc[biotype_ref_df[mask].index]
        attr_ref_df.sort_index(inplace=True)
        attr_ref_df['int_index'] = [i for i in range(len(attr_ref_df.index))]
        print(f"{len(attr_ref_df.index)} background genes are used. ")

        not_in_bg = gene_set.difference(set(attr_ref_df.index))
        if len(not_in_bg) > 0:
            gene_set = gene_set.difference(not_in_bg)
            warnings.warn(f"{len(not_in_bg)} genes in the enrichment set do not appear in the background genes. \n"
                          f"Enrichment will be run on the remaining {len(gene_set)}.")
        return attr_ref_df, gene_set

    def enrich_randomization_parallel(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None,
                                      fdr: float = 0.05, reps: int = 10000, biotype: str = 'protein_coding',
                                      background_genes=None, attr_ref_path: str = 'predefined',
                                      biotype_ref_path: str = 'predefined', save_csv: bool = False, fname=None,
                                      return_fig: bool = False, random_seed: int = None):

        """
        Calculates enrichment scores, p-values and adjusted p-values \
        for enrichment and depletion of selected attributes from an Attribute Reference Table using parallel processing. \
        Background set is determined by either the input variable 'background_genes', \
        or by the input variable 'biotype' and a Biotype Reference Table. \
        Parallel processing makes this function generally faster than FeatureSet.enrich_randomization. \
        To use it you must first start a parallel session, using rnalysis.general.start_parallel_session(). \
        P-values are calculated using a randomization test with the formula p = (successes + 1)/(repeats + 1). \
        P-values are corrected for multiple comparisons using \
        the Benjamini–Hochberg step-up procedure (original FDR method). \
        Enrichment/depletion is determined automatically by the calculated enrichment score: \
        if log2(enrichment score) is positive then enrichment is assumed, \
        and if log2(enrichment score) is negative then depletion is assumed. \
        In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

       :type attributes: str, int, iterable (list, tuple, set, etc) of str/int, or 'all'
       :param attributes: An iterable of attribute names or attribute numbers \
       (according to their order in the Attribute Reference Table). \
       If 'all', all of the attributes in the Attribute Reference Table will be used. \
       If None, a manual input prompt will be raised.
       :type fdr: float between 0 and 1
       :param fdr: Indicates the FDR threshold for significance.
       :type reps: int larger than 0
       :param reps: How many repetitions to run the randomization for. \
       10,000 is the default. Recommended 10,000 or higher.
       :type attr_ref_path: str or pathlib.Path (default 'predefined')
       :param attr_ref_path: the path of the Attribute Reference Table from which user-defined attributes will be drawn.
       :type biotype_ref_path: str or pathlib.Path (default 'predefined')
       :param biotype_ref_path: the path of the Biotype Reference Table. \
       Will be used to generate background set if 'biotype' is specified.
       :type biotype: str specifying a specific biotype, or 'all'. Default 'protein_coding'.
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
       r'C:\dir\file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
       :type return_fig: bool (default False)
       :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
       :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
       :return:
       a pandas DataFrame with the indicated attribute names as rows/index, and the columns 'log2_fold_enrichment'
       and 'pvalue'; and a matplotlib Figure, if 'return_figure' is set to True.

       .. figure::  enrichment_randomization.png
          :align:   center
          :scale: 40 %

          Example plot of enrich_randomization_parallel()
       """
        attr_ref_path = general._get_attr_ref_path(attr_ref_path)
        biotype_ref_path = general._get_biotype_ref_path(biotype_ref_path)
        attr_ref_df, gene_set = self._enrichment_get_reference(biotype=biotype, background_genes=background_genes,
                                                               attr_ref_path=attr_ref_path,
                                                               biotype_ref_path=biotype_ref_path)

        attributes = self._enrichment_get_attrs(attributes=attributes, attr_ref_path=attr_ref_path)
        fraction = lambda mysrs: (mysrs.shape[0] - mysrs.isna().sum()) / mysrs.shape[0]
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
        fraction_rep = list(repeat(fraction, k))
        reps_rep = list(repeat(reps, k))

        res = dview.map(FeatureSet._single_enrichment, gene_set_rep, attributes, attr_ref_df_rep, fraction_rep,
                        reps_rep)
        enriched_list = res.result()
        res_df = pd.DataFrame(enriched_list,
                              columns=['name', 'samples', 'n obs', 'n exp', 'log2_fold_enrichment',
                                       'pval'])
        res_df.replace(-np.inf, -np.max(np.abs(res_df['log2_fold_enrichment'].values)))
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.set_index('name', inplace=True)

        fig = self._plot_enrich_randomization(res_df, title=self.set_name)

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        if return_fig:
            return res_df, fig
        return res_df

    def enrich_randomization(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                             reps: int = 10000, biotype: str = 'protein_coding', background_genes=None,
                             attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                             save_csv: bool = False, fname=None, return_fig: bool = False, random_seed: int = None):

        """
        Calculates enrichment scores, p-values and adjusted p-values \
        for enrichment and depletion of selected attributes from an Attribute Reference Table. \
        Background set is determined by either the input variable 'background_genes', \
        or by the input variable 'biotype' and a Biotype Reference Table. \
        P-values are calculated using a randomization test with the formula p = (successes + 1)/(repeats + 1). \
        P-values are corrected for multiple comparisons using \
        the Benjamini–Hochberg step-up procedure (original FDR method). \
        Enrichment/depletion is determined automatically by the calculated enrichment score: \
        if log2(enrichment score) is positive then enrichment is assumed, \
        and if log2(enrichment score) is negative then depletion is assumed. \
        In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
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
        r'C:\dir\file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
       :type return_fig: bool (default False)
       :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return: a pandas DataFrame with the indicated attribute names as rows/index, and the columns 'log2_fold_enrichment'
        and 'pvalue'; and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  enrichment_randomization.png
           :align:   center
           :scale: 40 %

           Example plot of enrich_randomization()

        """
        attr_ref_path = general._get_attr_ref_path(attr_ref_path)
        biotype_ref_path = general._get_biotype_ref_path(biotype_ref_path)
        attr_ref_df, gene_set = self._enrichment_get_reference(biotype=biotype, background_genes=background_genes,
                                                               attr_ref_path=attr_ref_path,
                                                               biotype_ref_path=biotype_ref_path)
        attributes = self._enrichment_get_attrs(attributes=attributes, attr_ref_path=attr_ref_path)
        fraction = lambda mysrs: (mysrs.shape[0] - mysrs.isna().sum()) / mysrs.shape[0]
        enriched_list = []
        if random_seed is not None:
            assert isinstance(random_seed, int) and random_seed >= 0, f"random_seed must be a non-negative integer. " \
                                                                      f"Value {random_seed} invalid."
            random.seed(random_seed)

        for k, attribute in enumerate(attributes):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            print(f"Finished {k} attributes out of {len(attributes)}")
            df = attr_ref_df[[attribute, 'int_index']]
            srs = df[attribute]
            srs_int = (df.set_index('int_index', inplace=False))[attribute]
            obs_srs = srs.loc[gene_set]
            n = obs_srs.shape[0]
            expected_fraction = fraction(srs)
            observed_fraction = fraction(obs_srs)
            log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
            ind = set(srs_int.index)
            if log2_fold_enrichment >= 0:
                success = sum(
                    (fraction(srs_int.loc[random.sample(ind, n)]) >= observed_fraction
                     for _ in repeat(None, reps)))
            else:
                success = sum(
                    (fraction(srs_int.loc[random.sample(ind, n)]) <= observed_fraction
                     for _ in repeat(None, reps)))
            pval = (success + 1) / (reps + 1)

            enriched_list.append(
                (attribute, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval))

        res_df = pd.DataFrame(enriched_list,
                              columns=['name', 'samples', 'n obs', 'n exp', 'log2_fold_enrichment',
                                       'pval'])
        res_df.replace(-np.inf, -np.max(np.abs(res_df['log2_fold_enrichment'].values)))
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.set_index('name', inplace=True)

        fig = self._plot_enrich_randomization(res_df, title=self.set_name)

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        if return_fig:
            return res_df, fig
        return res_df

    def enrich_hypergeometric(self, attributes: Union[Iterable[str], str, Iterable[int], int] = None, fdr: float = 0.05,
                              biotype: str = 'protein_coding', background_genes=None,
                              attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                              save_csv: bool = False, fname=None, return_fig: bool = False):

        """
        Calculates enrichment scores, p-values and adjusted p-values \
        for enrichment and depletion of selected attributes from an Attribute Reference Table, \
        based on a hypergeometric test. \
        Background set is determined by either the input variable 'background_genes', \
        or by the input variable 'biotype' and a Biotype Reference Table. \
        P-values are calculated using a hypergeometric test: \
        Given M genes in the background set, n genes in the test set, \
        with N genes from the background set belonging to a specific attribute (or 'success') \
        and X genes from the test set belonging to that attribute. \
        If we were to randomly draw n genes from the background set (without replacement), \
        what is the probability of drawing X or more (in case of enrichment)/X or less (in case of depletion) \
        genes belonging to the given attribute? \
        P-values are corrected for multiple comparisons using \
        the Benjamini–Hochberg step-up procedure (original FDR method). \
        Enrichment/depletion is determined automatically by the calculated enrichment score: \
        if log2(enrichment score) is positive then enrichment is assumed, \
        and if log2(enrichment score) is negative then depletion is assumed. \
        In plots, for the clarity of display, complete depletion (linear enrichment = 0) \
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
        r'C:\dir\file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
       :type return_fig: bool (default False)
       :param return_fig: if True, returns a matplotlib Figure object in addition to the results DataFrame.
        :rtype: pd.DataFrame (default) or Tuple[pd.DataFrame, matplotlib.figure.Figure]
        :return:         a pandas DataFrame with the indicated attribute names as rows/index, and the columns 'log2_fold_enrichment'
        and 'pvalue'; and a matplotlib Figure, if 'return_figure' is set to True.

        .. figure::  enrichment_randomization.png
           :align:   center
           :scale: 40 %

           Example plot of enrich_hypergeometric()

        """
        attr_ref_path = general._get_attr_ref_path(attr_ref_path)
        biotype_ref_path = general._get_biotype_ref_path(biotype_ref_path)
        attr_ref_df, gene_set = self._enrichment_get_reference(biotype=biotype, background_genes=background_genes,
                                                               attr_ref_path=attr_ref_path,
                                                               biotype_ref_path=biotype_ref_path)
        attributes = self._enrichment_get_attrs(attributes=attributes, attr_ref_path=attr_ref_path)
        fraction = lambda mysrs: (mysrs.shape[0] - mysrs.isna().sum()) / mysrs.shape[0]
        enriched_list = []
        for k, attribute in enumerate(attributes):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            print(f"Finished {k} attributes out of {len(attributes)}")
            df = attr_ref_df[[attribute, 'int_index']]
            srs = df[attribute]
            obs_srs = srs.loc[gene_set]
            n = obs_srs.shape[0]
            expected_fraction = fraction(srs)
            observed_fraction = fraction(obs_srs)
            log2_fold_enrichment = np.log2(observed_fraction / expected_fraction) if observed_fraction > 0 else -np.inf
            pval = self._calc_hypergeometric_pval(bg_size=srs.shape[0], go_size=srs.notna().sum(),
                                                  de_size=obs_srs.shape[0], go_de_size=obs_srs.notna().sum())

            enriched_list.append(
                (attribute, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval))

        res_df = pd.DataFrame(enriched_list,
                              columns=['name', 'samples', 'n obs', 'n exp', 'log2_fold_enrichment',
                                       'pval'])
        res_df.replace(-np.inf, -np.max(np.abs(res_df['log2_fold_enrichment'].values)))
        significant, padj = multitest.fdrcorrection(res_df['pval'].values, alpha=fdr)
        res_df['padj'] = padj
        res_df['significant'] = significant
        res_df.set_index('name', inplace=True)

        fig = self._plot_enrich_randomization(res_df, title=self.set_name)

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        if return_fig:
            return res_df, fig
        return res_df

    @staticmethod
    def _plot_enrich_randomization(df: pd.DataFrame, title: str = ''):

        """
        Receives a DataFrame output from FeatureSet.enrich_randomization, and plots it in a bar plort \
        Static class method. \
        For the clarity of display, complete depletion (linear enrichment = 0) \
        appears with the smallest value in the scale.

        :param df: a pandas DataFrame created by FeatureSet.enrich_randomization.
        :param title: plot title.
        :return: a matplotlib.pyplot.bar instance

        """
        plt.style.use('seaborn-white')

        enrichment_names = df.index.values.tolist()
        enrichment_pvalue = df['padj']
        # set enrichment scores which are 'inf' or '-inf' to be the second highest/lowest enrichment score in the list
        enrichment_scores = df['log2_fold_enrichment'].values.copy()
        scores_no_inf = [i for i in enrichment_scores if i != np.inf and i != -np.inf and i < 0]
        if len(scores_no_inf) == 0:
            scores_no_inf.append(-1)
        for i in range(len(enrichment_scores)):
            if enrichment_scores[i] == -np.inf:
                enrichment_scores[i] = min(scores_no_inf)

        # get color values for bars
        data_color = [(i / 3) * 127.5 for i in enrichment_scores]
        data_color_norm = [i + 127.5 for i in data_color]
        data_color_norm_256 = [int(i) if i != np.inf and i != -np.inf else np.sign(i) * max(np.abs(scores_no_inf)) for i
                               in data_color_norm]
        my_cmap = plt.cm.get_cmap('coolwarm')
        colors = my_cmap(data_color_norm_256)

        # generate bar plot
        fig, ax = plt.subplots(constrained_layout=True, figsize=[6.4 * 0.5 + 0.5 * df.shape[0], 5.6])
        bar = ax.bar(x=range(len(enrichment_names)), height=enrichment_scores, color=colors, edgecolor='black',
                     linewidth=1)
        bar.tick_labels = enrichment_names
        # add horizontal line
        ax.axhline(color='black', linewidth=1)
        # add colorbar
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(3, -3))
        sm.set_array([])
        cbar = fig.colorbar(sm)
        cbar.set_label('Colorbar', fontsize=12)
        # apply xticks
        ax.set_xticks(range(len(enrichment_names)))
        ax.set_xticklabels(enrichment_names, fontsize=13, rotation=45)
        # ylabel and title
        ax.set_ylabel(r"$\log_2$(Fold Enrichment)", fontsize=14)
        ax.set_title(title, fontsize=16)
        # add significance asterisks
        for col, sig in zip(bar, enrichment_pvalue):
            fontweight = 'bold'
            if sig < 0.0001:
                asterisks = u'\u2217' * 4
            elif sig < 0.001:
                asterisks = u'\u2217' * 3
            elif sig < 0.01:
                asterisks = u'\u2217' * 2
            elif sig < 0.05:
                asterisks = u'\u2217'
            else:
                asterisks = 'ns'
                fontweight = 'normal'
            valign = 'bottom' if np.sign(col._height) == 1 else 'top'
            ax.text(x=col.xy[0] + 0.5 * col._width,
                    y=col._height, s=asterisks, fontname='DejaVu Sans', fontweight=fontweight,
                    fontsize=12, horizontalalignment='center', verticalalignment=valign)

        sns.despine()
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
            >>> d = filtering.Filter("tests/test_deseq.csv")
            >>> en = enrichment.FeatureSet(d)
            >>> en.biotypes(ref='tests/biotype_ref_table_for_tests.csv')
                            gene
            biotype
            protein_coding    26
            pseudogene         1
            unknown            1

        """

        ref = general._get_biotype_ref_path(ref)
        ref_df = general.load_csv(ref)
        general._biotype_table_assertions(ref_df)
        ref_df.columns = ref_df.columns.str.lower()
        not_in_ref = pd.Index(self.gene_set).difference(set(ref_df['gene']))
        if len(not_in_ref) > 0:
            warnings.warn(
                f'{len(not_in_ref)} of the features in the Filter object do not appear in the Biotype Reference Table. ')
            ref_df = ref_df.append(pd.DataFrame({'gene': not_in_ref, 'biotype': 'not_in_biotype_reference'}))
        return ref_df.set_index('gene', drop=False).loc[self.gene_set].groupby('biotype').count()


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
        elif issubclass(objs[obj].__class__, filtering.Filter):
            objs[obj] = objs[obj].index_set
        elif isinstance(objs[obj], FeatureSet):
            objs[obj] = objs[obj].gene_set
        elif isinstance(objs[obj], str):
            try:
                attr_table
            except NameError:
                pth = general._get_attr_ref_path(ref)
                attr_table = general.load_csv(pth, 0)
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
    fig = plt.figure()
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
    multi_ind = pd.MultiIndex.from_product([[True, False] for i in range(len(names))], names=names)[:-1]
    srs = pd.Series(index=multi_ind)
    for ind in multi_ind:
        sets = list(compress(names, ind))
        group_size = len(set.intersection(*[objs[s] for s in sets]))
        srs.loc[ind] = group_size
    return srs

# TODO: other types of plots
# TODO: heat map plot of multiple DESEQ files
