"""
This module receives a set of feature names, and can perform various enrichment analyses on them, such as \
GO enrichment, tissue/phenotype enrichment, enrichment and depletion of Reference Table attributes, etc. \
Furthermore, set operations (union, intersect, difference, symmetric difference) can be performed on two \
FeatureSet objects. \
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


class FeatureSet:
    """ receives a filtered gene set and preforms various enrichment analyses"""
    _go_dicts = {}

    def __init__(self, gene_set: set = None, set_name: str = ''):
        if gene_set is None:
            self.gene_set = general.parse_wbgene_string(input(
                "Please insert genomic features/indices separated by newline "
                "(example: \n'WBGene00000001\nWBGene00000002\nWBGene00000003')"))
        elif isinstance(gene_set, set):
            pass
        elif isinstance(gene_set, (list, tuple)):
            gene_set = set(gene_set)
        elif issubclass(gene_set.__class__, filtering.Filter):
            gene_set = gene_set.index_set
        else:
            raise TypeError(f"Error: 'gene_set' must be a set, list or tuple! Is a {type(gene_set)} instead. ")
        self.gene_set = gene_set
        self.set_name = set_name

    @staticmethod
    def _from_string(msg: str = '', del_spaces: bool = False, delimiter: str = '\n'):
        """
        Takes a manual string input from the user, and then splits it using a comma delimiter into a list of values. \
        Called when an FeatureSet instance is created without input, \
        or when FeatureSet.enrich_randomization is called without input.

        :param msg: a promprt to be printed to the user
        :param del_spaces: if True, will delete all spaces in each delimited value.
        :param delimiter: the delimiter used to separate the values. Default is '\n'
        :return:
        A list of the comma-seperated values the user inserted.
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
        :return:
        If inplace is True, returns a new instance of FeatureSet.
        """
        if inplace:
            self.gene_set = gene_set
        else:
            return FeatureSet(gene_set)

    def save_txt(self, fname):
        """
        Save the list of features in the FeatureSet object under the specified filename and path.

        :type fname: str, pathlib.Path
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

    def _set_ops(self, others, op):
        """
        Performs a given set operation on self and on another object (FeatureSet or set).
        :type other: FeatureSet, set or str
        :param other: Other object to perform set operation with.
        :param op: The set operation to be performed. \
        set.union, set.intersection, set.difference or set.symmetric_difference.
        :return:
        A set resulting from the set operation.
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
        :return:
        if inplace is False, returns a new instance of FeatureSet.
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
        :return:
        if inplace is False, returns a new instance of FeatureSet.
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
        :return:
        if inplace is False, returns a new instance of FeatureSet.
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
        :return:
        if inplace is False, returns a new instance of FeatureSet.
        """
        return self._inplace(self._set_ops([other], set.symmetric_difference), inplace)

    @staticmethod
    def _enrichment_save_csv(df: pd.DataFrame, fname: str):
        """
        Internal method, used to save enrichment results to .csv files. Static class method.

        :param df: pandas DataFrame to be saved.
        :param fname: Name and full path under which the DataFrame will be saved
        :param suffix: Suffix to add to the file name before the .csv.
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
        :return:
        a DataFrame which contains the significant enrichmenet terms

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
    def _single_enrichment(gene_set, attributes, big_table: pd.DataFrame, fraction, reps):
        attributes = [attributes] if not isinstance(attributes, list) else attributes
        for attribute in attributes:
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            df = big_table[[attribute, 'int_index']]
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
        try:
            big_table = general.load_csv(attr_ref_path, 0)
        except:
            raise ValueError(f"Invalid or nonexistent Attribute Reference Table path! path:'{attr_ref_path}'")

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

            big_table = big_table.loc[background_genes.intersection(set(big_table.index))]
            if len(big_table.index) < len(background_genes):
                warnings.warn(
                    f"{len(background_genes) - len(big_table.index)} indices from the requested "
                    f"background genes do not appear in the Attribute Reference Table, and are therefore ignored. \n"
                    f"This leaves a total of {len(big_table.index)} background genes. ")
        if biotype == 'all':
            pass
        else:
            biotype_ref = general.load_csv(biotype_ref_path, 0)
            if isinstance(biotype, (list, tuple, set)):
                mask = pd.Series(np.zeros_like(biotype_ref['bioType'].values, dtype=bool), biotype_ref['bioType'].index,
                                 name='bioType')
                for bio in biotype:
                    mask = mask | (biotype_ref['bioType'] == bio)
            else:
                biotype_ref = biotype_ref.loc[biotype_ref.index.intersection(big_table.index)]
                mask = biotype_ref['bioType'] == biotype
            big_table = big_table.loc[biotype_ref[mask].index]
        big_table.sort_index(inplace=True)
        big_table['int_index'] = [i for i in range(len(big_table.index))]
        print(f"{len(big_table.index)} background genes are used. ")

        not_in_bg = gene_set.difference(set(big_table.index))
        if len(not_in_bg) > 0:
            gene_set = gene_set.difference(not_in_bg)
            warnings.warn(f"{len(not_in_bg)} genes in the enrichment set do not appear in the background genes. \n"
                          f"Enrichment will be run on the remaining {len(gene_set)}.")
        return big_table, gene_set

    def enrich_randomization_parallel(self, attributes: list = None, fdr: float = 0.05, reps=10000,
                                      biotype: str = 'protein_coding', background_genes: set = None,
                                      attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                                      save_csv: bool = False, fname=None):
        """
        Calculates enrichment scores, p-values and q-values \
        for enrichment and depletion of selected attributes from an Attribute Reference Table using parallel processing. \
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
       :param attr_ref_path: the path of the Big Table file to be used as reference.
       :type biotype: str specifying a specific biotype, or 'all'. Default 'protein_coding'.
       :param biotype: the biotype you want your background genes to have. 'all' will include all biotypes, \
       'protein_coding' will include only protein-coding genes in the reference, etc. \
       Cannot be specified together with 'background_genes'.
       :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
       :param background_genes: a set of specific feature indices to be used as background genes. \
       Cannot be specified together with 'biotype'.
       :type save_csv: bool, default False
       :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
       :type fname: str or pathlib.Path
       :param fname: The full path and name of the file to which to save the results. For example: \
       r'C:\dir\file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
       :return:
       a pandas DataFrame with the indicated attribute names as rows/index, and the columns 'log2_fold_enrichment'
       and 'pvalue'.

       .. figure::  bigtable_en.png
          :align:   center
          :scale: 40 %

          Example plot of big table enrichment
       """
        attr_ref_path = general._get_attr_ref_path(attr_ref_path)
        biotype_ref_path = general._get_biotype_ref_path(biotype_ref_path)
        attributes = self._enrichment_get_attrs(attributes=attributes, attr_ref_path=attr_ref_path)
        big_table, gene_set = self._enrichment_get_reference(biotype=biotype, background_genes=background_genes,
                                                             attr_ref_path=attr_ref_path,
                                                             biotype_ref_path=biotype_ref_path)
        if attributes == ['all']:
            attributes = big_table.columns[:-1]
        fraction = lambda mysrs: (mysrs.shape[0] - mysrs.isna().sum()) / mysrs.shape[0]
        client = Client()
        dview = client[:]
        dview.execute("""import numpy as np
              import pandas as pd""")
        k = len(attributes)
        gene_set_rep = list(repeat(gene_set, k))
        big_table_rep = list(repeat(big_table, k))
        fraction_rep = list(repeat(fraction, k))
        reps_rep = list(repeat(reps, k))

        res = dview.map(FeatureSet._single_enrichment, gene_set_rep, attributes, big_table_rep, fraction_rep,
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

        self._plot_enrich_randomization(res_df, title=self.set_name)

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

        return res_df

    def enrich_randomization(self, attributes: list = None, fdr: float = 0.05, reps=10000,
                             biotype: str = 'protein_coding', background_genes: set = None,
                             attr_ref_path: str = 'predefined', biotype_ref_path: str = 'predefined',
                             save_csv: bool = False, fname=None):
        """
        Calculates enrichment scores, p-values and q-values \
        for enrichment and depletion of selected attributes from an Attribute Reference Table. \
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
        :type attr_ref_path: 'predefined' (default), str or pathlib.Path
        :param attr_ref_path: the path of the Attribute Reference Table file to be used as reference.
        :type biotype: str specifying a specific biotype, list/set of strings each specifying a biotype, or 'all'. \
        Default 'protein_coding'.
        :param biotype: the biotype you want your background genes to have. 'all' will include all biotypes, \
        'protein_coding' will include only protein-coding genes in the reference, etc. \
        Cannot be specified together with 'background_genes'.
        :type background_genes: set of feature indices, filtering.Filter object, or enrichment.FeatureSet object
        :param background_genes: a set of specific feature indices to be used as background genes. \
        Cannot be specified together with 'biotype'.
        :type save_csv: bool, default False
        :param save_csv: If True, will save the results to a .csv file, under the name specified in 'fname'.
        :type fname: str or pathlib.Path
        :param fname: The full path and name of the file to which to save the results. For example: \
        r'C:\dir\file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :return:
        a pandas DataFrame with the indicated attribute names as rows/index, and the columns 'log2_fold_enrichment'
        and 'pvalue'.

        .. figure::  bigtable_en.png
           :align:   center
           :scale: 40 %

           Example plot of enrich_randomization
        """
        attr_ref_path = general._get_attr_ref_path(attr_ref_path)
        biotype_ref_path = general._get_biotype_ref_path(biotype_ref_path)
        attributes = self._enrichment_get_attrs(attributes=attributes, attr_ref_path=attr_ref_path)
        big_table, gene_set = self._enrichment_get_reference(biotype=biotype, background_genes=background_genes,
                                                             attr_ref_path=attr_ref_path,
                                                             biotype_ref_path=biotype_ref_path)
        fraction = lambda mysrs: (mysrs.shape[0] - mysrs.isna().sum()) / mysrs.shape[0]
        enriched_list = []
        for k, attribute in enumerate(attributes):
            assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
            print(f"Finished {k} attributes out of {len(attributes)}")
            df = big_table[[attribute, 'int_index']]
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

        self._plot_enrich_randomization(res_df, title=self.set_name)

        if save_csv:
            self._enrichment_save_csv(res_df, fname)

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
        :return:
        a matplotlib.pyplot.bar instance
        """
        plt.style.use('seaborn-white')
        enrichment_names = df.index.values.tolist()
        enrichment_scores = df['log2_fold_enrichment'].values.copy()
        scores_no_inf = [i for i in enrichment_scores if i != np.inf and i != -np.inf and i < 0]
        if len(scores_no_inf) == 0:
            scores_no_inf.append(-1)
        for i in range(len(enrichment_scores)):
            if enrichment_scores[i] == -np.inf:
                enrichment_scores[i] = min(scores_no_inf)
        enrichment_pvalue = df['padj']
        abs_enrichment_scores = [abs(i) for i in enrichment_scores]
        data_color = [(i / 3) * 127.5 for i in enrichment_scores]
        data_color_norm = [i + 127.5 for i in data_color]
        data_color_norm_256 = [int(i) if i != np.inf and i != -np.inf else np.sign(i) * max(np.abs(scores_no_inf)) for i
                               in data_color_norm]
        my_cmap = plt.cm.get_cmap('coolwarm')
        colors = my_cmap(data_color_norm_256)
        fig, ax = plt.subplots()
        # ax.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        bar = ax.bar(x=range(len(enrichment_names)), height=enrichment_scores, color=colors)
        bar.tick_labels = enrichment_names
        # sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(*ax.get_ylim()))
        # absmax = max([abs(i) for i in ax.get_ylim()])
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(3, -3))
        sm.set_array([])
        cbar = fig.colorbar(sm)
        cbar.set_label('', rotation=90, labelpad=25, fontsize=26)
        plt.xticks(range(len(enrichment_names)), enrichment_names, fontsize=13, rotation='vertical')
        plt.ylabel("Log2 Fold Enrichment", fontsize=26)
        for col, sig in zip(bar, enrichment_pvalue):
            fontsize = 21
            if sig < 0.0001:
                asterisks = '****'
            elif sig < 0.001:
                asterisks = '***'
            elif sig < 0.01:
                asterisks = '**'
            elif sig < 0.05:
                asterisks = '*'
            else:
                asterisks = 'ns'
                fontsize = 16
            plt.text(x=col.xy[0] + 0.5 * col._width,
                     y=col._height + (max(enrichment_scores) - min(enrichment_scores)) / 50 * np.sign(col._height),
                     s=asterisks,
                     fontsize=fontsize, horizontalalignment='center', verticalalignment='center')

        sns.despine()
        plt.title(title)
        plt.tight_layout()
        plt.show()
        return bar

    def biotypes(self, ref: str = 'predefined'):
        """
        Returns a DataFrame of the biotypes in the gene set and their count.

        :type ref: str or pathlib.Path
        :param ref: Path of the reference file used to determine biotype. \
        Default is the path predefined in the settings file.
        """

        ref = general._get_biotype_ref_path(ref)
        ref_df = general.load_csv(ref)
        not_in_ref = pd.Index(self.gene_set).difference(set(ref_df['gene']))
        if len(not_in_ref) > 0:
            warnings.warn(
                f'{len(not_in_ref)} of the features in the Filter object do not appear in the biotype reference file. ')
            ref_df = ref_df.append(pd.DataFrame({'gene': not_in_ref, 'bioType': 'not_in_biotype_reference'}))
        return ref_df.set_index('gene', drop=False).loc[self.gene_set].groupby('bioType').count()


def _fetch_sets(objs: dict, ref: str = 'predefined'):
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


def upset_plot(objs: dict, ref: str = 'predefined'):
    """

    :param objs:
    :type objs:
    :param ref:
    :type ref:
    :return:
    :rtype:
    """

    upset_df = _generate_upset_srs(_fetch_sets(objs=objs, ref=ref))
    return upset.plot(upset_df)


def venn_diagram(objs: dict, ref: str = 'predefined', title: str = 'default', set_colors: tuple = ('r', 'g', 'b'),
                 alpha: float = 0.4, weighted: bool = True, lines: bool = True, linecolor: str = 'black',
                 linestyle='solid', linewidth=2.0,
                 normalize_to: float = 1.0):
    """

    :param objs:
    :type objs:
    :param ref:
    :type ref:
    :param title:
    :type title:
    :param set_colors:
    :type set_colors:
    :param alpha:
    :type alpha:
    :param weighted:
    :type weighted:
    :param lines:
    :type lines:
    :param linecolor:
    :type linecolor:
    :param linestyle:
    :type linestyle:
    :param linewidth:
    :type linewidth:
    :param normalize_to:
    :type normalize_to:
    :return:
    :rtype:
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

    :param objs:
    :type objs: dict of sets
    :return:
    :rtype:
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
