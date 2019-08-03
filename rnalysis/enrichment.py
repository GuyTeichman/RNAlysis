"""
This module receives a set of feature names, and can perform various enrichment analyses on them, such as \
GO enrichment, tissue/phenotype enrichment, enrichment and depletion of big table columns, etc. \
Furthermore, set operations (union, intersect, difference, symmetric difference) can be performed on two \
EnrichmentProcessing objects. \
Results of enrichment analyses can be saved to .csv files.
"""

import os
import numpy as np
import pandas as pd
from rnalysis import general, __gene_names_and_biotype__
import tissue_enrichment_analysis as tea
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from pathlib import Path
import statsmodels.stats.multitest as multitest


class EnrichmentProcessing:
    """ receives a filtered gene set and preforms various enrichment analyses"""
    _go_dicts = {}

    def __init__(self, gene_set: set = None, set_name: str = ''):
        if gene_set is None:
            self.gene_set = general.parse_wbgene_string(input(
                "Please insert WBGenes separated by newline "
                "(example: \n'WBGene00000001\nWBGene00000002\nWBGene00000003')"))
        elif isinstance(gene_set, set):
            pass
        elif isinstance(gene_set, (list, tuple)):
            gene_set = set(gene_set)
        else:
            raise TypeError(f"Error: 'gene_set' must be a set, list or tuple! Is a {type(gene_set)} instead. ")
        self.gene_set = gene_set
        self.set_name = set_name

    @staticmethod
    def _from_string(msg: str = '', del_spaces: bool = False, delimiter: str = '\n'):
        """
        Takes a manual string input from the user, and then splits it using a comma delimiter into a list of values. \
        Called when an EnrichmentProcessing instance is created without input, \
        or when EnrichmentProcessing.enrich_big_table is called without input.

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
        or create a new instance of the EnrichmentProcessing object.

        :param gene_set: The set of features resulting from the set operations
        :param inplace: bool. If True, gene_set will be saved to the current EnrichmentProcessing object. \
        If False, gene_set will be used to created a new instance of EnrichmentProcessing.
        :return:
        If inplace is True, returns a new instance of EnrichmentProcessing.
        """
        if inplace:
            self.gene_set = gene_set
        else:
            return EnrichmentProcessing(gene_set)

    @staticmethod
    def _get_ref_path(ref):
        if ref == 'predefined':
            return general.read_reference_table_path()
        else:
            return ref

    def save_txt(self, fname):
        """
        Save the list of features in the EnrichmentProcessing object under the specified filename and path.

        :type fname: str, pathlib.Path
        :param fname: full filename/path for the output file.
        """
        assert isinstance(fname, (str, Path)), "fname must be str or pathlib.Path!"
        with open(fname, 'w') as f:
            for gene in self.gene_set:
                f.write(gene + '\n')

    def _set_ops(self, others, op):
        """
        Performs a given set operation on self and on another object (EnrichmentProcessing or set).
        :type other: EnrichmentProcessing, set or str
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
            elif isinstance(other, EnrichmentProcessing):
                others[i] = other.gene_set
            elif isinstance(other, str):
                others[i] = general.parse_wbgene_string(other)
            else:
                raise TypeError("'other' must be an EnrichmentProcessing object or a set!")
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
         Calculates the set union of the WBGene indices from multipple EnrichmentProcessing objects \
        (the indices that exist in at least one of the EnrichmentProcessing objects).

        :type others: EnrichmentProcessing, set or str
        :param others: The objects against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of EnrichmentProcessing. \
        If False, returns a new instance of EnrichmentProcessing.
        :return:
        if inplace is False, returns a new instance of EnrichmentProcessing.
        """
        return self._inplace(self._set_ops(others, set.union), inplace)

    def intersection(self, *others, inplace: bool = True):
        """
        Calculates the set intersection of the WBGene indices from multiple EnrichmentProcessing objects \
        (the indices that exist in ALL of the EnrichmentProcessing objects).

        :type others: EnrichmentProcessing, set or str
        :param others: The objects against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of EnrichmentProcessing. \
        If False, returns a new instance of EnrichmentProcessing.
        :return:
        if inplace is False, returns a new instance of EnrichmentProcessing.
                """
        return self._inplace(self._set_ops(others, set.intersection), inplace)

    def difference(self, *others, inplace: bool = True):
        """
        Calculates the set difference of the WBGene indices from multiple EnrichmentProcessing objects \
        (the indices that appear in the first EnrichmentProcessing object but NOT in the other objects).

        :type others: EnrichmentProcessing, set or str
        :param others: The objects against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of EnrichmentProcessing. \
        If False, returns a new instance of EnrichmentProcessing.
        :return:
        if inplace is False, returns a new instance of EnrichmentProcessing.
        """
        return self._inplace(self._set_ops(others, set.difference), inplace)

    def symmetric_difference(self, other, inplace: bool = True):
        """
        Calculates the set symmetric difference of the WBGene indices from two EnrichmentProcessing objects \
        (the indices that appear in EXACTLY ONE of the EnrichmentProcessing objects, and not both/neither). \
        A-symmetric difference-B is equivalent to (A-difference-B)-union-(B-difference-A).

        :type other: EnrichmentProcessing, set or str
        :param other: A second  object against which the current object will be compared.
        :type inplace: bool
        :param inplace: If True (default), modifies the current instance of EnrichmentProcessing. \
        If False, returns a new instance of EnrichmentProcessing.
        :return:
        if inplace is False, returns a new instance of EnrichmentProcessing.
        """
        return self._inplace(self._set_ops([other], set.symmetric_difference), inplace)

    @staticmethod
    def _enrichment_save_csv(df: pd.DataFrame, fname: str, suffix: str = ''):
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
        general.save_to_csv(df, filename=fname + '.csv', suffix=suffix)

    def go_enrichment(self, mode: str = 'go', alpha: float = 0.05, save_csv: bool = False, fname: str = None):
        """
        Analyzes GO, Tissue and/or Phenotype enrichment of the given group of features. \
        Uses the the Anatomy, Phenotype and Gene Ontology annotations for C. elegans. \
        Corrected p-values are calculated using hypergeometric statistics. \
        For more details see GitHub page of the developers: https://github.com/dangeles/TissueEnrichmentAnalysis

        :param mode: the enrichment you wish to perform. 'go' for gene ontology enrichment, \
        'tissue' for tissue enrichment, 'phenotype' for phenotype enrichment.
        :param alpha: float. Significance threshold. Default is 0.05
        :param save_csv: bool. False by default. If True, save the result to a csv.
        :param fname: Name and path in which to save the results. Must be filled if save_csv is True.
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
        if mode == 'all':
            d = []
            df_comb = pd.DataFrame()
            for k, arg in enumerate(('go', 'tissue', 'phenotype')):
                print(f'Calculating... {100 * k / 3 :.2f}% done')
                if arg in EnrichmentProcessing._go_dicts:
                    d.append(EnrichmentProcessing._go_dicts[arg])
                else:
                    d.append(tea.fetch_dictionary(arg))
                    EnrichmentProcessing._go_dicts[arg] = d[-1]
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

    def enrich_big_table(self, attributes: list = None, fdr: float = 0.05, reps=10000, biotype: str = 'protein_coding',
                         big_table_pth: str = 'predefined', save_csv: bool = False, fname=None):
        """
        Calculates enrichment scores, p-values and q-values \
        for enrichment and depletion of selected attributes from the Big Table. \
        P-values are calculated using a randomization test, and corrected for multiple comparisons using \
        the Benjamini–Hochberg step-up procedure (original FDR method). \
        Enrichment/depletion is determined automatically by the calculated enrichment score: \
        if log2(enrichment score) is positive then enrichment is assumed, \
        and if log2(enrichment score) is negative then depletion is assumed.

        :param attributes: An iterable of attribute names (strings). If None, a manual input prompt will be raised.
        :param fdr: float. Indicates the FDR threshold for significance.
        :param reps: How many repetitions to run the randomization for. \
        10,000 is the default. Recommended 10,000 or higher.
        :param big_table_pth: the path of the Big Table file to be used as reference.
        :param biotype: the biotype you want your reference to have. 'all' will include all biotypes, \
        'protein_coding' will include only protein-coding genes in the reference, etc.
        :param save_csv: bool. If True, will save the results to a .csv file, under the name specified in 'fname'.
        :param fname: str/Path. The full path and name of the file to which to save the results. For example: \
        r'C:\dir\file'. No '.csv' suffix is required. If None (default), fname will be requested in a manual prompt.
        :return:
        a pandas DataFrame with the indicated attribute names as rows/index, and the columns 'log2_fold_enrichment'
        and 'pvalue'.

        .. figure::  bigtable_en.png
           :align:   center
           :scale: 40 %

           Example plot of big table enrichment
        """
        big_table_pth = self._get_ref_path(big_table_pth)
        if attributes is None:
            attributes = self._from_string(
                "Please insert attributes separated by newline "
                "(for example: \n'epigenetic_related_genes\nnrde-3 targets\nALG-3/4 class small RNAs')")
        elif isinstance(attributes, str):
            attributes = [attributes]
        else:
            assert isinstance(attributes, (list, tuple, set)), "'attributes' must be a list, tuple or set!"

        try:
            big_table = general.load_csv(big_table_pth, 0, drop_gene_names=False)
        except:
            raise ValueError("Invalid or nonexistent big table path!")

        assert (isinstance(biotype, (str, list, set, tuple)))
        if biotype == 'all':
            pass
        else:
            biotype_ref = general.load_csv(__gene_names_and_biotype__, 0, drop_gene_names=False)
            if isinstance(biotype, (list, tuple, set)):
                mask = pd.Series(np.zeros_like(biotype_ref['bioType'].values, dtype=bool), biotype_ref['bioType'].index,
                                 name='bioType')
                for bio in biotype:
                    mask = mask | (biotype_ref['bioType'] == bio)
            else:
                mask = biotype_ref['bioType'] == biotype
            big_table = big_table.loc[biotype_ref[mask].index]

            big_table['int_index'] = [int(i[6:14]) for i in big_table.index]
            fraction = lambda mysrs: (mysrs.shape[0] - mysrs.isna().sum()) / mysrs.shape[0]
            enriched_list = []
            for k, attribute in enumerate(attributes):
                assert isinstance(attribute, str), f"Error in attribute {attribute}: attributes must be strings!"
                print(f"Finished {k} attributes out of {len(attributes)}")
                df = big_table[[attribute, 'int_index']]
                srs = df[attribute]
                srs_int = (df.set_index('int_index', inplace=False))[attribute]
                obs_srs = srs.loc[self.gene_set]
                expected_fraction = fraction(srs)
                observed_fraction = fraction(obs_srs)
                log2_fold_enrichment = np.log2((observed_fraction + 0.0001) / (expected_fraction + 0.0001))
                success = sum((fraction(srs_int.loc[np.random.choice(srs_int.index, obs_srs.shape[0],
                                                                     replace=False)]) >= observed_fraction
                               if log2_fold_enrichment >= 0 else fraction(
                    srs_int.loc[
                        np.random.choice(srs_int.index, obs_srs.shape[0], replace=False)]) <= observed_fraction
                               for _ in range(reps)))
                pval = (success + 1) / (reps + 1)
                n = obs_srs.shape[0]
                enriched_list.append(
                    (attribute, n, int(n * observed_fraction), n * expected_fraction, log2_fold_enrichment, pval))

            enriched_df = pd.DataFrame(enriched_list,
                                       columns=['name', 'samples', 'n obs', 'n exp', 'log2_fold_enrichment',
                                                'pvals'])
            significant, padj = multitest.fdrcorrection(enriched_df['pvals'].values, alpha=fdr)
            enriched_df['padj'] = padj
            enriched_df['significant'] = significant
            enriched_df.set_index('name', inplace=True)

            self._plot_enrich_big_table(enriched_df, title=self.set_name)

            if save_csv:
                self._enrichment_save_csv(enriched_df, fname)
            print(enriched_df)

            return enriched_df

    @staticmethod
    def _plot_enrich_big_table(df: pd.DataFrame, title: str = ''):
        """
        Receives a DataFrame output from EnrichmentProcessing.enrich_big_table, and plots it in a bar plort \
        Static class method.

        :param df: a pandas DataFrame created by EnrichmentProcessing.enrich_big_table.
        :param title: plot title.
        :return:
        a matplotlib.pyplot.bar instance
        """
        enrichment_names = df.index.values.tolist()
        enrichment_scores = df['log2_fold_enrichment']
        enrichment_pvalue = df['padj']
        abs_enrichment_scores = [abs(i) for i in enrichment_scores]
        data_color = [float(i / max(abs_enrichment_scores)) for i in enrichment_scores]
        my_cmap = plt.cm.get_cmap('coolwarm')
        colors = my_cmap(data_color)
        fig, ax = plt.subplots()
        # ax.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        bar = ax.bar(x=range(len(enrichment_names)), height=enrichment_scores, color=colors)
        bar.tick_labels = enrichment_names
        # sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(*ax.get_ylim()))
        absmax = max([abs(i) for i in ax.get_ylim()])
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(absmax, -absmax))
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

# TODO: other types of plots
# TODO: heat map plot of multiple DESEQ files
# TODO: function that prints all biotypes in the sample
# TODO: option to give a list/set of background genes (enrichment)
# TODO: option to remove specific genes from the background genes (enrichment)
