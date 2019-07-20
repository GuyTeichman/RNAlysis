"""
This module can append WBGene indices to feature names in a DESeq/HTCount .csv file, \
according to a reference table (such as the Big Table).
"""

import os
import pandas as pd
from pathlib import Path
from rnalysis import general, __gene_names_and_biotype__


class GeneNameTranslator:
    """
    A class that takes a file with gene names that does not have WBGene indices, and a reference file that includes
    gene names and their corresponding WBGene indices, and adds WBGene indices to the original file.
    """

    def __init__(self, target_filename,
                 reference_filename: str = __gene_names_and_biotype__,
                 ref_wbgene_col: str = 'gene', ref_gene_col: str = 'gene ID', ref_seq_col: str = 'Sequence ID',
                 ref_other_col: str = 'Other IDs', target_gene_col: str = 'genes'):
        """

        :param target_filename: name of the file for which WBGene indices need to be added.
        :param reference_filename: name of the reference file to take WBGene indices from.
        :param ref_wbgene_col: The name of the column in the reference file in which WBGene indices appear
        :param ref_gene_col: the name of the column in the reference file in which gene symbols appear (like 'rde-4')
        :param ref_seq_col: the name of the column in the reference file in which sequence names appear
        :param ref_other_col: the name of the column in the reference file in which other gene names appear
        :param target_gene_col: the name of the column in the target file
        """
        assert isinstance(target_filename, (str, Path))
        assert isinstance(reference_filename, (str, Path))
        assert (isinstance(ref_wbgene_col, str))

        self.target_filename = target_filename
        self.reference_filename = reference_filename
        self.target_df = None
        self.not_identified = []

        self.gene_symbol_dict = dict()
        self.sequence_dict = dict()
        self.other_id_dict = dict()
        self.target_wbgene = None

        self.ref_wbgene_col = ref_wbgene_col
        self.ref_gene_col = ref_gene_col
        self.ref_seq_col = ref_seq_col
        self.ref_other_col = ref_other_col
        self.target_gene_col = target_gene_col

    def read_reference(self):
        """
        load the reference file into a pandas DataFrame
        """
        ref = general.load_csv(self.reference_filename, drop_gene_names=False)
        self.gene_symbol_dict = {name: wbgene for name, wbgene in zip(ref[self.ref_gene_col], ref[self.ref_wbgene_col])
                                 if not pd.isna(name)}
        self.sequence_dict = {name: wbgene for name, wbgene in zip(ref[self.ref_seq_col], ref[self.ref_wbgene_col])
                              if not pd.isna(name)}

        split_other_id = ref[self.ref_other_col].str.split(pat=";")
        for namelst, wbgene in zip(split_other_id, ref[self.ref_wbgene_col]):
            if isinstance(namelst, list):
                for name in namelst:
                    self.other_id_dict[name] = wbgene

    def import_target(self):
        """
        import the target filename into a pandas DataFrame
        """
        self.target_df = general.load_csv(self.target_filename, drop_gene_names=False)
        self.target_wbgene = [''] * self.target_df.shape[0]

    def generate_target_dict(self):
        """
        Generate the feature name to WBGene dictionary,
        and translate the imported pandas DataFrame feature names to WBGene.
        """
        for i, tgt_gene in enumerate(self.target_df[self.target_gene_col]):
            if tgt_gene in self.sequence_dict:
                self.target_wbgene[i] = self.sequence_dict[tgt_gene]
            elif tgt_gene in self.gene_symbol_dict:
                self.target_wbgene[i] = self.gene_symbol_dict[tgt_gene]
            elif tgt_gene in self.other_id_dict:
                self.target_wbgene[i] = self.other_id_dict[tgt_gene]
            else:
                self.not_identified.append((i, tgt_gene))

        self.target_df.insert(loc=0, column='WBGene', value=self.target_wbgene)
        self.target_df.set_index('WBGene', inplace=True)

    def rewrite_target(self):
        """
        save the target file as a csv with the suffix "_WBGene".
        """
        general.save_to_csv(self.target_df, self.target_filename, '_WBGene')

    def run(self):
        """
        load the specified file, translate its feature names to WBGenes, and save it as a csv with the "_WBGene" suffix.
        """
        self.read_reference()
        self.import_target()
        self.generate_target_dict()
        self.rewrite_target()


def run_all(path):
    """
   run the translator on a list of files
   """
    target_filenames = [path]
    ref_filename = __gene_names_and_biotype__
    ref_wbgene_col = 'gene'
    ref_gene_col = 'gene ID'
    ref_seq_col = 'Sequence ID'
    ref_other_col = 'Other IDs'
    target_gene_col = 'genes'

    tran = []
    for fname in target_filenames:
        tran.append(GeneNameTranslator(fname, ref_filename, ref_wbgene_col, ref_gene_col, ref_seq_col, ref_other_col,
                                       target_gene_col))
        tran[-1].run()
    return tran


def name_to_wbgene_folder(pth: str):
    """
    iterate over all files in a directory, translate their feature names to WBGene
    and save them as csv in the same folder.
    :param pth: path of the folder on which to iterate
    :return:
    a list of all GeneNameTranslator objects
    """
    tran = []
    mypth = Path(pth)
    for item in mypth.iterdir():
        if item.is_file():
            tran.append(GeneNameTranslator(item))
            tran[-1].run()
    return tran

# TODO: organize this!!!
