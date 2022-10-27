#! /usr/bin/env python

##############################################################################
##  pyvolve: Python platform for simulating evolutionary sequences.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) 
##############################################################################

'''
    This module contains various genetic definitions and mappings used throughout pyvolve.
'''
from collections import defaultdict

from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA


class Genetics():
    '''
        Molecular alphabet objects.
    '''
    
    def __init__(self, gencode=2):
        '''
            Set up internally-used genetic code lists.
        '''
        codontable: NCBICodonTableDNA = CodonTable.unambiguous_dna_by_id[gencode]

        self.pyrims       = ["C", "T"]
        self.purines      = ["A", "G"]
        self.nucleotides  = ["A", "C", "G", "T"]
        self.amino_acids  = list(codontable.protein_alphabet)
        self.genetic_code = self.__create_codon_table(codontable)
        self.codon_dict   = codontable.forward_table
        self.codons       = sorted(codontable.forward_table.keys())
        self.stop_codons  = codontable.stop_codons

        # PAML order of amino acids
        self.paml_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    @staticmethod
    def __create_codon_table(codontable: NCBICodonTableDNA):
        aa2codons = defaultdict(set)
        for cdn, aa in codontable.forward_table.items():
            aa2codons[aa].add(cdn)
        genetic_code = [sorted(aa2codons[aa]) for aa in codontable.protein_alphabet]
        return genetic_code

Genetics()
