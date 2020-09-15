# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:37:10 2020

@author: nicolas
"""

import os
import argparse


class Param:
    """
    Wrapper of all input files and arguments
    """
    def __init__(self, args):
        self.fasta_fname = args.fna
        self.gff_fname = args.gff
        if 'CDS' not in args.type:
            args.type.append('CDS')
        self.types = list(set(args.type))
        self.orf_len = args.orf_len
        self.co_ovp = args.co_ovp

        self.outpath = args.out + '/'
        os.makedirs(self.outpath, exist_ok=True)
        self.default_basename = 'mapping_orf_'
        self.default_mainname = '.'.join(os.path.basename(self.fasta_fname).split('.')[:-1])
        self.outfile = self.outpath + self.default_basename + self.default_mainname


def get_args():
    """
    Returns:
        Parameters
    """

    parser = argparse.ArgumentParser(description='Genomic mapping of pseudo-ORF')
    parser.add_argument("-fna", required=True, nargs="?",
                        help="Genomic fasta file (.fna) ")
    parser.add_argument("-gff", required=True, nargs="?",
                        help="GFF annotation file (.gff)")
    parser.add_argument("-type", required=False, nargs="+", default=['CDS'],
                        help="Type attribute(s) a flag is desired for")
    parser.add_argument("-orf_len", required=False, nargs="?", default=60, type=int,
                        help="Minimum number of nucleotides required to define a sequence between two consecutive stop codons\
                         as an ORF sequence (60 nucleotides by default).")
    parser.add_argument("-co_ovp", required=False, nargs="?", default=0.7, type=float,
                        help="Cutoff defining the minimum CDS overlapping ORF fraction required to label on ORF as a CDS.\
                             By default, an ORF sequence will be tagged as a CDS if at least 70 per cent of its sequence overlap\
                             with the CDS sequence.")
    parser.add_argument("-out", required=False, nargs="?", default='./', type=str,
                        help="Output directory")

    args = parser.parse_args()
    
    return args
