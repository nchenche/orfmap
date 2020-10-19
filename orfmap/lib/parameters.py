# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:37:10 2020

@author: nicolas
"""

import os
import argparse
from orfmap.lib import logHandler

logger = logHandler.Logger(name=__name__)


class Param:
    """
    Wrapper of all input files and arguments
    """
    def __init__(self, args):
        self.fasta_fname = args.fna
        self.gff_fname = args.gff
        self.chr = args.chr
        if 'CDS' not in args.type:
            args.type.append('CDS')
        self.types = list(set(args.type))
        self.include = args.o_include
        self.exclude = args.o_exclude
        self.orf_len = args.orf_len
        self.co_ovp = args.co_ovp
        self.bool_types = args.bool_types
        self.bool_chrs = args.bool_chrs

        self.outpath = args.out + '/'
        os.makedirs(self.outpath, exist_ok=True)
        self.default_basename = 'mapping_orf_'
        self.default_mainname = '.'.join(os.path.basename(self.fasta_fname).split('.')[:-1])
        self.outfile = self.outpath + self.default_basename + self.default_mainname

    def description(self):
        """

        A formatted string describing all key index positions stored.

        Returns:
            object: str

        """
        chr = self.chr if self.chr else 'None'
        logger.info('')
        logger.info('Parameters description:')
        logger.info('- fasta filename: ' + self.fasta_fname)
        logger.info('- gff filename: ' + self.gff_fname)
        logger.info('- chr: ' + chr)
        logger.info('- types: ' + ', '.join(self.types))
        logger.info('- o_include: ' + ', '.join(self.include))
        logger.info('- o_exclude: ' + ', '.join(self.exclude))
        logger.info('- orf_len: ' + str(self.orf_len))
        logger.info('- co_ovp : ' + str(self.co_ovp))
        logger.info('- outfile: ' + self.outfile)
        logger.info('- bool_types: ' + str(self.bool_types))
        logger.info('- bool_chrs: ' + str(self.bool_chrs))
        logger.info('')


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
    parser.add_argument("-chr", required=False, nargs="?", type=str, default=None,
                        help="Chromosome name")
    parser.add_argument("-type", required=False, nargs="+", default=['CDS'],
                        help="Type feature(s) a flag is desired for ('CDS' in included by default).")
    parser.add_argument("-o_include", required=False, nargs="+", default=['all'],
                        help="Type feature(s) and/or Status attribute(s) desired to be written in the output (all by default).")
    parser.add_argument("-o_exclude", required=False, nargs="+", default=[],
                        help="Type feature(s) and/or Status attribute(s) desired to be excluded (None by default).")
    parser.add_argument("-orf_len", required=False, nargs="?", default=60, type=int,
                        help="Minimum number of nucleotides required to define a sequence between two consecutive stop codons\
                         as an ORF sequence (60 nucleotides by default).")
    parser.add_argument("-co_ovp", required=False, nargs="?", default=0.7, type=float,
                        help="Cutoff defining the minimum CDS overlapping ORF fraction required to label on ORF as a CDS.\
                             By default, an ORF sequence will be tagged as a CDS if at least 70 per cent of its sequence overlap\
                             with the CDS sequence.")
    parser.add_argument("-out", required=False, nargs="?", default='./', type=str,
                        help="Output directory")
    parser.add_argument('--show-types', action='store_true', default=False,
                        dest='bool_types',
                        help='Print all type features')
    parser.add_argument('--show-chrs', action='store_true', default=False,
                        dest='bool_chrs',
                        help='Print all chromosome names')


    args = parser.parse_args()
    
    return args
