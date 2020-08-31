# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:37:10 2020

@author: nicolas
"""

import argparse

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
    parser.add_argument("-out", required=False, nargs="?", default='./', type=str,
                        help="Output directory")
    
    return parser.parse_args()
