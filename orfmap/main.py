# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import os
import sys

from orfmap.lib import orfmap
from orfmap.lib import logHandler
from orfmap.lib import fasta_parser
from orfmap.lib import gff_parser
from orfmap.lib import parameters
from orfmap.lib import seqio
from orfmap.lib import inspect


def main():
    # gets arguments
    param = parameters.Param(args=parameters.get_args())
    logger = logHandler.Logger(name='orfmap', outpath=param.outpath)


    # parses fasta & gff by chromosomes
    logger.info('# Parsing GFF and fasta input files #', decoration=True)
    fasta_hash = fasta_parser.parse(fasta_filename=param.fasta_fname)
    gff_data = gff_parser.parse(gff_fname=param.gff_fname, fasta_hash=fasta_hash)
    inspect.check_types(gff_data=gff_data, types=param.types) # checking if type(s) given in argument is(are) valid

    # ORFs mapping
    logger.info('# Mapping ORFs (stop-to-stop codons) #', decoration=True)
    all_orfs = orfmap.mapping(gff_data=gff_data, fasta_hash=fasta_hash, param=param)

    # writes outputs
    logger.info('# Writing GFF and fasta output files #', decoration=True)
    seqio.write_orfs(all_orfs=all_orfs, param=param)
    
if __name__ == '__main__':
    sys.exit(main())

