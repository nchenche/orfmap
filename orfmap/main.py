# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import os
import sys

from orfmap import orfmapper
from orfmap.lib import logHandler
from orfmap.lib import fasta_parser
from orfmap.lib import gff_parser
from orfmap.lib import parameters
from orfmap.lib import seqio
from orfmap.lib import inspect


def main():
    # gets arguments
    param = parameters.Param(args=parameters.get_args())
    logger = logHandler.get_logger(name='__name__', outpath=param.outpath)

    # parses fasta & gff by chromosomes
    logger.info('Parsing input files')
    fasta_hash = fasta_parser.parse(fasta_filename=param.fasta_fname)
    gff_data = gff_parser.parse(gff_fname=param.gff_fname, fasta_hash=fasta_hash)

    # checking if type(s) given in argument is(are) valid
    inspect.check_types(gff_data=gff_data, types=param.types)
    sys.exit(0)

    all_orfs = orfmap.mapping(gff_data=gff_data, fasta_hash=fasta_hash, param=param)

    # writes outputs
    seqio.write_orfs(all_orfs=all_orfs, param=param)
    
if __name__ == '__main__':
    sys.exit(main())

