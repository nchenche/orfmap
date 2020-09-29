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
    logger = logHandler.Logger(name='main', outpath=param.outpath)
    logo(logger)
    param.description()

    # parses fasta & gff by chromosomes
    logger.title('# Parsing GFF and fasta input files #')
    fasta_hash = fasta_parser.parse(fasta_filename=param.fasta_fname)
    gff_data = gff_parser.parse(param=param, fasta_hash=fasta_hash)

    # if param.write_types:
    #     gff_data.write_types(outpath=param.outpath)
    #     sys.exit(0)
    # elif param.show_types:
    #     gff_data.show_types(outpath=param.outpath)

    # checking if type(s) given in argument is(are) valid
    inspect.check_types(gff_data=gff_data, types=param.types)

    # ORFs mapping (scans genome for stop-to-stop sequences and assigns them a status)
    logger.title('# Mapping ORFs (stop-to-stop codons) #')
    all_orfs = orfmap.mapping(gff_data=gff_data, fasta_hash=fasta_hash, param=param)

    # writes outputs
    # logger.title('# Writing GFF and fasta output files #')
    # seqio.write_orfs(all_orfs=all_orfs, param=param)


def logo(logger):
    logger.info('')
    logger.info('   ___    ___   ___   __  __               ')
    logger.info('  / _ \  | _ \ | __| |  \/  |  __ _   _ __ ')
    logger.info(' | (_) | |   / | _|  | |\/| | / _` | | ')
    logger.info('  \___/  |_|_\ |_|   |_|  |_| \__,_| | .__/')
    logger.info('                                     |')


if __name__ == '__main__':
    sys.exit(main())

