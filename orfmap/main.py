# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import sys
from orfmap.lib import orfmap
from orfmap.lib import logHandler
from orfmap.lib import fasta_parser
from orfmap.lib import gff_parser
from orfmap.lib import parameters
from orfmap.lib import inspect
from orfmap.lib import tools
import os

def main():
    # gets arguments
    param = parameters.Param(args=parameters.get_args())

    if param.bool_chrs:
        tools.get_infos(_input=param.gff_fname, option='chrs')
        sys.exit(0)
    elif param.bool_types:
        tools.get_infos(_input=param.gff_fname, option='types')
        sys.exit(0)

    logger = logHandler.Logger(name='main', outpath=param.outpath)
    logo(logger)
    param.description()

    # parses fasta & gff by chromosomes
    logger.title('# Parsing GFF and fasta input files #')
    fasta_hash = fasta_parser.parse(fasta_filename=param.fasta_fname)
    gff_data = gff_parser.parse(param=param, fasta_hash=fasta_hash, chr_id=param.chr)

    # with open(param.outpath + os.path.basename(param.gff_fname)[:-4] + '.faa', 'w') as faa_file:
    #     for chr_id in sorted(gff_data):
    #         gff_chr = gff_data[chr_id]
    #         for fasta in gff_chr.proteins_fasta():
    #             faa_file.write(fasta)

    # for chr_id in sorted(gff_data)[:1]:
    #     proteins = gff_data[chr_id].group_cds()
    #
    #     for protein in sorted(proteins):
    #         logger.info('## {}'.format(protein))
    #         for cds in proteins[protein]:
    #             logger.info('  - cds: {} {}/{} {}'.format(cds.id_, cds.idx, cds.nb, cds.get_coors()))
    # sys.exit(0)

    # checking if type(s) given in argument is(are) valid
    # inspect.check_types(gff_data=gff_data, types=param.types)

    # ORFs mapping (scans genome for stop-to-stop sequences and assigns them a status)
    logger.title('# Mapping ORFs (stop-to-stop codons) #')
    orfmap.mapping(gff_data=gff_data, param=param)


def logo(logger):
    logger.info('')
    logger.info('   ___    ___   ___   __  __               ')
    logger.info('  / _ \  | _ \ | __| |  \/  |  __ _   _ __ ')
    logger.info(' | (_) | |   / | _|  | |\/| | / _` | | ')
    logger.info('  \___/  |_|_\ |_|   |_|  |_| \__,_| | .__/')
    logger.info('                                     |')


if __name__ == '__main__':
    sys.exit(main())
