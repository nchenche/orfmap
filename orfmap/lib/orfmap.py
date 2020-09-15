import os
import sys

from orfmap.lib import logHandler
from orfmap.lib import fasta_parser
from orfmap.lib import gff_parser
from orfmap.lib import parameters
from orfmap.lib import seqio
from orfmap.lib import inspect


logger = logHandler.get_logger(name='__name__')


def mapping(gff_data, fasta_hash, param):

    all_orfs = []
    for chr_id in sorted(fasta_hash):

        log = '# Mapping ORFs (stop-to-stop codons) for chromosome {} #'.format(chr_id)
        log_deco = '-' * len(log)
        logger.info(log_deco)
        logger.info(log)
        logger.info(log_deco)
        logger.info('')

        gff_chr = gff_data[chr_id]
        orfs = gff_parser.get_orfs(gff_chr=gff_chr, orf_len=param.orf_len)


        log = '# Assigning status for all {} ORFs found #'.format(len(orfs))
        log_deco = '-' * len(log)
        logger.info(log_deco)
        logger.info(log)
        logger.info(log_deco)
        logger.info('')

        for orf in sorted(orfs, key=lambda x: (x.seqid, x.start)):
            elements = gff_chr.get_elements(coors=orf.get_coors(), strand=orf.strand, types=param.types)
            orf.assignment(elements, co_ovp=param.co_ovp)
            all_orfs.append(orf)

    return all_orfs