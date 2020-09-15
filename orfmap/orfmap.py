# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import os
import sys

from orfmap.lib import logHandler
from orfmap.lib import fasta_parser
from orfmap.lib import gff_parser
from orfmap.lib import parameters
from orfmap.lib import seqio
from orfmap.lib import inspect


def main():
    # gets arguments
    param = parameters.Param(args=parameters.get_args())
    logger = logHandler.get_logger(name='orfmap', outpath=param.outpath)

    # parses fasta & gff by chromosomes
    logger.info('Parsing input files')
    fasta_hash = fasta_parser.parse(fasta_filename=param.fasta_fname)
    gff_data = gff_parser.parse(gff_fname=param.gff_fname, fasta_hash=fasta_hash)

    # checking if type(s) given in argument is(are) valid
    inspect.check_types(gff_data=gff_data, types=param.types)
    sys.exit(0)

    all_orfs = []
    for chr_id in sorted(fasta_hash):
        print('\n-----------------------')
        print('#Reading chromosome: {}'.format(chr_id))    
        if chr_id in (gff_data):
            gff_chr = gff_data[chr_id]
        else:
            print('*** Warning: Chromosome {} not in {}. It will not be considered. ***\n'.format(chr_id, os.path.basename(gff_fname)))
            pass
    
        # searching for ORFs        
        print('#Searching for ORFs')
        orfs = gff_parser.get_orfs(gff_chr=gff_chr, orf_len=param.orf_len)

        # assigning ORFs
        print('#Assigning ORFs status')
        for orf in sorted(orfs, key=lambda x: (x.seqid, x.start)):
            elements = gff_chr.get_elements(coors=orf.get_coors(), strand=orf.strand, types=param.types)
            orf.assignment(elements, co_ovp=param.co_ovp)
            all_orfs.append(orf)

    # writes outputs
    seqio.write_orfs(all_orfs=all_orfs, param=param)
    
if __name__ == '__main__':
    sys.exit(main())

