# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import os
import sys

#sys.path.append('/home/nicolas/spyder_workspace/ORFMap/orfmap')

from orfmap.lib import fasta_parser
from orfmap.lib import gff_parser
from orfmap.lib import parameters

def main():
    # gets arguments
    args = parameters.get_args()

    outpath = args.out
    fasta_fname = args.fna
    gff_fname = args.gff
    outfile = outpath+'mapping_orf_' + '.'.join(os.path.basename(fasta_fname).split('.')[:-1])

    # parses fasta & gff by chromosomes    
    fasta_hash = fasta_parser.parse(fasta_filename=fasta_fname)
    gff_data = gff_parser.parse(gff_fname=gff_fname, fasta_hash=fasta_hash)
    
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
        orfs = gff_parser.get_orfs(gff_chr=gff_chr, orf_len=60)

        # assigning ORFs
        print('#Assigning ORFs status')
        for orf in sorted(orfs, key=lambda x: (x.seqid, x.start)):
            elements = gff_chr.get_elements(coors=orf.get_coors(), strand=orf.strand, types=['CDS'])
            orf.assignment(elements)
            all_orfs.append(orf)

    print('\n-----------------------')
    print('#Writing output files (gff and fasta)')
    with open(outfile+'.gff', "w") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(gff_fname))
        out_gff.write(header)
        with open(outfile+'.fa', "w") as out_fasta:
            for orf in all_orfs:
                out_gff.write(orf.get_gffline())
                out_fasta.write(orf.get_fastaline())

    

if __name__ == '__main__':
    sys.exit(main())

