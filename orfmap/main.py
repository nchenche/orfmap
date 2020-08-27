# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:52:13 2020

@author: nicolas
"""

import os
import sys
sys.path.append('/home/nicolas/spyder_workspace/ORFMap/orfmap')

from lib import fasta_parser
from lib import gff_parser
from lib import tools

#==============================================================================
# 
#==============================================================================
outpath = '/home/nicolas/spyder_workspace/ORFMap/orfmap/data/output/'
filenames = ['GCF_000146045.2_R64_genomic',
             'Acupr',
             'Hjila',
             'Hmedi']

filename = filenames[-2]

fasta_fname = '/home/nicolas/spyder_workspace/ORFMap/orfmap/data/'+filename+'.fa'
gff_fname = '/home/nicolas/spyder_workspace/ORFMap/orfmap/data/'+filename+'.gff'
outfile = outpath+'mapping_orf_'+filename

#==============================================================================
# MAIN TMP
#==============================================================================
fasta_hash = fasta_parser.parse(fasta_filename=fasta_fname)
gff_data = gff_parser.parse(gff_fname=gff_fname)

all_orfs = []
orf_sequences = []
for chr_id in sorted(fasta_hash):
    print('#Dealing with: {}'.format(chr_id))
    gff_chr = gff_data[chr_id]

    fasta_chr = fasta_hash[chr_id]

    print('#Searching for ORFs and assignment...')
    orf_coors, orfs = tools.get_orf_coors(fasta_chr=fasta_chr, orf_len=60, gff_chr=gff_chr)
    
    for orf in sorted(orfs, key=lambda x: x.start):
        sequence = '>'+orf._id+'\n'
        sequence += fasta_chr.translate(start=orf.start, end=orf.end, strand=orf.strand, phase='.')
        sequence += '\n'
        orf_sequences.append(sequence)
        all_orfs.append(orf)

# writes gff
print('#Writing gff file of ORFs: {}'.format(outfile+'_v2.gff'))
with open(outfile+'_v2.gff', "w") as out_file:
    header = '# Input genomic fasta file: {}\n'.format(os.path.basename(fasta_fname))
    header += '# Input gff file: {}\n'.format(os.path.basename(gff_fname))
    out_file.write(header)
    for orf in sorted(all_orfs, key=lambda x: x.start):
        out_file.write('\t'.join(orf.gff_line)+'\n')

#writes fasta
print('#Writing fasta file of ORFs: {}'.format(outfile+'_v2.gff'))
with open(outfile+'_v2.fa', "w") as out_file:
    for sequence in orf_sequences:
        out_file.write(sequence)
           
sys.exit(1)

[ (x.start, x._id) for x in sorted(orfs[:10], key=lambda x: x.start) ]

#==============================================================================
# DEBUG
#==============================================================================
chr_id = sorted(fasta_hash)[-1]

gff_chr1 = gff_data[chr_id]
fasta_chr1 = fasta_hash[chr_id]


orf = [ x for x in all_orfs if x.start == 2372 ][0]
orf.frame

cds = [ x for x in gff_chr1.gff_elements if x.start == 2372 ][-1]
cds.frame

elements = gff_chr1.get_elements(coors=orf.get_coors(), strand='-', types=['CDS'])
[(x._id, x.get_coors(), x.frame, x.strand) for x in elements]



tools.assign_orf_type(start=2372, end=6937, strand='-', gff_chr=gff_chr1, frame=0)
#1    2368 
#2372    6835 
fasta_chr1.get_seq(132800, 132821, strand='+', phase=0)
fasta_chr1.translate(2372, 6835, strand='-', phase=0)
fasta_chr1.translate(2372, 6937, strand='-', phase='.')

fasta_chr1.get_seq(5530, 5544, strand='-', phase=0)
fasta_chr1.translate(5530, 5623,strand='-', phase='.')
fasta_chr1.translate(5530, 5622,strand='-', phase='.')

orfs_chr1 = [ x for x in all_orfs if chr_id in x.parent ]
orfs_chr1_f0 = [ x for x in orfs_chr1 if x.frame == 0 and x.strand == '+']
orfs_chr1_f1 = [ x for x in orfs_chr1 if x.frame == 1 and x.strand == '+']
orfs_chr1_f2 = [ x for x in orfs_chr1 if x.frame == 2 and x.strand == '+']

fasta_chr1.translate(1, orfs_chr1_f0[0].get_coors()[1], strand='+', phase='.')
fasta_chr1.translate(2, orfs_chr1_f1[0].get_coors()[1], strand='+', phase='.')
fasta_chr1.translate(3, orfs_chr1_f2[0].get_coors()[1], strand='+', phase='.')


orfs_chr1_f0_minus = [ x for x in orfs_chr1 if x.frame == 0 and x.strand == '-']
orfs_chr1_f1_minus = [ x for x in orfs_chr1 if x.frame == 1 and x.strand == '-']
orfs_chr1_f2_minus = [ x for x in orfs_chr1 if x.frame == 2 and x.strand == '-']


fasta_chr1.translate(orfs_chr1_f0_minus[-1].get_coors()[0], fasta_chr1.nucid_max-0, strand='-', phase='.')
fasta_chr1.translate(orfs_chr1_f1_minus[-1].get_coors()[0], fasta_chr1.nucid_max-1, strand='-', phase='.')
fasta_chr1.translate(orfs_chr1_f2_minus[-1].get_coors()[0], fasta_chr1.nucid_max-2, strand='-', phase='.')

fasta_chr1.translate(229565-6, 229816+3*6, strand='-', phase='.')


sequence = fasta_chr1.get_seq(orfs_chr1_f2_minus[-1].get_coors()[1], fasta_chr1.nucid_max, strand='+', phase=0)
frame = 0
codons = [ sequence[i:i+3].upper() for i in range(frame, len(sequence), 3) ]

#1024681	1025532
#1024543 1025532
tools.get_overlap([1024543, 1025532], [1024681, 1025532])
offset = (3 - ((1025532-1024681+1-0)%3) ) % 3

#==============================================================================
# CHECK
#==============================================================================
def read_fasta(fasta_fname=None):
    """
    - Input: fasta file
    - Output: dictionary with 1st header part as keys and sequences as values
    """
    fasta = {}
    with open(fasta_fname, 'r') as sequences_file:
        for l in sequences_file:
            if l.startswith('>'):
                header = l.strip()
            else:
                sequence = l.strip().replace('*', '')
                fasta[sequence] = header
            
    return fasta


fasta_orfmap = '/home/nicolas/spyder_workspace/ORFMap/orfmap/data/mapping_orf_Hjila.fa'
fasta_strali = '/home/nicolas/spyder_workspace/ORFMap/orfmap/data/output/control/Hjila_IGORF.pmultifasta'

seq_orfmap = read_fasta(fasta_orfmap)
seq_strali = read_fasta(fasta_strali)

len(sorted(seq_orfmap))
len(sorted(seq_strali))

uniq_strali = []
for key in sorted(seq_strali):
    if key not in seq_orfmap:
        uniq_strali.append(seq_strali[key])
        print(seq_strali[key])

len(uniq_strali)
val = sorted(uniq_strali)[0:10][-1]

[key for key in sorted(seq_strali) if seq_strali[key] == val]


