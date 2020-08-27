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
gff_data = gff_parser.parse(gff_fname=gff_fname, fasta_hash=fasta_hash)

all_orfs = []
orf_sequences = []
for chr_id in sorted(fasta_hash):
    print('#Dealing with: {}'.format(chr_id))
    gff_chr = gff_data[chr_id]

    print('#Searching for ORFs and assignment...')
    orfs = tools.get_orf_coors(gff_chr=gff_chr, orf_len=60)
    
    for orf in sorted(orfs, key=lambda x: x.start):
        sequence = '>'+orf._id+'\n'
        sequence += orf.translate(start=orf.start, end=orf.end, strand=orf.strand, phase=orf.phase)
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

gff_chr = gff_data[chr_id]

#==============================================================================
# BUILDING
#==============================================================================


def get_orfs(gff_chr, orf_len=60):
    orfs = []              
    sequence = gff_chr.sequence()

    # loops on each possible frame (the negative frame is defined in "frame_rev")
    for frame in range(3):
        # list of codons in frame "frame"
        codons = [ sequence[i:i+3].upper() for i in range(frame, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    
        start_pos = frame + 1
        end_pos = None
        
        frame_rev = (gff_chr.end%3 - frame) % 3
        start_pos_rev = None
        end_pos_rev = None
       
        for pos, codon in enumerate(codons):
            if codon in ['TAG', 'TGA', 'TAA']:
                end_pos = pos*3 + 1 + 2 + frame
                if end_pos - start_pos + 1 >= orf_len:
                    orf = gff_parser.GffElement()
                    orf.fasta_chr = gff_chr.fasta_chr
                    orf.seqid = gff_chr._id
                    orf.source = gff_chr.source
                    orf.strand = '+'
                    orf.frame = frame
                    if start_pos == frame + 1:
                        orf.start = start_pos
                    else:
                        orf.start = start_pos + 3
                    orf.end = end_pos
                    orfs.append(orf)

                start_pos = end_pos - 2
                    
            elif codon in ['CTA', 'TCA', 'TTA']:
                if start_pos_rev is None:
                    start_pos_rev = pos*3 + 1 + frame
                else:
                    end_pos_rev = pos*3 + 1 + 2 + frame
                    if end_pos_rev - start_pos_rev + 1 >= orf_len:
                        orf = gff_parser.GffElement(fasta_chr=gff_chr.fasta_chr)
                        orf.fasta_chr = gff_chr.fasta_chr
                        orf.seqid = gff_chr._id
                        orf.source = gff_chr.source
                        orf.strand = '-'
                        orf.frame = frame_rev
                        orf.start = start_pos_rev
                        orf.end = end_pos_rev - 3
                        orfs.append(orf)
                        
                    start_pos_rev = end_pos_rev - 2
        
        # adds coordinates of ORF in extremities
        if end_pos_rev:
            start_pos_rev = end_pos_rev - 2
            end_pos_rev = pos*3 + 1 + 2 + frame
            if end_pos_rev - start_pos_rev + 1 >= orf_len:
                orf = gff_parser.GffElement(fasta_chr=gff_chr.fasta_chr)
                orf.fasta_chr = gff_chr.fasta_chr
                orf.seqid = gff_chr._id
                orf.source = gff_chr.source
                orf.strand = '-'
                orf.frame = frame_rev
                orf.start = start_pos_rev
                orf.end = end_pos_rev
                orfs.append(orf)
                
    return orfs

def assign_orf_type(orfs, gff_chr):
    for orf in orfs:        
        elements = gff_chr.get_elements(coors=orf.get_coors(), strand=orf.strand, types=['CDS'])
        assignment(orf, elements)
        
def assignment(orf, elements):
    orf_ovp_max = -1
    if elements:
        for element in elements:
            orf_ovp, element_ovp = tools.get_overlap(orf_coors=orf.get_coors(), other_coors=element.get_coors())
            if element_ovp == 1.0 or orf_ovp >= 0.7:
                if isinstance(element.phase, int) and element.frame == orf.frame:
                    if element not in orf.ovp_phased:
                        orf.ovp_phased.append(element)
                else:
                    if element not in orf.ovp_unphased:
                        if orf_ovp > orf_ovp_max:
                            orf_ovp_max = orf_ovp
                            orf.ovp_unphased.append(element)
                
    orf._set_attributes()


#==============================================================================
chr_id = sorted(fasta_hash)[-1]
gff_chr = gff_data[chr_id]


orfs = get_orfs(gff_chr=gff_chr, orf_len=60)
len(orfs)
assign_orf_type(orfs, gff_chr)
orf = orfs[0]
len([x.type for x in orfs if x.type == 'ORF_CDS'])
len([x.type for x in orfs if x.type == 'ORF_nc_ovp-CDS'])
len([x.type for x in orfs if x.type == 'ORF_nc_intergenic'])

[ x.get_coors() for x in sorted(orfs[0:10], key=lambda x: (x.seqid, x.start)) ]

for orf in sorted(orfs, key=lambda x: (x.seqid, x.start)):
    if orf.type == 'ORF_CDS':
        break
    print(orf.get_gffline())
    
element = orf.ovp_phased[0]
tmp = [element, element]
orf.start, element.start


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
                if 'nc_' in header:
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


#==============================================================================
# 
#==============================================================================
nuc_comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def rev_comp(sequence):
    """
    Returns the reverse completary sequence of a given nucleotide sequence
    
    Arguments:
        - sequence (str): nucleotide sequence
        
    Returns:
        - reverse_complement (str)        
    """
    
    seq = list(sequence)
    seq.reverse()
    reverse_complement = ''
    for nuc in seq:             
        reverse_complement += nuc_comp[nuc.upper()]
        
    return reverse_complement        


sequence = 'TTTACGCCCGGGCTAGGGCCCTTTCCCCTAGGGTTTGGGCATCCCTTTAAACTAGGG'
sequence = 'TTTCCCCTAGGGTTTGGGCATCCCTTTAAACTAGGG'
sequence = sequence[21:29+1]
'-'.join([ sequence[i:i+3] for i in range(0, len(sequence), 3) ])

seq_rev = rev_comp(sequence)
codons = [ seq_rev[i:i+3] for i in range(0, len(seq_rev), 3) ]
'-'.join(codons)
protein_sequence = ''.join([ fasta_parser.Genecode[x.upper()] for x in codons if len(x) == 3 ])
'FKGPFG'

























