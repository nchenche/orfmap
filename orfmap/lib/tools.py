# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 18:27:31 2020

@author: nicolas
"""

#import os
from orfmap.lib import gff_parser

def set_orf_color(orf_type):
    if 'nc' not in orf_type:
        color = '#ff4d4d'
    else:
        if 'intergenic' in orf_type:
            color = '#3366ff'
        else:
            color = '#2eb82e'
#            color = 'ffe01a'
            
    return color

def set_gffline(gff_chr, start, end, strand, frame):
    gff_line = gff_chr._id
    gff_line += '\tRefSeq'
    orf_type = assign_orf_type(gff_chr, start, end, strand, frame)
    gff_line += '\t' + orf_type
    gff_line += '\t' + str(start)
    gff_line += '\t' + str(end)
    gff_line += '\t' + '.'
    gff_line += '\t' + strand
    gff_line += '\t' + '.'
    
    # attributes
    gff_line += '\tID=' + '_'.join([gff_chr._id, strand, str(start)+'-'+str(end), str(frame), orf_type])
    gff_line += ';Parent=' + gff_chr._id+'_'+'1-'+str(gff_chr.fasta_chr.nucid_max)
    gff_line += ';Status=' + orf_type.split('_')[-1]
    gff_line += ';color=' + set_orf_color(orf_type)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    
    return gff_line

def assign_orf_type(gff_chr, start, end, strand, frame):
    coors = (start, end)
    orf_type = None
    orf_status = []
    
    elements = gff_chr.get_elements(coors=coors, strand=strand, types=['CDS'])
    if elements:
        for element in elements:
            orf_ovp, element_ovp = get_overlap(orf_coors=coors, other_coors=element.get_coors())
            if element_ovp == 1.0 or orf_ovp >= 0.7:
                if isinstance(element.phase, int) and element.frame == frame:
                        orf_type = 'ORF_'+element.type
                else:
                    orf_type = 'ORF_nc_ovp-'+element.type
            else:
                orf_type = 'ORF_nc_intergenic'            
            
            orf_status.append((orf_type, orf_ovp, start, end, strand, frame,
                               element.get_coors()))
                
        orf_types = set([ x[0] for x in orf_status ])
        if 'ORF_CDS' in orf_types:
            orf_type = 'ORF_CDS'
#            if len([ x[0] for x in orf_status if x[0] == 'ORF_CDS']) > 1:
#                print('\n'+10*'*')
#                for orf in [ (x[0], x[1:]) for x in orf_status if x[0] == 'ORF_CDS']:
#                    print('#- {}'.format(orf))
        elif len(orf_types) == 1:
            orf_type = list(orf_types)[0]
        else:
            orf_type = sorted(orf_status, key = lambda x: x[1])[-1][0]
            
#        orf_types_test = [ x[0] for x in orf_status if x[0] != 'ORF_nc_intergenic' ]
#        if len(orf_types_test) > 1:
#            print(start, end, strand, frame, orf_types_test)
            
    else: 
        orf_type = 'ORF_nc_intergenic'
        
    return orf_type        


def get_orf_coors(gff_chr, orf_len=60):
    """
    Gets coordinates for ORFs in all possible frames. An ORF is defined from
    stop codon to stop codon.
    
    Argument:
        - orf_len: int value defining the minimum length for an ORF to be accepted
    
    Returns:
        - coors_stop{'frame_0': [[int, int], ...],
                 'frame_1': [[int, int], ...],
                 'frame_2': [[int, int], ...],
                 'frame_rev_0': [[int, int], ...],
                 'frame_rev_1': [[int, int], ...],
                 'frame_rev_2': [[int, int], ...]}
    
    """
    
    # gets nucleotide sequence of a chromosome
    sequence = gff_chr.sequence()
#    coors_stop = {}
    orfs = []
    
    # loops on each possible frame (the negative frame is defined in "frame_rev")
    for frame in range(3):
        # list of codons in frame "frame"
        codons = [ sequence[i:i+3].upper() for i in range(frame, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    
        start_pos = frame + 1
        end_pos = None
        
        frame_rev = (gff_chr.len_chr%3 - frame) % 3
        start_pos_rev = None
        end_pos_rev = None
       
        for pos, codon in enumerate(codons):
            if codon in ['TAG', 'TGA', 'TAA']:
                end_pos = pos*3 + 1 + 2 + frame
                if end_pos - start_pos + 1 >= orf_len:
                    if start_pos == frame + 1:
                        gff_line = set_gffline(gff_chr=gff_chr, start=start_pos, end=end_pos, strand='+', frame=frame)
                    else:
                        gff_line = set_gffline(gff_chr=gff_chr, start=start_pos+3, end=end_pos, strand='+', frame=frame)
                    orfs.append(gff_parser.GffElement(gff_line, gff_chr.fasta_chr))
                    
                start_pos = end_pos - 2
                    
            elif codon in ['CTA', 'TCA', 'TTA']:
                if start_pos_rev is None:
                    start_pos_rev = pos*3 + 1 + frame
                else:
                    end_pos_rev = pos*3 + 1 + 2 + frame
                    if end_pos_rev - start_pos_rev + 1 >= orf_len:
                        gff_line = set_gffline(gff_chr=gff_chr, start=start_pos_rev, end=end_pos_rev-3, strand='-', frame=frame_rev)
                        orfs.append(gff_parser.GffElement(gff_line, gff_chr.fasta_chr))
                        
                    start_pos_rev = end_pos_rev - 2
        
        # adds coordinates of ORF in extremities
        if end_pos_rev:
            start_pos_rev = end_pos_rev - 2
            end_pos_rev = pos*3 + 1 + 2 + frame
            if end_pos_rev - start_pos_rev + 1 >= orf_len:
                gff_line = set_gffline(gff_chr=gff_chr, start=start_pos_rev, end=end_pos_rev, strand='-', frame=frame_rev)
                orfs.append(gff_parser.GffElement(gff_line, gff_chr.fasta_chr))
        
    return orfs


def get_overlap(orf_coors=[], other_coors=[]):
    """
    Function defining if ORF coordinates overlap with another genomic element
    coordinates.
    
    Arguments:
        - orf_coors: start and end coordinates of an ncORF (list)
        - other_coors: start and end coordinates of a genomic element (list)
        
    Returns:
        - a tuple of the overlapping fraction between the ORF and the other element
        (float if overlap, None otherwise)
    """
    orf_ovp, other_ovp = 0, 0
    x_max = max(orf_coors[0], other_coors[0])
    y_min = min(orf_coors[1], other_coors[1])
    if x_max < y_min:
        len_ovp = y_min - x_max + 1
        len_orf = orf_coors[1] - orf_coors[0] + 1
        len_other = other_coors[1] - other_coors[0] + 1
        orf_ovp = round(len_ovp / float(len_orf), 2)
        other_ovp = round(len_ovp / float(len_other), 2)
        
    return (orf_ovp, other_ovp)
