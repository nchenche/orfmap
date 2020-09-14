# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:59:28 2020

@author: nicolas
"""
import sys
from orfmap.lib import inspect


class GffElement:

    def __init__(self, gff_line=None, fasta_chr=None):
        self.gff_line = gff_line.split() if gff_line else None
        self.fasta_chr = fasta_chr if fasta_chr else None
        self.len_chr = fasta_chr.nucid_max if fasta_chr else None

        self.seqid = self.gff_line[0] if gff_line else None
        self.source = self.gff_line[1] if gff_line else None
        self.type = self.gff_line[2] if gff_line else None
        self.start = int(self.gff_line[3]) if gff_line else None
        self.end = int(self.gff_line[4]) if gff_line else None
        self.score = self.gff_line[5] if gff_line else '.'
        self.strand = self.gff_line[6] if gff_line else None

        if gff_line:
            if self.gff_line[7] != '.':
                self.phase = int(self.gff_line[7])
            else:
                self.phase = self.gff_line[7]

            self.frame = self._get_frame()
            self._get_attributes()
        else:
            self.phase = '.'
            self.frame = None
            self.init_attributes()

        self.ovp_phased = []
        self.ovp_unphased = []
        self.suborfs = []

    def init_attributes(self):
        self._id = None
        self.name = None
        self.parent = None
        self.status = None
        self.color = None

    def _get_attributes(self):
        attributes = self.gff_line[8:][0]

        if 'ID' in attributes:
            self._id = self._parse_attributes(key='ID')
        else:
            self._id = '_'.join([self.type, str(self.start), str(self.end), self.strand])
        if 'Name' in attributes:
            self.name = self._parse_attributes(key='Name')
        else:
            self.name = self._id

        if 'Parent' in attributes:
            self.parent = self._parse_attributes(key='Parent')
        else:
            self.parent = None
        self.status = None
        self.color = None

    def _parse_attributes(self, key):
        attributes_col = self.gff_line[8:][0].split(';')
        attribute = [ x for x in attributes_col if key in x ][0]

        return attribute.split('=')[-1]

    def _get_frame(self):
        if isinstance(self.phase, int) or 'ORF' in self.type:
            if self.strand == '+':
                return (self.get_coors()[0] - 1) % 3
            else:
                return (self.len_chr - self.get_coors()[1]) % 3
        else:
            return None

    def get_coors(self):
        return self._coors_adjusted()

    def _coors_adjusted(self):
        start = self.start
        end = self.end
        if isinstance(self.phase, int):
            offset = (3 - ((self.end-self.start+1-self.phase)%3) ) % 3

            if self.strand == '+':
                start = self.start + self.phase
                end = self.end+offset if self.end+offset <= self.len_chr else self.len_chr
            elif self.strand == '-':
                start = self.start-offset if self.start-offset > 0 else self.start
                end = self.end-self.phase

        return (start, end)

    def get_len(self):
        return self.get_coors()[1] - self.get_coors()[0] + 1

    def sequence(self):
        phase = 0 if not isinstance(self.phase, int) else self.phase
        return self.fasta_chr.get_seq(start=self.start, end=self.end,
                                 strand=self.strand, phase=phase)

    def translate(self):
        return self.fasta_chr.translate(start=self.start, end=self.end,
                                 strand=self.strand, phase=self.phase)

    def get_fastaline(self):
        fastaline = '>'+self._id+'\n'+self.translate()+'\n'

        if self.suborfs:
            for suborf in self.suborfs:
                fastaline += suborf.get_fastaline()
        return fastaline

    def get_gffline(self):
        if self.gff_line and self.type not in ['ORF_nc_5-CDS', 'ORF_nc_3-CDS']:
            return '\t'.join(self.gff_line)
        else:
            gff_line = self.seqid
            gff_line += '\t'+self.source
            gff_line += '\t' + self.type
            gff_line += '\t' + str(self.start)
            gff_line += '\t' + str(self.end)
            gff_line += '\t' + self.score
            gff_line += '\t' + self.strand
            gff_line += '\t' + self.phase
            gff_line += '\tID=' + self.format_id()
            # gff_line += ';Parent=' + self.seqid+'_'+'1-'+str(self.len_chr)
            gff_line += ';Parent=' + self.parent
            gff_line += ';Status=' + self.status
            gff_line += ';color=' + self.color

            if self.ovp_phased:
                gff_line += ';Ovp_with=' + '|'.join([ x.format_id() for x in self.ovp_phased ])
            elif self.ovp_unphased:
                gff_line += ';Ovp_with=' + '|'.join([ x.format_id() for x in self.ovp_unphased])

            if self.suborfs:
                gff_line += '\n'
                for suborf in self.suborfs:
                    gff_line += suborf.get_gffline()
                # print(gff_line)

            return gff_line + '\n' if not self.suborfs else gff_line

    def assignment(self, elements, co_ovp):
        self.check_ovp(elements, co_ovp)
        self._set_attributes()

    def check_ovp(self, elements, co_ovp=0.7):
        orf_ovp_max = -1
        if elements:
            for element in elements:
                orf_ovp, element_ovp = get_overlap(orf_coors=self.get_coors(), other_coors=element.get_coors())
                if element_ovp == 1.0 or orf_ovp >= co_ovp:
                    if isinstance(element.phase, int) and element.frame == self.frame:
                        if element not in self.ovp_phased:
                            self.ovp_phased.append(element)
                    else:
                        if element not in self.ovp_unphased:
                            if orf_ovp > orf_ovp_max:
                                orf_ovp_max = orf_ovp
                                self.ovp_unphased.append(element)

    def _set_attributes(self):
        self._set_type()
        self._id = self.format_id()
        self._set_parent()
        # self.parent = self.seqid+'_'+'1-'+str(self.len_chr)
        self._set_color()
        self._set_status()

        if self.ovp_phased:
            self._get_suborfs()

    def _set_type(self):
        if self.ovp_phased:
            self.type = 'ORF_CDS'
        elif self.ovp_unphased:
            self.type = 'ORF_nc_ovp-'+self.ovp_unphased[-1].type
        else:
            self.type = 'ORF_nc_intergenic'

    def format_id(self):
        return '_'.join([self.seqid, self.strand,
                         str(self.start)+'-'+str(self.end),
                         str(self.frame), self.type])

    def _set_parent(self):
        if self.ovp_phased:
            self.parent = '|'.join([x._id for x in self.ovp_phased])
        elif self.ovp_unphased:
            self.parent = '|'.join([x._id for x in self.ovp_unphased])
        else:
            self.parent = self.seqid+'_'+'1-'+str(self.len_chr)

    def _set_color(self):
        if 'ORF' in self.type:
            if 'nc' not in self.type:
                self.color = '#ff0000'  # ff4d4d
            else:
                if 'intergenic' in self.type:
                    self.color = '#3366ff'
                else:
                    self.color = '#2eb82e'
        else:
            print('Warning: this function is only made for ORF.\n')

    def _set_status(self):
        if self.ovp_phased:
            self.status = 'coding'
        else:
            self.status = 'non-coding'

    def _get_suborfs(self):
        # if len(self.ovp_phased) > 1:
        #     print(self.get_coors(), [ (x.get_coors(), x.parent) for x in self.ovp_phased ])
        gff_line = self.get_gffline()
        for element in self.ovp_phased:
            if self.strand == '+':
                if element.get_coors()[0]-1 - self.start + 1 >= 60:
                    suborf = GffElement(gff_line=gff_line, fasta_chr=self.fasta_chr)
                    suborf.end = element.get_coors()[0] - 1
                    suborf.type = 'ORF_nc_5-CDS'
                    suborf.status = 'non-coding'
                    suborf.color = '#FFFF00' #D67229
                    suborf._id = suborf.format_id()
                    self.suborfs.append(suborf)
                if self.end - element.get_coors()[1]+1 + 1 >= 60:
                    suborf = GffElement(gff_line=gff_line, fasta_chr=self.fasta_chr)
                    suborf.start = element.get_coors()[1] + 1
                    suborf.type = 'ORF_nc_3-CDS'
                    suborf.status = 'non-coding'
                    suborf.color = '#FFFF00'
                    suborf._id = suborf.format_id()
                    self.suborfs.append(suborf)
            else:
                if self.end - element.get_coors()[1]+1 + 1 >= 60:
                    suborf = GffElement(gff_line=gff_line, fasta_chr=self.fasta_chr)
                    suborf.start = element.get_coors()[1] + 1
                    suborf.type = 'ORF_nc_5-CDS'
                    suborf.status = 'non-coding'
                    suborf.color = '#FFFF00'
                    suborf._id = suborf.format_id()
                    self.suborfs.append(suborf)
                if element.get_coors()[0]-1 - self.start + 1 >= 60:
                    suborf = GffElement(gff_line=gff_line, fasta_chr=self.fasta_chr)
                    suborf.end = element.get_coors()[0] - 1
                    suborf.type = 'ORF_nc_3-CDS'
                    suborf.status = 'non-coding'
                    suborf.color = '#FFFF00'
                    suborf._id = suborf.format_id()
                    self.suborfs.append(suborf)


class Chromosome():

    def __init__(self, _id, fasta_chr):
        self._id = _id
        self.fasta_chr = fasta_chr
        self.start = 1
        self.end = fasta_chr.nucid_max
        self.source = None
        self.coors_intervals = self._set_intervals()
        self.gff_elements = []
        
    def _set_intervals(self, value=10000):
        return { (x, x + value - 1): [] for x in range(1, self.end, value) }
        
    def _get_intervals(self, coors):
        """
        Returns a list of coordinates that overlap with the element coordinates.
        """
        return [ x for x in self.coors_intervals if get_overlap(x, coors)[1] ]
        
    def _add_to_intervals(self):
        last_element = self.gff_elements[-1]
        for interval in self._get_intervals(coors=last_element.get_coors()):
            self.coors_intervals[interval].append(self.gff_elements.index(last_element))
        
    def add(self, gff_element):
        """
        Adds a Gff_element instance in self.gff_elements.
        
        Arguments:
            - gff_element: instance of Gff_element()
        """
        
        self.gff_elements.append(gff_element)
        self._add_to_intervals()
        
    def get_elements_in_intervals(self, coors):
        intervals_mx = [ self.coors_intervals[x] for x in self._get_intervals(coors=coors) ]
        intervals_flat = set([val for sublist in intervals_mx for val in sublist])

        return [ self.gff_elements[x] for x in intervals_flat ]
        
    def get_elements(self, coors=None, frame=None, strand='+', types=None):
        """
        Returns a list of Gff_element instances of CDS type. If the frame is
        given, only CDS in this frame will be returned, all CDS otherwise.
        
        Arguments:
            - frame: None or 0, 1 or 2
            
        Returns:
            - list of Gff_element instances        
        """
        if types:
            if coors:
                elements = [ x for x in self.get_elements_in_intervals(coors) if x.type in types and x.strand == strand ]
            else:
                elements = [ x for x in self.gff_elements if x.type in types and x.strand == strand]
        else:
            if coors:
                elements = [ x for x in self.get_elements_in_intervals(coors) if x.strand == strand ]
            else:
                elements = [ x for x in self.gff_elements if x.strand == strand]

        if frame == None:
            return elements
        else:
            return [ x for x in elements if x.frame == frame ]
        
    def get_types(self):
        return set([ x.type for x in self.gff_elements ])
        
    def sequence(self):
        return self.fasta_chr.get_seq(start=self.start, end=self.end,
                                      strand='+', phase=0)
                                      
    def rev_comp(self):
        return self.fasta_chr.rev_comp(self.sequence())


def get_overlap(orf_coors=(), other_coors=None):
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


def get_orfs(gff_chr, orf_len=60):
    orfs = []              
    sequence = gff_chr.sequence()
    pos = 0

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
                    orf = GffElement(fasta_chr=gff_chr.fasta_chr)
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
                        orf = GffElement(fasta_chr=gff_chr.fasta_chr)
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
                orf = GffElement(fasta_chr=gff_chr.fasta_chr)
                orf.seqid = gff_chr._id
                orf.source = gff_chr.source
                orf.strand = '-'
                orf.frame = frame_rev
                orf.start = start_pos_rev
                orf.end = end_pos_rev
                orfs.append(orf)
                
    return orfs


GFF_DESCR = {}


def set_gff_descr(gff_fname):
    global GFF_DESCR
    GFF_DESCR = {}
    
    with open(gff_fname, 'r') as gff_file:
        line = gff_file.readline()
        while line:
            if not line.startswith('#'):
                name = line.split()[0]
                pos_chr = gff_file.tell()-len(line)
                if name not in GFF_DESCR:
                    GFF_DESCR[name] = pos_chr
                
            line = gff_file.readline()


def parse(gff_fname, fasta_hash, chr_id=None):
    if not GFF_DESCR:
        set_gff_descr(gff_fname)

    chrs_common = inspect.check_chrids(chrs_gff=sorted(GFF_DESCR), chrs_fasta=sorted(fasta_hash))
    chr_ids = sorted(chrs_common) if not chr_id else [chr_id]
    if chr_id and chr_id not in chrs_common:
        print('Error: wrong chromosome id\n')
        sys.exit(1)
    # sys.exit(0)

    gff_data = {}
    with open(gff_fname, 'r') as gff_file:
        for chr_id in chr_ids:
            gff_file.seek(GFF_DESCR[chr_id], 0)
            line = gff_file.readline()
            chr_name = line.split()[0]
            while chr_name == chr_id:
                if chr_name not in gff_data:
                    gff_data[chr_name] = Chromosome(_id=chr_name, fasta_chr=fasta_hash[chr_id])
                    chromosome = gff_data[chr_name]
                    chromosome.source = line.split()[1]

                chromosome.add(gff_element=GffElement(gff_line=line, fasta_chr=fasta_hash[chr_id]))

                line = gff_file.readline()
                chr_name = line.split()[0]
            
    return gff_data

