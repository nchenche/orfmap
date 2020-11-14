# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:59:28 2020

@author: nicolas
"""
import sys
from orfmap.lib import logHandler
from orfmap.lib import inspect


logger = logHandler.Logger(name=__name__)


class GffElement:

    def __init__(self, gff_line=None, fasta_chr=None):
        self.gff_line = gff_line.split('\t') if gff_line else None
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
            self._id = None
            self.name = None
            self.parent = None
            self.status = None
            self.color = None

        self.ovp_phased = []
        self.ovp_unphased = []
        self.suborfs = []

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
        attribute = [x for x in attributes_col if key in x][0]

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
            offset = (3 - ((self.end-self.start+1-self.phase) % 3)) % 3

            if self.strand == '+':
                start = self.start + self.phase
                end = self.end+offset if self.end+offset <= self.len_chr else self.len_chr
            elif self.strand == '-':
                start = self.start-offset if self.start-offset > 0 else self.start
                end = self.end-self.phase

        return start, end

    def get_len(self):
        return self.get_coors()[1] - self.get_coors()[0] + 1

    def sequence(self):
        phase = 0 if not isinstance(self.phase, int) else self.phase

        return self.fasta_chr.sequence(start=self.start, end=self.end, strand=self.strand, phase=phase)

    def translate(self):
        return self.fasta_chr.translate(start=self.start, end=self.end, strand=self.strand, phase=self.phase)

    def get_fastaline(self):
        fastaline = '>'+self._id+'\n'+self.translate()+'\n'

        return fastaline

    def get_gffline(self):
        if self.gff_line and self.type not in ['nc_5-CDS', 'nc_3-CDS']:
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
            gff_line += ';Parent=' + self.parent
            gff_line += ';Status=' + self.status
            gff_line += ';color=' + self.color

            if self.ovp_phased:
                gff_line += ';Ovp_with=' + '|'.join([x.name for x in self.ovp_phased])
            elif self.ovp_unphased:
                gff_line += ';Ovp_with=' + '|'.join([x.name for x in self.ovp_unphased])

            return gff_line + '\n'

    def set_type(self):
        if self.ovp_phased:
            self.type = 'c_CDS'
        elif self.ovp_unphased:
            if 'CDS' in [x.type for x in self.ovp_unphased]:
                ovp_unphased_CDS = [x for x in self.ovp_unphased if x.type == 'CDS']

                ovp_elements_same_strand = [x for x in ovp_unphased_CDS if x.strand == self.strand]
                if ovp_elements_same_strand:
                    highest_overlapping_element = ovp_elements_same_strand[-1]
                else:
                    highest_overlapping_element = ovp_unphased_CDS[-1]
            else:
                ovp_elements_same_strand = [x for x in self.ovp_unphased if x.strand == self.strand]
                if ovp_elements_same_strand:
                    highest_overlapping_element = ovp_elements_same_strand[-1]
                else:
                    highest_overlapping_element = self.ovp_unphased[-1]

            self.type = 'nc_ovp-' + highest_overlapping_element.type
            if highest_overlapping_element.strand != self.strand:
                self.type += '-opp'
        else:
            self.type = 'nc_intergenic'

    def format_id(self):
        return '_'.join([self.seqid, self.strand,
                         str(self.start)+'-'+str(self.end),
                         str(self.frame), self.type])

    def set_parent(self):
        self.parent = self.seqid+'_'+'1-'+str(self.len_chr)

    def set_color(self):
        if self.type == 'c_CDS':
            self.color = '#ff4d00'  # ff4d4d
        else:
            if self.type == 'nc_intergenic':
                self.color = '#005073'
            else:
                if 'opp' in self.type:
                    self.color = '#71c7ec'
                else:
                    self.color = '#2935d6'

    def set_status(self):
        if self.ovp_phased:
            self.status = 'coding'
        else:
            self.status = 'non-coding'


class Chromosome:

    def __init__(self, _id, fasta_chr):
        self._id = _id
        self.fasta_chr = fasta_chr
        self.start = 1
        self.end = fasta_chr.nucid_max
        self.source = None
        self.coors_intervals = self._set_intervals()
        self.gff_elements = []
        
    def _set_intervals(self, value=50000):
        return {(x, x + value - 1): [] for x in range(1, self.end, value)}
        
    def _get_intervals(self, coors):
        """
        Returns a list of coordinates that overlap with the element coordinates.
        """
        return (x for x in self.coors_intervals if get_overlap(x, coors)[2])
        
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
        intervals_mx = (self.coors_intervals[x] for x in self._get_intervals(coors=coors))
        intervals_flat = sorted(set([val for sublist in intervals_mx for val in sublist]))

        return (self.gff_elements[x] for x in intervals_flat)
        
    def get_elements(self, coors=None, frame=None, strand=None, types=None):
        """
        Returns a list of Gff_element instances of a given type. If the frame is
        given, only CDS in this frame will be returned, all CDS otherwise.
        
        Arguments:
            - frame: None or int in 0, 1 or 2
            - strand: None or str in '+' or '-'
            - coors: None or list of coors in the form [(104, 395), ...]
            - types: None or list of str (e.g. ['CDS', 'tRNA']
            
        Returns:
            - list of Gff_element instances        
        """

        if types:
            if coors:
                if strand:
                    elements = (x for x in self.get_elements_in_intervals(coors) if x.type in types and x.strand == strand)
                else:
                    elements = (x for x in self.get_elements_in_intervals(coors) if x.type in types)
            else:
                if strand:
                    elements = (x for x in self.gff_elements if x.type in types and x.strand == strand)
                else:
                    elements = (x for x in self.gff_elements if x.type in types)
        else:
            if coors:
                if strand:
                    elements = (x for x in self.get_elements_in_intervals(coors) if x.strand == strand)
                else:
                    elements = (x for x in self.get_elements_in_intervals(coors))
            else:
                if strand:
                    elements = (x for x in self.gff_elements if x.strand == strand)
                else:
                    elements = (x for x in self.gff_elements)

        if not frame:
            return elements
        else:
            return (x for x in elements if x.frame == frame)
        
    def get_types(self):
        return set([x.type for x in self.gff_elements])
        
    def sequence(self, start: int, end: int):
        start = start if start else self.start
        end = end if end else self.end
        return self.fasta_chr.sequence(start=start, end=end, strand='+', phase=0)
                                      
    def rev_comp(self):
        return self.fasta_chr.reverse_complement(self.sequence())


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
    is_overlap = False
    orf_ovp, other_ovp = 0, 0
    x_max = max(orf_coors[0], other_coors[0])
    y_min = min(orf_coors[1], other_coors[1])
    if x_max < y_min:
        is_overlap = True
        len_ovp = y_min - x_max + 1
        len_orf = orf_coors[1] - orf_coors[0] + 1
        len_other = other_coors[1] - other_coors[0] + 1
        orf_ovp = len_ovp / float(len_orf)
        other_ovp = len_ovp / float(len_other)
        
    return orf_ovp, other_ovp, is_overlap


def get_orfs(gff_chr, orf_len=60):
    orfs = []
    sequence = gff_chr.sequence()
    pos = 0

    # loops on each possible frame (the negative frame is defined in "frame_rev")
    for frame in range(3):
        # list of codons in frame "frame"
        codons = [sequence[i:i+3].upper() for i in range(frame, len(sequence), 3) if len(sequence[i:i+3]) == 3]

        start_pos = frame + 1

        frame_rev = (gff_chr.end % 3 - frame) % 3
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


def is_orf_asked(orf=None, param=None):
    if 'all' in param.include:
        if not param.exclude:
            return True
        else:
            if not is_orf_exclude(orf=orf, exclude=param.exclude):
                return True
            else:
                return False
    else:
        if is_orf_include(orf=orf, include=param.include):
            if not is_orf_exclude(orf=orf, exclude=param.exclude):
                return True
            else:
                return False
        else:
            return False


def is_orf_include(orf=None, include=None):
    return True in [x in [orf.type, orf.status] for x in include]


def is_orf_exclude(orf=None, exclude=None):
    return True in [x in [orf.type, orf.status] for x in exclude]


GFF_DESCR = {}


def set_gff_descr(gff_fname):
    global GFF_DESCR
    GFF_DESCR = {}
    
    with open(gff_fname, 'r') as gff_file:
        line = gff_file.readline()
        while line:
            if not line.startswith('#'):
                name = line.split('\t')[0]
                pos_chr = gff_file.tell()-len(line)
                if name not in GFF_DESCR:
                    GFF_DESCR[name] = pos_chr
                
            line = gff_file.readline()


def parse(param=None, fasta_hash=None, chr_id=None):
    logger.title('# Parsing GFF file')

    gff_fname = param.gff_fname
    fasta_hash = fasta_hash
    chr_id = chr_id

    if not GFF_DESCR:
        set_gff_descr(gff_fname)

    logger.info('Checking chromosome IDs consistency between GFF and fasta file...')
    logger.info('')
    chrs_common = inspect.check_chrids(chrs_gff=sorted(GFF_DESCR), chrs_fasta=sorted(fasta_hash))
    chr_ids = sorted(chrs_common) if not chr_id else [chr_id]
    if chr_id and chr_id not in chrs_common:
        logger.error('Error: wrong chromosome ID')
        logger.error('')
        sys.exit(1)

    gff_data = {}
    with open(gff_fname, 'r') as gff_file:
        eof = gff_file.seek(0, 2)
        for chr_id in chr_ids:
            gff_file.seek(GFF_DESCR[chr_id], 0)
            line = gff_file.readline()
            chr_name = line.split('\t')[0]
            while chr_name == chr_id:
                if chr_name not in gff_data:
                    logger.debug('  - Reading chromosome: ' + chr_name)
                    gff_data[chr_name] = Chromosome(_id=chr_name, fasta_chr=fasta_hash[chr_id])
                    chromosome = gff_data[chr_name]
                    chromosome.source = line.split('\t')[1]

                element_type = line.split('\t')[2]
                if element_type not in ['chromosome', 'region']:

                    if not param.types_except and not param.types_only:
                        chromosome.add(gff_element=GffElement(gff_line=line, fasta_chr=fasta_hash[chr_id]))
                    else:
                        if param.types_except:
                            if element_type not in param.types_except:
                                chromosome.add(gff_element=GffElement(gff_line=line, fasta_chr=fasta_hash[chr_id]))
                        elif param.types_only:
                            if element_type in param.types_only:
                                chromosome.add(gff_element=GffElement(gff_line=line, fasta_chr=fasta_hash[chr_id]))

                line = gff_file.readline()
                if gff_file.tell() == eof:
                    break
                else:
                    chr_name = line.split('\t')[0]

    return gff_data
