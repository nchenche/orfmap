# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:59:14 2020

@author: nicolas
"""


Genecode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
                'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
                }


class Fasta_hash():
    """
    Class build to parse a genomic (i.e. nucleotide) fasta file with no need 
    to store the whole file in memory.
    
    Arguments:    
        - descriptor = {'fasta_fname': str,
                      'header_id': int,
                      'header_start_pos': int,
                      'header_len': int,
                      'seqline_len': int,
                      'seqfile_end_pos': int}
    
    Constructor:
        - filename  (str):  the fasta filename with its path 
                            (-> "/path/to/fasta_filename.fa")
        - _id       (str):  ID of the nucleotide sequence (e.g. "NC_037310.1")
        - id_pos    (int):  cursor position pointing at the first character of 
                            _id
        - id_len    (int):  number of character (including "\n") in the line 
                            containing _id
        - seq_len   (int):  number of character (including "\n") in a nucleotide
                            sequence
        - seqfile_end_pos   (int):  number of character (including "\n") in the 
                                    whole nucleotide sequence
    """
    
    nuc_comp = {'A': 'T', 'T': 'A',
             'G': 'C', 'C': 'G'}
             
    def __init__(self, descriptor):
        self.filename = descriptor['fasta_fname']
        self._id = descriptor['header_id']
        self.id_pos = descriptor['header_start_pos']
        self.id_len = descriptor['header_len']
        self.seq_len = descriptor['seqline_len']
        self.seqfile_end_pos = descriptor['seqfile_end_pos']
        self.init_nucid_max()
#        self.sequence = None
        
    def init_nucid_max(self):
        start_pos = self.id_pos
        end_pos = self.seqfile_end_pos - self.id_len - start_pos
        n_line = int(end_pos / float(self.seq_len))
        n_last_nuc = end_pos % self.seq_len
        self.nucid_max = n_line * (self.seq_len-1) + n_last_nuc
      
    def rev_comp(self, sequence):
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
            reverse_complement += self.nuc_comp[nuc.upper()]
            
        return reverse_complement        
        
    def get_seq(self, start=1, end=10, phase=0, strand='+'):
        """
        Returns the nucleotide sequence corresponding to its given coordinates
        
        Arguments:
            - start (int):  coordinate of the first nucleotide in the sequence
            - end (int):    coordinate of the last nucleotide in the sequence
            - strand ('+' or '-'):  strand of the nucleotide sequence
            - phase (0, 1 or2): number of nucleotide to "remove" in the first codon (http://gmod.org/wiki/GFF#Nesting_Features)
                            
        Returns:
            - seq (str)
        """
        end = self.nucid_max if end > self.nucid_max else end

        _from_tmp = start / (self.seq_len-1)
        if _from_tmp % 1 == 0:
            _from = int(_from_tmp) - 1
        else:
            _from = int(_from_tmp)

        to_tmp = end / (self.seq_len-1)
        if to_tmp % 1 == 0:
            to = int(to_tmp) - 1
        else:
            to = int(to_tmp)

        len_seq = end-start

        sequence = self._get_lines(_from=_from, to=to)

        start_off = self._get_offset(start)
        end_off = start_off+len_seq

        if strand == '+':
            seq = ''.join(sequence)[(start_off-1+phase):end_off]
        else:
            seq = self.rev_comp(''.join(sequence)[(start_off-1):end_off-phase])

        return seq

    def _get_line(self, n=0):
        """
        Returns the fasta file line coresponding to a given nucleotide 
        coordinate.
        
        Arguments:
            - n (int): index of the required line
            
        Returns:
            - (str)
        """
        
        if n >= 0:
            seekpos_line = self.id_pos + self.id_len + n*self.seq_len
            
            with open(self.filename, 'rb') as _file:
                _file.seek(seekpos_line, 0)
                line = _file.readline()
                
            return line.decode().strip()
        else:
            print('Error for \'n\' value: must be an integer above zero')
            
    def _get_lines(self, _from=1, to=1):
        """
        Returns all nucleotide lines between the given coordinate indexes.
        
        Arguments:
            - _from (int): index of the first required line
            - to (int): index of the last required line
            
        Returns:
            - line (list): list of str
        """
        
        line_nb = to - _from + 1
        
        lines = []
        for i in range(line_nb):
            line = self._get_line(_from + i)
            lines.append(line)
            
        return lines

    def _get_offset(self, start):
        """
        Returns the position of the nucleotide to get acccording to a given starting
        coordinate. It also corresponds to the position the cursor file has to 
        be to get in the line to get the required nucleotide.
        
        Explanation:
        Let's consider the two first lines of an imaginary fasta file:
        line 0: 'acguacguac\n' -> 10 nt + '\n' = 11 characters
        line 1: 'gucaggacgu\n' -> 10 nt + '\n' = 11 characters
        
        If you want the 13th nucleotide (c), you have to know:
            1) at which line this nucleotide is
            2) at which position this nucleotide is in the line
            
        1) idx_line_pos = int(13/10) = 1
        2) offset = 13%11 = 2 
        
        So the 13th nucleotide is at the line 1, and in this line, at the index
        (starting at 0) 2. 
        
        Arguments:
            - start (int): coordinate of the required nucleotide
            
        Returns:
            - (int): the index of the nucleotide in its line        
        """
        
        idx_line_pos = start/(self.seq_len-1)
        if idx_line_pos % 1 == 0: 
            offset = (start % self.seq_len)
            offset += int(idx_line_pos) - 1
        else:
            offset = (start % self.seq_len)
            offset += int(idx_line_pos)
            
        return offset % self.seq_len

    def translate(self, start=1, end=10, strand='+', phase=0):
        """
        Translates a nucleotide sequence from its coordinates
        
        Arguments:
            - start (int):          coordinate of the first nucleotide in the sequence
            - end (int):            coordinate of the last nucleotide in the sequence
            - strand ('+' or '-'):  strand of the nucleotide sequence
            - phase (0, 1 or 2):    number of nucleotide to "remove" in the starting codon (http://gmod.org/wiki/GFF#Nesting_Features)
                            
        Returns:
            - protein_sequence (str)
        """
        
        if isinstance(phase, int):
            offset = (3 - ((end-start+1-phase)%3) ) % 3
            
            if strand == '+':
                end = end+offset if end+offset <= self.nucid_max else end
            elif strand == '-':                
                start = start-offset if start-offset > 0 else start
                
        else:
            phase = 0    
            
        sequence = self.get_seq(start=start, end=end, phase=phase, strand=strand)
        codons = [ sequence[i:i+3] for i in range(0, len(sequence), 3) ]
        protein_sequence = ''.join([ Genecode[x.upper()] for x in codons if len(x) == 3 ])

        return protein_sequence
        
class Fasta:
    """

    Class used to read/get informations of a chromosome from a genomic fasta file

    """

    base_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    def __init__(self):
        """

        All key index file positions to read a chromosome from a genomic fasta file.

        fasta_fname (str): fasta filename
        chr: chromosome id (str)
        curpos_start: index file position of the first nucleotide (int)
        curpos_end: index file position of the last nucleotide (int)
        seq_len: length of a nucleotide line in the fasta file (int)
        line_len: length of a nucleotide line with extra-characters (int)
        off_char: difference between seq_len and line_len (int)

        """
        self.fasta_fname = None
        self.chr = None
        self.curpos_start = None
        self.curpos_end = None
        self.seq_len = None
        self.line_len = None
        self.off_char = None
        self.nucid_max = None

    def _init_nucid_max(self):
        """

        Returns the last nucleotide ID position of the chromosome.

        Returns:
            object: int

        """
        len_pos = self.curpos_end - self.curpos_start + 1
        n_line = int(len_pos / self.line_len)

        self.nucid_max = n_line * self.seq_len + len_pos % self.line_len

    def sequence(self, start=1, end=5, phase=0, strand='+'):
        """

        Returns the sequence corresponding to the given coordinates.
        The start position must be at least 1. The end position must not exceed nucid_max.

        Args:
            start (int): start position of the desired sequence
            end (int): end position of the desired sequence (int)
            phase (int): number of nucleotide to "remove" in the first codon (0, 1 or 2)
            (http://gmod.org/wiki/GFF#Nesting_Features)
            strand (str): strand of the nucleotide sequence ('+' or '-')

        Returns:
            object: str

        """
        start = start if start >= 1 else 1
        end = end if end <= self.nucid_max else self.nucid_max

        line_start = self.get_line_nucindex(index=start)
        line_end = self.get_line_nucindex(index=end)
        lines = self.get_lines(_from=line_start, to=line_end)

        len_sequence = end - start + 1
        start_pos = start%self.seq_len if start%self.seq_len else self.seq_len
        end_pos = start_pos - 1 + len_sequence

        if strand == '+':
            return ''.join(lines)[start_pos - 1 + phase:end_pos]
        else:
            return self.reverse_complement(''.join(lines)[start_pos - 1:end_pos - phase])

    def get_line_nucindex(self, index=1):
        """

        Defines the line position where to find a given nucleotide index

        Args:
            index: nucleotide number (int)

        Returns:
            line_idx: line position where to find the nucleotide (int)

        """
        line_idx_raw = index / (self.seq_len)
        line_idx = int(line_idx_raw) if line_idx_raw % 1 == 0 else int(line_idx_raw) + 1

        return line_idx

    def get_lines(self, _from=1, to=1):
        """

        Returns all lines included between the given indexes

        Args:
            _from: index of the desired starting line (int)
            to: index of the desired last line (int)

        Returns:
            lines: a list of strings

        """
        n_lines = to - _from + 1
        lines = []
        with open(self.fasta_fname, 'rb') as fasta_file:
            fasta_file.seek(self.curpos_start + self.line_len * (_from-1))
            for n in range(n_lines):
                lines.append(fasta_file.readline().decode().strip())

        return lines

    def reverse_complement(self, sequence: str) -> str:
        """

        Returns the reverse complementary sequence of a given nucleotide sequence

        Arguments:
            - sequence: nucleotide sequence (str)

        Returns:
            object: str

        """
        sequence = list(sequence)
        sequence.reverse()

        return ''.join([self.base_complement[x.upper()] for x in sequence])

    def translate(self, start=1, end=10, strand='+', phase=0):
        """

        Translates a nucleotide sequence from its coordinates

        Arguments:
            start (int): start position of the desired sequence
            end (int): end position of the desired sequence (int)
            phase (int): number of nucleotide to "remove" in the first codon (0, 1 or 2)
            (http://gmod.org/wiki/GFF#Nesting_Features)
            strand (str): strand of the nucleotide sequence ('+' or '-')

        Returns:
            - object: str

        """

        if isinstance(phase, int):
            offset = (3 - ((end - start + 1 - phase) % 3)) % 3

            if strand == '+':
                end = end + offset if end + offset <= self.nucid_max else end
            elif strand == '-':
                start = start - offset if start - offset > 0 else start

        else:
            phase = 0

        sequence = self.sequence(start=start, end=end, phase=phase, strand=strand)
        codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]

        return ''.join([Genecode[x.upper()] for x in codons if len(x) == 3])

    def index_resume(self):
        """

        A formatted string describing all key index positions stored.

        Returns:
            object: str

        """
        tab = '{:12}' * 7
        return tab.format(self.chr, self.curpos_start, self.curpos_end,
                                       self.seq_len, self.line_len, self.off_char, self.nucid_max)


def parse_fasta(fasta_filename):
    """
    Reads a fasta file and stores key index positions to parse it without the need to store the whole file in memory.

    Args:
        fasta_filename: fasta filename (str)

    Returns:
        A list of FastaIndex instances.

    """

    chr_indexes = []
    with open(fasta_filename, 'rb') as fasta_file:
        for line in fasta_file:
            if line.startswith(b'>'):
                chr_index = Fasta()
                chr_index.fasta_fname = fasta_filename
                chr_index.chr = line.decode().strip().split('>')[-1]
                chr_index.curpos_start = fasta_file.tell()

                seqline = fasta_file.readline()
                chr_index.line_len = len(seqline)
                chr_index.seq_len = len(seqline.decode().strip())
                chr_index.off_char = chr_index.line_len - chr_index.seq_len
                fasta_file.seek(-len(seqline), 1)

                if chr_indexes:
                    chr_indexes[-1].curpos_end = chr_index.curpos_start - len(line) - chr_index.off_char - 1
                    chr_indexes[-1]._init_nucid_max()

                chr_indexes.append(chr_index)

        chr_indexes[-1].curpos_end = fasta_file.tell() - chr_indexes[-1].off_char - 1
        chr_indexes[-1]._init_nucid_max()

    return { x.chr: x for x in chr_indexes }


def fasta_descriptors(fasta_fname):
    """
    Function returning key elements required to parse/read a genomic fasta file
    without the need to store the whole file in memory.
    
    Argument:
        - fasta_fname: fasta filename (str)
        
    Returns:
        - fasta_descriptor: a list of a descriptor dictionary, one per chromosome:
                - descriptor = {'fasta_fname': str,
                              'header_id': int,
                              'header_start_pos': int,
                              'header_len': int,
                              'seqline_len': int,
                              'seqfile_end_pos': int}
    
                with:
                    - fasta_fname:  the fasta filename with its path 
                                        (-> "/path/to/fasta_filename.fa")
                    - header_id:  ID of the nucleotide sequence (e.g. "NC_037310.1")
                    - header_start_pos:  cursor position pointing at the first character of 
                                        _id
                    - header_len:  number of character (including "\n") in the line 
                                    containing _id
                    - seqline_len:  number of character (including "\n") in a nucleotide
                                    sequence
                    - seqfile_end_pos:  number of character (including "\n") in the 
                                        whole nucleotide sequence
    """
    fasta_descriptor = []

    with open(fasta_fname, 'rb') as sequences_file:
        line = sequences_file.readline()
        line = line.decode()

        while line:
            if line.startswith('>'):
                descriptor = {'fasta_fname': None,
                              'header_id': None,
                              'header_start_pos': None,
                              'header_len': None,
                              'seqline_len': None,
                              'seqfile_end_pos': None}

                descriptor['fasta_fname'] = fasta_fname
                descriptor['header_id'] = line.split()[0].split('>')[-1]
                descriptor['header_len'] = len(line)
                descriptor['header_start_pos'] = sequences_file.tell() - descriptor['header_len']
                seq = sequences_file.readline()
                descriptor['seqline_len'] = len(seq)

                if fasta_descriptor:
                    fasta_descriptor[-1]['seqfile_end_pos'] = descriptor['header_start_pos'] - 1

                fasta_descriptor.append(descriptor)

                sequences_file.seek(descriptor['header_start_pos']+descriptor['header_len'], 0)

            line = sequences_file.readline()
            line = line.decode()
            
        fasta_descriptor[-1]['seqfile_end_pos'] = sequences_file.tell()
        
    return fasta_descriptor

def parse(fasta_filename):
    """
    Function generating a descriptor of a genomic fasta file and using it to create
    a dictionary of Fasta_hash instances, one per chromosome.
    
    Argument:
        - fasta_fname: fasta filename (str)
        
    Returns:
        - fasta_hash: a dictionary of Fasta_hash instances
    """
    fasta_descriptor = fasta_descriptors(fasta_filename)

    fasta_hash = {}
    for descriptor in fasta_descriptor:
        fasta_hash[descriptor['header_id']] = Fasta_hash(descriptor)
    return fasta_hash
