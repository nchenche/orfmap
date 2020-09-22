# import unittest
import sys
import random

path2data = '/home/nicolas/python_workspace/git_workspace/ORFMap/ORFMap/orfmap/data/'
fasta_filename = path2data + 'Saccer.fna'

class Index:

    fasta_filename = path2data + 'Saccer.fna'

    def __init__(self):
        self.chr = None
        self.curpos_start = None
        self.curpos_end = None
        self.seq_len = None
        self.line_len = None
        self.off_char = None

    def nucid_max(self):
        len_pos = self.curpos_end - self.curpos_start + 1
        n_line = int(len_pos / self.line_len)

        return n_line * self.seq_len + len_pos % self.line_len

    def sequence(self, start=1, end=5):
        line_start = self.get_line_nucindex(index=start)
        line_end = self.get_line_nucindex(index=end)
        lines = self.get_lines(_from=line_start, to=line_end)

        len_sequence = end - start + 1
        start_pos = start%self.seq_len if start%self.seq_len else self.seq_len
        end_pos = start_pos - 1 + len_sequence

        return ''.join(lines)[start_pos - 1:end_pos]

    def get_lines(self, _from=1, to=1):
        n_lines = to - _from + 1
        lines = []
        with open(fasta_filename, 'rb') as fasta_file:
            fasta_file.seek(self.curpos_start + self.line_len * (_from-1))
            for n in range(n_lines):
                lines.append(fasta_file.readline().decode().strip())

        return lines

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

    def resume(self):
        tab = '{:12}' * 7
        return tab.format(self.chr, self.curpos_start, self.curpos_end,
                                       self.seq_len, self.line_len, self.off_char, self.nucid_max())


def fasta_parse(fasta_filename):
    chr_indexes = []
    with open(fasta_filename, 'rb') as fasta_file:
        for line in fasta_file:
            if line.startswith(b'>'):
                chr_index = Index()
                chr_index.chr = line.decode().strip().split('>')[-1]
                chr_index.curpos_start = fasta_file.tell()

                seqline = fasta_file.readline()
                chr_index.line_len = len(seqline)
                chr_index.seq_len = len(seqline.decode().strip())
                chr_index.off_char = chr_index.line_len - chr_index.seq_len
                fasta_file.seek(-len(seqline), 1)

                if chr_indexes:
                    chr_indexes[-1].curpos_end = chr_index.curpos_start - len(line) - chr_index.off_char - 1

                chr_indexes.append(chr_index)

        chr_indexes[-1].curpos_end = fasta_file.tell() - chr_indexes[-1].off_char - 1

    return chr_indexes


class Fasta:
    def __init__(self, fasta_filename):
        self.fasta_dict = parse_fasta(fasta_filename)

    def get_chromosomes(self):
        return sorted(self.fasta_dict)

    def get_line(self, chr=None, n=1):
        return self.fasta_dict[chr][n-1]

    def get_lines(self, chr=None, _from=1, to=10):
        return [ self.get_line(chr, n=x) for x in range(_from, to+1) ]

    def sequence(self, chr=None, start=None, end=None):
        if start and end:
            return ''.join(self.fasta_dict[chr])[start-1:end]
        elif start and not end:
            return ''.join(self.fasta_dict[chr])[start-1:]
        elif end and not start:
            return ''.join(self.fasta_dict[chr])[1-1:end]
        else:
            return ''.join(self.fasta_dict[chr])

def parse_fasta(fasta_filename):
    fasta_dict = {}
    with open(fasta_filename, 'rb') as fasta_file:
        for line in fasta_file:
            line = line.decode()

            if line.startswith('>'):
                chr_id = line.split('>')[-1].strip()
                fasta_dict[chr_id] = []

            else:
                fasta_dict[chr_id].append(line.strip())

    return fasta_dict


def test_getlines(fasta, chr_index):
    for i in range(500):
        rand_from = random.randint(1, 500)
        len_sequence = random.randint(2, 50)

        lines_safe = fasta.get_lines(chr='chrI', _from=rand_from, to=rand_from+len_sequence)
        lines_exp = chr_index.get_lines(_from=rand_from, to=rand_from+len_sequence)


        if lines_safe == lines_exp:
            print(rand_from, rand_from+len_sequence)
            print(lines_safe)
            print(lines_exp)
            print('')


def test_sequence(fasta, chr_index):
    for i in range(500):
        rand_start = random.randint(1, 10000)
        len_sequence = random.randint(10, 50)

        seq_safe = fasta.sequence('chrI', start=rand_start, end=rand_start+len_sequence)
        seq_exp = chr_index.sequence(start=rand_start, end=rand_start+len_sequence)


        if seq_safe != seq_exp:
            print(rand_start, rand_start+len_sequence)
            print(seq_safe)
            print(seq_exp)
            print('')


fasta_indexes = fasta_parse(fasta_filename)
for index in fasta_indexes:
    print(index.resume())
# fasta_safe = Fasta(fasta_filename)
# test_getlines(fasta, chr_indexes[0])
# test_sequence(fasta, chr_indexes[0])

# print(fasta.sequence('chrI', start=6046, end=6059))
# print(chr_indexes[0].sequence(start=6046, end=6059))
# print(chr_indexes[0].get_lines(101, 101))
#
# line_idx_raw = 6046 / 60
# line_idx = int(line_idx_raw) if line_idx_raw % 1 == 0 else int(line_idx_raw) + 1
#
# start_pos = 6046%60 if 6046%60 else 60

# class fastaParser(unittest.TestCase):
#     def setUp(self):
#         with open(fasta_filename, 'rb') as fasta_file:
#             self.line = fasta_file.readline()
#
#     def test_read(self):
#         print(self.line)
#         self.assertIn('chrI', self.line)
#
# class MyTestCase(unittest.TestCase):
#     def test_something(self):
#         self.assertEqual(True, False)
#
#
# if __name__ == '__main__':
#     unittest.main()

