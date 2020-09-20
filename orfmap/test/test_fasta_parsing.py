# import unittest
import sys

path2data = '/home/nicolas/python_workspace/git_workspace/ORFMap/ORFMap/orfmap/data/'
fasta_filename = path2data + 'Saccer.fna'

class Index:
    def __init__(self):
        self.chr = None
        self.curpos_start = None
        self.curpos_end = None
        self.seq_len = None
        self.line_len = None
        self.off_char = None


    def resume(self):
        tab = '{:12}' * 6
        return tab.format(self.chr, self.curpos_start, self.curpos_end,
                                       self.seq_len, self.line_len, self.off_char)

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

    fasta_file.seek(237896)
    print(fasta_file.readline())

for index in chr_indexes:
    print(index.resume())

sys.exit(0)

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

fasta = Fasta(fasta_filename)
print(len(fasta.sequence('chrI')))
sys.exit(0)
print(fasta.get_lines('chrI', _from=2, to=3))
print(fasta.get_line('chrI', n=2))
print('')
# print(fasta.sequence('chrI', 1, 5))
# print(fasta.sequence('chrI', 1, 10))


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

