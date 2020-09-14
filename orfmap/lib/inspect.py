import sys

def check_types(gff_data=None, types=[]):
    unconsistent_types = []
    all_types = get_types(gff_data)

    for _type in types:
        if _type not in all_types:
            unconsistent_types.append(_type)

    if unconsistent_types:
        print('Warning: wrong type(s) has(have) been given:')
        print('\n'.join([' - ' + x for x in unconsistent_types]))
        print('\nYou can choose amongst the types listed below:')
        print('\n'.join([' - ' + x for x in all_types])+'\n')
        sys.exit(1)


def get_types(gff_data):
    all_types_2d = [ gff_data[x].get_types() for x in (sorted(gff_data)) ]
    all_types_flatten = (sorted(set([val for sublist in all_types_2d for val in sublist])))

    return all_types_flatten

def check_chrids(chrs_gff=[], chrs_fasta=[]):
    chr_common = set(chrs_gff).intersection(chrs_fasta)

    if chr_common:
        if len(chr_common) != len(chrs_fasta):
            print('\nWarning: all chromosomes are not shared between gff and fasta files.\n')
            table_chrs(chrs_gff, chrs_fasta)
        else:
            print('\nAll chromosomes are shared between gff and fasta files.\n')
            table_chrs(chrs_gff, chrs_fasta)
        return chr_common
    else:
        print('\nError: chromosomes are not consistent between gff and fasta files.\n')
        table_chrs(chrs_gff, chrs_fasta)
        sys.exit(1)

def table_chrs(chrs_gff, chrs_fasta):
    table_header = ["Chromosome ids", "in GFF", "in fasta"]
    spacer_len = 20
    table_border = spacer_len * len(table_header) * '-'
    row_format = ('{:>'+str(spacer_len)+'}') * (len(table_header))
    print(row_format.format(*table_header))
    print(table_border)

    all_chrs = sorted(set(chrs_gff + chrs_fasta))
    for chr in all_chrs:
        if chr in chrs_gff and chr in chrs_fasta:
            print(row_format.format(chr, 'X', 'X'))
        elif chr in chrs_gff and chr not in chrs_fasta:
            print(row_format.format(chr, 'X', '-'))
        elif chr in chrs_fasta and chr not in chrs_gff:
            print(row_format.format(chr, '-', 'X'))
    print(table_border+'\n')




