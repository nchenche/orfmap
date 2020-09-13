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


def check_chrids(gff_data=None, fasta_hash=None), param=None:
    chr_fasta_not_in_gff = set(sorted(fasta_hash)).difference(set(sorted(gff_data)))
    chr_gff_not_in_fasta = set(sorted(gff_data)).difference(set(sorted(fasta_hash)))
    chr_not_in_both = chr_fasta_not_in_gff.intersection(chr_gff_not_in_fasta)

    if chr_not_in_both:
        pass #print unconsistent chromosome ids, check your files
    if chr_fasta_not_in_gff:
        print('\nWarning: the following chromosome id(s) in {} are absent in {}:'.format(chr_id,
                                                                                         param.fasta_fname,
                                                                                         param.gff_fname)
        print('\n'.join([' - ' + x for x in sorted(chr_fasta_not_in_gff)]))
        sys.exit(1)



