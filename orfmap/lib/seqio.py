import os


def write_orfs(all_orfs: list, param=None):
    """
    Writes fasta and gff files for orfs.

    Args:
        all_orfs: list of lib.gff_parser.Gff_element instances
        param: instance of lib.parameters.Param

    Returns:
        None

    """

    with open(param.outfile + '.gff', "w") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(param.fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(param.gff_fname))
        out_gff.write(header)
        with open(param.outfile + '.fa', "w") as out_fasta:
            for orf in [ x for x in all_orfs if is_orf_asked(orf=x, param=param)]:
                out_gff.write(orf.get_gffline())
                out_fasta.write(orf.get_fastaline())


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

