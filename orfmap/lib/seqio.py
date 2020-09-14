import os


def write_orfs(all_orfs=[], param=None):
    """
    Writes fasta and gff files for orfs.

    Args:
        all_orfs: list of lib.gff_parser.Gff_element instances
        param: instance of lib.parameters.Param

    Returns:
        None

    """

    print('\n-----------------------')
    print('#Writing output files (gff and fasta)')
    with open(param.outfile + '.gff', "w") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(param.fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(param.gff_fname))
        out_gff.write(header)
        with open(param.outfile + '.fa', "w") as out_fasta:
            for orf in all_orfs:
                out_gff.write(orf.get_gffline())
                out_fasta.write(orf.get_fastaline())