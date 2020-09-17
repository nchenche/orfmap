from orfmap.lib import logHandler
from orfmap.lib import gff_parser

logger = logHandler.Logger(name=__name__)


def mapping(gff_data, fasta_hash, param):

    all_orfs = []
    for chr_id in sorted(fasta_hash):
        logger.info('Reading chromosome {} ...'.format(chr_id))
        gff_chr = gff_data[chr_id]

        logger.info(' - ORF mapping')
        orfs = gff_parser.get_orfs(gff_chr=gff_chr, orf_len=param.orf_len)

        logger.info(' - Assigning status for all {} ORFs found'.format(len(orfs)))
        logger.info('')
        for orf in sorted(orfs, key=lambda x: (x.seqid, x.start)):
            # elements are defined genomic sequences used as a reference to assign (non-)overlapping ORFs (CDS by default)
            elements = gff_chr.get_elements(coors=orf.get_coors(), types=param.types)
            orf.assignment(elements, co_ovp=param.co_ovp)
            all_orfs.append(orf)

    return all_orfs
