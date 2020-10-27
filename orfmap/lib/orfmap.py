import os
from orfmap.lib import logHandler
from orfmap.lib import gff_parser


logger = logHandler.Logger(name=__name__)


def mapping(gff_data, param):

    if os.path.exists(param.outfile + '.gff'):
        os.remove(param.outfile + '.gff')
    if os.path.exists(param.outfile + '.fa'):
        os.remove(param.outfile + '.fa')

    with open(param.outfile + '.gff', "a+") as out_gff:
        header = '# Input genomic fasta file: {}\n'.format(os.path.basename(param.fasta_fname))
        header += '# Input gff file: {}\n'.format(os.path.basename(param.gff_fname))
        out_gff.write(header)
        with open(param.outfile + '.fa', "a+") as out_fasta:

            for chr_id in sorted(gff_data):
                logger.info('Reading chromosome {} ...'.format(chr_id))
                gff_chr = gff_data[chr_id]

                logger.info(' - ORF mapping and assignment')
                get_orfs(gff_chr=gff_chr, param=param, outfiles=[out_gff, out_fasta])
                logger.info('')


def get_orfs(gff_chr, param, outfiles: list):
    max_subsequence_length = 1999998
    out_gff = outfiles[0]
    out_fasta = outfiles[1]
    orf_len = param.orf_len
    pos = 0

    # loops on each possible frame (the negative frame is defined in "frame_rev")
    for frame in range(3):
        subsequences = (gff_chr.sequence(start=i, end=i+max_subsequence_length) for i in
                        range(1+frame, gff_chr.end, max_subsequence_length))

        frame_rev = (gff_chr.end % 3 - frame) % 3
        start_pos = None
        start_pos_rev = None
        end_pos_rev = None

        for n, sequence in enumerate(subsequences):

            logger.debug('    - reading frame {} of subsequence {}'.format(str(frame), str(n)))
            # list of codons in frame "frame"
            codons = (sequence[i:i + 3].upper() for i in range(0, len(sequence), 3) if len(sequence[i:i + 3]) == 3)

            start_pos = start_pos if start_pos else frame + 1

            for pos, codon in enumerate(codons, start=int(n*max_subsequence_length/3)):
                if codon in ['TAG', 'TGA', 'TAA']:
                    end_pos = pos * 3 + 1 + 2 + frame
                    if end_pos - start_pos + 1 >= orf_len:
                        orf = gff_parser.GffElement(fasta_chr=gff_chr.fasta_chr)
                        orf.seqid = gff_chr._id
                        orf.source = gff_chr.source
                        orf.strand = '+'
                        orf.frame = frame
                        if start_pos == frame + 1:
                            orf.start = start_pos
                        else:
                            orf.start = start_pos + 3
                        orf.end = end_pos

                        suborfs = assignment(orf=orf, gff_chr=gff_chr, param=param)
                        out_gff.write(orf.get_gffline())
                        out_fasta.write(orf.get_fastaline())
                        if suborfs:
                            for suborf in suborfs:
                                out_gff.write(suborf.get_gffline())
                                out_fasta.write(suborf.get_fastaline())

                    start_pos = end_pos - 2

                elif codon in ['CTA', 'TCA', 'TTA']:
                    if start_pos_rev is None:
                        start_pos_rev = pos * 3 + 1 + frame
                    else:
                        end_pos_rev = pos * 3 + 1 + 2 + frame
                        if end_pos_rev - start_pos_rev + 1 >= orf_len:
                            orf = gff_parser.GffElement(fasta_chr=gff_chr.fasta_chr)
                            orf.seqid = gff_chr._id
                            orf.source = gff_chr.source
                            orf.strand = '-'
                            orf.frame = frame_rev
                            orf.start = start_pos_rev
                            orf.end = end_pos_rev - 3

                            suborfs = assignment(orf=orf, gff_chr=gff_chr, param=param)
                            out_gff.write(orf.get_gffline())
                            out_fasta.write(orf.get_fastaline())
                            if suborfs:
                                for suborf in suborfs:
                                    out_gff.write(suborf.get_gffline())
                                    out_fasta.write(suborf.get_fastaline())

                        start_pos_rev = end_pos_rev - 2

        # adds coordinates of ORF in extremities
        if end_pos_rev:
            start_pos_rev = end_pos_rev - 2
            end_pos_rev = pos * 3 + 1 + 2 + frame
            if end_pos_rev - start_pos_rev + 1 >= orf_len:
                orf = gff_parser.GffElement(fasta_chr=gff_chr.fasta_chr)
                orf.seqid = gff_chr._id
                orf.source = gff_chr.source
                orf.strand = '-'
                orf.frame = frame_rev
                orf.start = start_pos_rev
                orf.end = end_pos_rev

                suborfs = assignment(orf=orf, gff_chr=gff_chr, param=param)
                out_gff.write(orf.get_gffline())
                out_fasta.write(orf.get_fastaline())
                if suborfs:
                    for suborf in suborfs:
                        out_gff.write(suborf.get_gffline())
                        out_fasta.write(suborf.get_fastaline())


def assignment(orf, gff_chr, param):
    elements = gff_chr.get_elements(coors=orf.get_coors())
    check_ovp(orf=orf, elements=elements, co_ovp=param.co_ovp)
    suborfs = set_attributes(orf=orf, orf_len=param.orf_len)

    return suborfs


def check_ovp(orf, elements, co_ovp=0.7):
    """

    Adds overlapping elements to the ORF to either orf.ovp_phased or orf.ovp_unphased lists.

    Args:
        orf (GffElement): a sequence from stop to stop codons
        elements (list[GffElement]): elements of the GFF input files
        co_ovp (float): cutoff used to define an overlapping element with an ORF

    Returns:
        None

    """
    orf_ovp_max = -1
    if elements:
        for element in elements:
            orf_ovp, element_ovp = gff_parser.get_overlap(orf_coors=orf.get_coors(), other_coors=element.get_coors())
            if element_ovp == 1.0 or orf_ovp >= co_ovp:
                if isinstance(element.phase, int) and element.frame == orf.frame and element.strand == orf.strand:
                    if element not in orf.ovp_phased:
                        orf.ovp_phased.append(element)
                else:
                    if element not in orf.ovp_unphased:
                        if orf_ovp > orf_ovp_max:
                            orf_ovp_max = orf_ovp
                            orf.ovp_unphased.append(element)


def set_attributes(orf, orf_len):
    suborfs = []
    orf._set_type()
    orf._id = orf.format_id()
    orf._set_parent()
    orf._set_color()
    orf._set_status()

    if orf.ovp_phased:
        suborfs = get_suborfs(orf=orf, orf_len=orf_len)

    return suborfs


def get_suborfs(orf: gff_parser.GffElement, orf_len):
    """

    Checks if the sequence extremities of an ORF overlapping with a CDS are long enough to be considered as nc_ORF.
    Adds them to the orfs list if the criteria is met.

    Args:
        orf (list): list of GffElement instances
        orf_len (int): cutoff length of the sequence to be an ORF

    Returns:
        None

    """
    gff_line = orf.get_gffline()
    suborfs = []

    for element in orf.ovp_phased:
        if is_5ter_ok(orf, element, orf_len=orf_len):
            suborf_5ter = gff_parser.GffElement(gff_line=gff_line, fasta_chr=orf.fasta_chr)
            suborf_5ter.type = 'nc_5-CDS'

            if orf.strand == '+':
                suborf_5ter.end = element.get_coors()[0] - 1
            else:
                suborf_5ter.start = element.get_coors()[1] + 1

            suborf_5ter.frame = orf.frame
            suborf_5ter.status = 'non-coding'
            suborf_5ter.color = '#ffc100'
            suborf_5ter._id = suborf_5ter.format_id()
            suborf_5ter.parent = orf.format_id()
            suborfs.append(suborf_5ter)

        if is_3ter_ok(orf, element, orf_len=orf_len):
            suborf_3ter = gff_parser.GffElement(gff_line=gff_line, fasta_chr=orf.fasta_chr)
            suborf_3ter.type = 'nc_3-CDS'

            if orf.strand == '+':
                suborf_3ter.start = element.get_coors()[1] + 1
            else:
                suborf_3ter.end = element.get_coors()[0] - 1

            suborf_3ter.frame = orf.frame
            suborf_3ter.status = 'non-coding'
            suborf_3ter.color = '#ffc100'
            suborf_3ter._id = suborf_3ter.format_id()
            suborf_3ter.parent = orf.format_id()
            suborfs.append(suborf_3ter)

    return suborfs


def is_5ter_ok(orf, element, orf_len=60):
    if orf.strand == '+':
        return element.get_coors()[0] - 1 - orf.start + 1 >= orf_len
    else:
        return orf.end - element.get_coors()[1] + 1 + 1 >= orf_len


def is_3ter_ok(orf, element, orf_len=60):
    if orf.strand == '+':
        return orf.end - element.get_coors()[1] + 1 + 1 >= orf_len
    else:
        return element.get_coors()[0] - 1 - orf.start + 1 >= orf_len


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
