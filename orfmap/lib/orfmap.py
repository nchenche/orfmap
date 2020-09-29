from orfmap.lib import logHandler
from orfmap.lib import gff_parser

logger = logHandler.Logger(name=__name__)


def mapping(gff_data, fasta_hash, param):

    all_orfs = []
    for chr_id in sorted(gff_data):
        logger.info('Reading chromosome {} ...'.format(chr_id))
        gff_chr = gff_data[chr_id]

        logger.info(' - ORF mapping and assignment')
        orfs = get_orfs(gff_chr=gff_chr, param=param)
        logger.info('')

        for orf in sorted(orfs, key=lambda x: (x.seqid, x.start)):
            if is_orf_asked(orf=orf, param=param):
                all_orfs.append(orf)

    return all_orfs


def get_orfs(gff_chr, param):
    orf_len = param.orf_len
    orfs = []
    sequence = gff_chr.sequence()
    pos = 0

    # loops on each possible frame (the negative frame is defined in "frame_rev")
    for frame in range(3):
        # list of codons in frame "frame"
        codons = [sequence[i:i + 3].upper() for i in range(frame, len(sequence), 3) if len(sequence[i:i + 3]) == 3]

        start_pos = frame + 1

        frame_rev = (gff_chr.end % 3 - frame) % 3
        start_pos_rev = None
        end_pos_rev = None

        for pos, codon in enumerate(codons):
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
                    orfs.append(orf)
                    assignment(orfs=orfs, gff_chr=gff_chr, param=param)


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
                        orfs.append(orf)
                        assignment(orfs=orfs, gff_chr=gff_chr, param=param)

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
                orfs.append(orf)
                assignment(orfs=orfs, gff_chr=gff_chr, param=param)

    return orfs


def assignment(orfs, gff_chr, param):
    orf = orfs[-1]
    elements = gff_chr.get_elements(coors=orf.get_coors(), types=param.types)
    check_ovp(orf=orf, elements=elements, co_ovp=param.co_ovp)
    set_attributes(orfs=orfs, orf_len=param.orf_len)


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


def set_attributes(orfs, orf_len):
    orf = orfs[-1]
    orf._set_type()
    orf._id = orf.format_id()
    orf._set_parent()
    orf._set_color()
    orf._set_status()

    if orf.ovp_phased:
        get_suborfs(orfs=orfs, orf_len=orf_len)


def get_suborfs(orfs, orf_len):
    """

    Checks if the sequence extremities of an ORF overlapping with a CDS are long enough to be considered as nc_ORF.
    Adds them to the orfs list if the criteria is met.

    Args:
        orfs (list): list of GffElement instances
        orf_len (int): cutoff length of the sequence to be an ORF

    Returns:
        None

    """
    orf = orfs[-1]
    gff_line = orf.get_gffline()
    suborf = None
    for element in orf.ovp_phased:
        is_ter = False
        if is_5ter_ok(orf, element, orf_len=orf_len):
            is_ter = True
            suborf = gff_parser.GffElement(gff_line=gff_line, fasta_chr=orf.fasta_chr)
            suborf.type = 'nc_5-CDS'

            if orf.strand == '+':
                suborf.end = element.get_coors()[0] - 1
            else:
                suborf.start = element.get_coors()[1] + 1

        if is_3ter_ok(orf, element, orf_len=orf_len):
            is_ter = True
            suborf = gff_parser.GffElement(gff_line=gff_line, fasta_chr=orf.fasta_chr)
            suborf.type = 'nc_3-CDS'

            if orf.strand == '+':
                suborf.start = element.get_coors()[1] + 1
            else:
                suborf.end = element.get_coors()[0] - 1

        if is_ter:
            suborf.frame = orf.frame
            suborf.status = 'non-coding'
            suborf.color = '#ffc100'
            suborf._id = suborf.format_id()
            suborf.parent = orf.format_id()
            orfs.append(suborf)


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
