#!usr/bin/python3

""" * fasta2pickle.py *
"""

# Stdlib Imports
import os
import pickle
import argparse

# Constants
EXT_NAME = ".fapick"

__author__ = "Pierre Bertin"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Pierre Bertin"
__email__ = "pierre.bertin@i2bc.paris-saclay.fr"
__status__ = "Development"


def main(args):
    """Main of the program. """
    print("\n** \033[4;36mfasta2pickle\033[0m **\n")
    chromosomes_sequence = {}
    sizes = {}

    output_file = os.path.join(
        args.out_dir,
        args.genome_name + EXT_NAME
    )


    for f in args.fasta_files:
        print("    -> Parsing \033[3;33m%s\033[0m" % (f))
        # chromosomes_sequence will be populated
        fasta_to_dictionary(f, chromosomes_sequence, sizes)

    output_data = {
        "sequences":chromosomes_sequence,
        "original_fasta":args.fasta_files,
        "genome_name":args.genome_name,
        "sizes":sizes
    }

    print(" -> Saving output file: %s" % output_file)
    pickle.dump(output_data, open(output_file, "wb"))

    print("\033[3;32mDone\033[0m\n")


def fasta_to_string(fasta_file):
    """ Return the complete fasta file in a string. """
    fasta_string = ""

    with open(fasta_file, "r") as fasta:
        for line in fasta:
            fasta_string += line
    return fasta_string

def fasta_to_dictionary(fasta_file, fasta_dic, sizes):
    """ Return the fasta file content in a dictionary. """

    fasta_string = fasta_to_string(fasta_file)
    splitted = fasta_string.split(">")
    for chunk in splitted:
        chunk_split = chunk.split("\n")
        #>ref|NC_001133| [org=Saccharomyces cerevisiae]
        # [strain=S288C] [moltype=genomic] [chromosome=I]
        # Here, need to implement a way to looking for chromosome ID.
        # in chunk_split[0]
        seq_id = chunk_split[0]
        sequence = "".join(chunk_split[1:]).strip()
        if seq_id:
            sizes[seq_id] = len(sequence)
            fasta_dic[seq_id] = sequence
            print("    *\033[3;36m%10s\033[0m len:%10d" % (
                seq_id, len(sequence)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=""" ** fasta2pickle **
        The basic use of this program is to
        load genomes sequences into pickle dictionary.
        The loading can be done from a multifasta file
        with all chromosomes sequences, or from several
        files containing fasta sequences. The name of the
        chromosome needs to be explicitly given in the header
        of each sequences. The name of the genome to rename output
        file have to be given as input.
        """
    )
    parser.add_argument(
        "fasta_files",
        metavar = "file.fasta",
        nargs = "+",
        help = "fasta files containing sequences."
    )
    parser.add_argument(
        "-n", "--genome_name",
        help = "Name to give to the output pickle file (.fapick)",
        required = True
    )
    parser.add_argument(
        "-o", "--out_dir",
        metavar = "DIR",
        default = os.getcwd(),
        help = """output directory where the results files will be written.
        The output names are generated from the input files names.
        (Default = %s)""" % (os.getcwd())
    )

    args = parser.parse_args()
    main(args)
# launch 
#python3.5 /mnt/workspace/Maxime_Stage_M2/software/strali/programs/fasta2pickle.py /mnt/workspace/Maxime_Stage_M2/data/Saccer3.fasta -n Saccer3 -o /mnt/workspace/Maxime_Stage_M2/data/