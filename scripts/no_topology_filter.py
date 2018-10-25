from __future__ import division
import sys
from Bio import SeqIO
import numpy as np
import subprocess
import re
import math


def hydrophobicity_calculation(sequence):
    '''
    Calculates the hydrophobicity of a string of amino acids"
    '''
    sequence = list(sequence)
    hydrophobicitiy_conversion = {
        'A': 1.8,
        'C': 2.5,
        'D': - 3.5,
        'E': - 3.5,
        'F': 2.8,
        'G': - 0.4,
        'H': - 3.2,
        'I': 4.5,
        'K': - 3.9,
        'L': 3.8,
        'M': 1.9,
        'N': - 3.5,
        'P': - 1.6,
        'Q': - 3.5,
        'R': - 4.5,
        'S': - 0.8,
        'T': - 0.7,
        'V': 4.2,
        'W': - 0.9,
        'Y': - 1.3,
        'X': np.nan

    }
    residue_hydrophobicities = []
    for residue in sequence:
        residue_hydrophobicities.append(
            hydrophobicitiy_conversion[str(residue)])
    return np.mean(residue_hydrophobicities)


def disorder_calculation(sequence):
    '''
    Calculates the disorder of a string of amino acids. Rune and Linding disorder propensity : GlobPlot : NAR (2003) 31:3701
    '''
    sequence = list(sequence)
    disorder_conversion = {
        'A': - 0.26154,
        'C': - 0.01515,
        'D':  0.22763,
        'E': - 0.20469,
        'F': - 0.22557,
        'G':  0.43323,
        'H': - 0.00122,
        'I': - 0.42224,
        'K': - 0.10009,
        'L': - 0.33793,
        'M': - 0.22590,
        'N':  0.22989,
        'P': 0.55232,
        'Q': - 0.18768,
        'R': - 0.17659,
        'S': 0.14288,
        'T': 0.00888,
        'V': - 0.38618,
        'W': - 0.24338,
        'Y': - 0.20751,
        'X': np.nan

    }
    residue_disorder = []
    for residue in sequence:
        residue_disorder.append(disorder_conversion[str(residue)])
    return np.mean(residue_disorder)


def entropy(string):
    '''
    Use the following code for a custom command.
    via "Shannon's entropy equation is the standard method of calculation.
    Here is a simple implementation in Python, shamelessly copied from the Revelation codebase, and thus GPL licensed:"
    '''
    # get probability of chars in string
    prob = [float(string.count(c)) / len(string)
            for c in dict.fromkeys(list(string))]
    # calculate the entropy
    entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])
    return entropy

input_file = str(sys.argv[1])
flank_length = 5
# This works with uniprot filetype. From the seqIO biopython wiki:
# Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also
# the UniProt XML format if something goes wrong.
input_format = "swiss"
# For future modification, this can be used to look for any annotation in
# the file.
feature_type = "TRANSMEM"
signal_feature = "SIGNAL"
# Simply the output name, can be anything as it is written in binary (not
# file-type specific language).
output_filename_fasta = str(str(input_file) + "-TMD-filtered.fasta")

# For each file, a table is generated for each of the flank lengths set.

# File output name
output_filename = input_file.replace(".txt", ".csv")

# The header row in the file.
with open(output_filename, 'w') as my_file:
    my_file.write("Name and description, ID, tmh start location, tmh end location, tmh length, full protein sequence, tmh sequence, N flank sequence, C flank sequence, Hydrophobicity of TMH, Hydrophobicity of TMH and flanks, Disorder of TMH, Disorder of TMH and flanks, Entropy of TMH, Entropy of TMH and flanks \n")
my_file.closed

# We need to check against nearby features to prevent overlapping
# flanking regions. Note here we want to avoid clashing with INTRAMEM
# regions, since their flanking regions may be similar, however INTRAMEM regions
# should not be included in the logged transmembrane features since
# they may have extremely short transmembrane sequences which would
# disrupt very much so the alignment of the flanking regions.
unknown = 0
reliable_flank_length = flank_length

tail_errors = 0
signal_errors = 0
multipass = 0

# We iterate through each record, parsed by biopython.
for record in SeqIO.parse(input_file, input_format):
    for i, f in enumerate(record.features):
        if f.type == feature_type:
            # Calls the parsing functions as strings for later use.
            id_of_record = record.id
            name_of_record = str(
                record.description).replace(",", "")

            # Some transmembrane annotations have unknown sequence
            # positions where it is ambiguous. These features are
            # discounted.
            if "UnknownPosition" in str(f.location):
                pass
                unknown = unknown + 1
                print(id_of_record, "had an unknown position TMH.")
            else:

                full_sequence = str(record.seq)
                tmh_start = int(f.location.start)
                tmh_stop = int(f.location.end)
                tmh_sequence = str(
                    record.seq[(f.location.start):(f.location.end)])
                # These are not c or n terminal for sure. This is just an
                # assumption we make since we are not formally filtering
                # anything in this list.
                n_terminal_flank = record.seq[(
                    f.location.start + 1 - 5):(f.location.start)]
                c_terminal_flank = record.seq[(
                    f.location.end):(f.location.end + 5)]

                # This is the information that will be written for the record.
                # +/-1s are used since slices originally call how many steps to iterate rather than the sequence postion. This matches the Uniprot sequence numbering
                tmh_record = [name_of_record, id_of_record, tmh_start + 1, tmh_stop, abs(tmh_start - tmh_stop) - 1, full_sequence, tmh_sequence, n_terminal_flank, c_terminal_flank, hydrophobicity_calculation(tmh_sequence), hydrophobicity_calculation(str(
                    c_terminal_flank + tmh_sequence + n_terminal_flank)), disorder_calculation(tmh_sequence), disorder_calculation(str(c_terminal_flank + tmh_sequence + n_terminal_flank)), entropy(tmh_sequence), entropy(str(c_terminal_flank + tmh_sequence + n_terminal_flank))]

                # Some filters need to be applied to check for signal anchors
                # and erroneously long tail anchors.
                tail_length = abs(tmh_stop - len(full_sequence))

                if tail_length <= 25:
                    signal = False
                    singlepass = False
                    tmd_count = 0

                    for i, f in enumerate(record.features):
                        if f.type == feature_type:
                            tmd_count = tmd_count + 1
                    if tmd_count == 1:

                        for i, f in enumerate(record.features):
                            if f.type == signal_feature:
                                signal = True
                        if signal == False:
                            with open(output_filename, 'a') as my_file:
                                for i in tmh_record:
                                    my_file.write(str(i))
                                    my_file.write(",")
                                my_file.write("\n")
                            with open(output_filename_fasta, 'a') as filtered_fasta_file:
                                # This prevents several Fasta entries for the
                                # same record if splice isoforms exist.
                                fasta_written = True
                                fasta_record = str(
                                    ">" + str(record.id) + "\n" + str(record.seq) + "\n")
                                filtered_fasta_file.write(fasta_record)
                        elif signal == True:
                            signal_errors = signal_errors + 1
                            print("Signal error in", record.id, ". Record ",
                                  record.id, "contained a signal sequence.")
                    elif tmd_count > 1:
                        multipass = multipass + 1
                        print(tmd_count, " TMHs in Record ", record.id)

                elif tail_length > 25:
                    tail_errors = tail_errors + 1
                    print("Tail error in", record.id, ". Tail length should be 25 or lower. ",
                          record.id, "had a tail-length of ", tail_length, ".")

print(tail_errors, " proteins exceeded the tail length restriction and were excluded from the .csv and .fasta output..")
print(signal_errors, " proteins contained a signal sequence and were excluded from the .csv and .fasta output.")
print(multipass, " proteins contained a multiple TMHs and were excluded from the .csv and .fasta output.")
