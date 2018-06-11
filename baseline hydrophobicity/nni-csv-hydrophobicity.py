from __future__ import division
import sys
from Bio import SeqIO
import numpy as np
import subprocess
import re
import math

# List format:
# Name and description, ID, N terminal inside/outside, tmh start location, tmh end location, full protein sequence, tmh sequence, N flank sequence, C flank sequence, transmembrane helix sequential number, number of transmembrane helices in protein
list_of_files = sys.argv[1:]

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


print("Factor,Test, Test-statistic, P value, Bahadur slope")

list_of_files = sys.argv[1:]
disorder_sets = []
flanks_disorder_sets = []

hydrophobicity_sets = []
flanks_hydrophobicity_sets = []

entropy_sets = []
flanks_entropy_sets = []

for file_number, file in enumerate(list_of_files):
      # This generates an empty list for each potential possition. Values of
      # hydrophobicity will be added to this later and will contribute to the
      # average.
    results = []

    hydrophobicity_set = []
    flanks_hydrophobicity_set = []


    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))
    for entry in results:
        if entry == results[0]:
            pass
        else:
            name = str(entry[0])
            id = str(entry[1])
            n_terminal_start = str(entry[2])
            tmh_start_location = int(entry[3])
            tmh_end_location = int(entry[4])
            sequence = str(entry[5])
            tmh_sequence = str(entry[6])
            n_flank_sequence = str(entry[7])
            c_flank_sequence = str(entry[8])
            tmh_number = int(entry[9])
            total_tmd_count = int(entry[10])

            if total_tmd_count == 1:
                tmh_hydrophobicity = hydrophobicity_calculation(str(tmh_sequence))
                tmh_flank_hydrophobicity = hydrophobicity_calculation(str(n_flank_sequence)+str(tmh_sequence)+str(c_flank_sequence))
                hydrophobicity_set.append(tmh_hydrophobicity)
                flanks_hydrophobicity_set.append(tmh_flank_hydrophobicity)

    hydrophobicity_sets.append(hydrophobicity_set)
    flanks_hydrophobicity_sets.append(flanks_hydrophobicity_set)

print("Mean Values And Error")
for n, dataset in enumerate(list_of_files):
    print(list_of_files[n])
    print("Factor, TMH average, TMH and flanks average, TMH standard deviation, TMH and flanks standard deviation")
    print("Hydrophobicity,", np.mean(hydrophobicity_sets[n]), ",", np.mean(flanks_hydrophobicity_sets[
          n]), ",", np.std(hydrophobicity_sets[n]), ",", np.std(flanks_hydrophobicity_sets[n]))
