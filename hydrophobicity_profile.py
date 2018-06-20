import numpy as np
import scipy
from scipy import stats
import sys
import numpy as np

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
for file_number, file in enumerate(list_of_files):
      # This generates an empty list for each potential possition. Values of
      # hydrophobicity will be added to this later and will contribute to the
      # average.
    results = []

    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))

    #To align correctly, we need to know the longest sequence of n+tmh/2.
    length_compensation = 0
    max_length = 0
    for entry in results:
        if entry == results[0]:
            pass
        else:
            name = str(entry[0])
            id = str(entry[1])
            tmh_start_location = int(entry[2])
            tmh_end_location = int(entry[3])
            tmh_length = int(entry[4])
            sequence = str(entry[5])
            tmh_sequence = str(entry[6])
            n_flank_sequence = str(entry[7])
            c_flank_sequence = str(entry[8])
            tmh_hydrophobicity = float(entry[9])
            tmh_flank_hydrophobicity = float(entry[10])
            tmh_disorder = float(entry[11])
            tmh_flank_disorder = float(entry[12])
            tmh_entropy = float(entry[13])
            tmh_flank_entropy = float(entry[14])

            if len(str(n_flank_sequence+tmh_sequence/2)) > length_compensation:
                length_compensation=len(str(n_flank_sequence+tmh_sequence/2))
            if len(str(n_flank_sequence+tmh_sequence+c_flank_sequence)) > max_length:
                max_length=len(str(n_flank_sequence+tmh_sequence+c_flank_sequence))


    # Now we will build an empty list based on the max length.
    hydrophobicity_in_position = []
    for n in range(max_length):
        hydrophobicity_in_position.append([])

    for entry in results:
        if entry == results[0]:
            pass
        else:
            name = str(entry[0])
            id = str(entry[1])
            tmh_start_location = int(entry[2])
            tmh_end_location = int(entry[3])
            tmh_length = int(entry[4])
            sequence = str(entry[5])
            tmh_sequence = str(entry[6])
            n_flank_sequence = str(entry[7])
            c_flank_sequence = str(entry[8])
            tmh_hydrophobicity = float(entry[9])
            tmh_flank_hydrophobicity = float(entry[10])
            tmh_disorder = float(entry[11])
            tmh_flank_disorder = float(entry[12])
            tmh_entropy = float(entry[13])
            tmh_flank_entropy = float(entry[14])

            hydrophobic_profile=[]
            for i in str(n_flank_sequence+tmh_sequence+c_flank_sequence):
                hydrophobic_profile.append(hydrophobicity_calculation(i))
            for n in range(length_compensation-len(str(n_flank_sequence+tmh_sequence/2))):
                insert(0, np.nan)
            for n, i in enumerate(hydrophobic_profile):
                hydrophobicity_in_position[n].append(i)

    hydrophobicity_in_position_averages=[]
    for list in hydrophobicity_in_position:
        hydrophobicity_in_position_averages.append(np.mean(list))
    print(file)
    print(hydrophobicity_in_position_averages)
