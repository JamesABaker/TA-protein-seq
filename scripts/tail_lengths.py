import sys
import re
import numpy as np


def maximum(array):
    maximum = 0
    for i in array:
        if i > maximum:
            maximum = i
    return(maximum)


list_of_files = sys.argv[1:]

for file_number, file in enumerate(list_of_files):
      # This generates an empty list for each potential possition. Values of
      # hydrophobicity will be added to this later and will contribute to the
      # average.
    tail_length_list = []
    results = []
    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))
    for entry in results:
        if entry == results[0]:
            pass
        else:
            name = entry[0]
            id = entry[1]
            tmh_start_location = int(entry[2])
            tmh_end_location = int(entry[3])
            tmh_length = int(entry[4])
            sequence = entry[5]
            tmh_sequence = entry[6]
            n_flank_sequence = entry[7]
            c_flank_sequence = entry[8]
            tmh_hydrophobicity = float(entry[9])
            tmh_flank_hydrophobicity = float(entry[10])
            tmh_disorder = float(entry[11])
            tmh_flank_disorder = float(entry[12])
            tmh_entropy = float(entry[13])
            tmh_flank_entropy = float(entry[14])

            tail_length = abs(tmh_end_location - len(sequence))

            tail_length_list.append(tail_length)

            if tail_length > 50:
                print(id, tail_length)

    print("Mean:", np.mean(tail_length_list), ", S.D:", np.std(
        tail_length_list), ", Max:", maximum(tail_length_list))
