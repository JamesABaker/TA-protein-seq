from __future__ import division
import sys
from Bio import SeqIO
import numpy as np
import subprocess
import re
import math

list_of_files = ["reference_list.csv"]
input_files = sys.argv[1:]
flank_length = 5
# This works with uniprot filetype. From the seqIO biopython wiki:
# Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also
# the UniProt XML format if something goes wrong.
input_format = "swiss"
# For future modification, this can be used to look for any annotation in
# the file.

results = []
for file_number, reference_file in enumerate(list_of_files):
    with open(reference_file) as reference_file_check:
        for line in reference_file_check:
            results.append(line.strip().split(','))

for input_file_number, input_file in enumerate(input_files):
    disorder_set = []
    hydrophobicity_set = []
    flanks_disorder_set = []
    flanks_hydrophobicity_set = []
    entropy_set = []
    flanks_entropy_set = []
    # We iterate through each record, parsed by biopython.
    for record in SeqIO.parse(input_file, input_format):
        for entry in results:
            if entry == results[0]:
                pass
            else:
                name = str(entry[0])
                id = str(entry[1])
                tmh_start_location = entry[2]
                tmh_end_location = entry[3]
                tmh_length = entry[4]
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
                source = str(entry[15])

                if id == str(record.id):
                    hydrophobicity_set.append(tmh_hydrophobicity)
                    flanks_hydrophobicity_set.append(tmh_flank_hydrophobicity)
                    disorder_set.append(tmh_disorder)
                    flanks_disorder_set.append(tmh_flank_disorder)
                    entropy_set.append(tmh_entropy)
                    flanks_entropy_set.append(tmh_flank_entropy)

    print("\n\n", input_file)
    print(len(hydrophobicity_set), " records captured.")
    print("Factor, TMH average, TMH and flanks average, TMH standard deviation, TMH and flanks standard deviation")
    print("Hydrophobicity,", np.mean(hydrophobicity_set), ",", np.mean(flanks_hydrophobicity_set), ",", np.std(hydrophobicity_set), ",", np.std(flanks_hydrophobicity_set))
    print("Disorder,", np.mean(disorder_set), ",", np.mean(flanks_disorder_set), ",", np.std(disorder_set), ",", np.std(flanks_disorder_set))
    print("Sequence Entropy,", np.mean(entropy_set), ",", np.mean(flanks_entropy_set), ",", np.std(entropy_set), ",", np.std(flanks_entropy_set))
