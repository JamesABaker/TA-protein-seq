from __future__ import division
import subprocess
import sys
import re
import matplotlib.pyplot as plt
import numpy as np



header = ">Token_header"

base_input_file = str(sys.argv[1])
comparison_input_file = str(sys.argv[2])

list_of_files = [base_input_file, comparison_input_file]

list_of_scales = ["hessa.pl", "ww.pl",
                  "eisenberg.pl", "kd.pl"]

for scale in list_of_scales:
    for file in list_of_files:
        output_filename = "%s_%s_hydrophobicity.csv" % (scale, file)

        results = []

        with open(file) as inputfile:
            for line in inputfile:
                results.append(line.strip().split(','))

        maximum_tmd_length = 0
        for entry in results:
            if entry == results[0]:
                pass
            else:
                tmh_sequence = entry[6]
                if len(tmh_sequence) > maximum_tmd_length:
                    maximum_tmd_length = len(tmh_sequence)

        max_sequence_length = maximum_tmd_length + 40

        for entry in results:
            if entry == results[0]:
                pass
            else:
                name = entry[0]
                id = entry[1]
                tmh_start_location = entry[2]
                tmh_end_location = entry[3]
                sequence = entry[4]
                tmh_sequence = entry[5]
                n_flank_sequence = entry[6]
                c_flank_sequence = entry[7]

                correction_number = 0

                #if "Outside" in str(n_terminal_start):
                #    tmh_unaltered_sequence = str(
                #        n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                #    tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                #    correction_number = (
                #        (maximum_tmd_length - len(tmh_sequence)) / 2 + (20 - len(c_flank_sequence)))
                #    tmh_segment = tmh_reversed_sequence

                #if "Inside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(
                    n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                correction_number = (
                    (maximum_tmd_length - len(tmh_sequence)) / 2 + (20 - len(n_flank_sequence)))
                tmh_segment = tmh_unaltered_sequence

                sequence = tmh_segment

                with open('KD_calc_in.txt', 'wb') as temp_fasta:
                    temp_fasta.write(header)
                    temp_fasta.write("\n")
                    temp_fasta.write(sequence)

                var = "/"
                pipe = subprocess.Popen(["perl", scale, var])
                pipe.wait()

                with open('KDcalc_out.txt', 'rb') as totalkd:
                    lines = totalkd.readlines()
                    kd_line = lines[4]
                    totalkd.close()

                    kd = str(kd_line)

                    kd = re.sub("[KD=]", '', kd)
                    result = kd
                    result = re.sub("  ", ', ', result)
                    result = re.sub(", , , ,", '', result)
                    result = re.sub(" , ", '', result)
                    result = re.sub("\n", '', result)

                result = result.replace(" -", ", -")
                result = result.replace(",, -", ", -")
                result = result.replace(",", "")
                hydrophobicity=result
                hydrophobicity = hydrophobicity.split()
                for i in range(int(correction_number)):
                    hydrophobicity.insert(0, np.nan)
                print hydrophobicity
                output_line = [id, correction_number, hydrophobicity]
                with open(output_filename, 'a') as my_file:
                    for i in output_line:
                        my_file.write(str(i))
                        my_file.write(",")
                if file == base_input_file:
                    color = "gray"
                else:
                    color = "orange"
                plt.plot(list(hydrophobicity), linestyle='-', marker='', linewidth=1, color=color)

plt.ylabel('Relative position.')
plt.show()
my_file.write("\n")
