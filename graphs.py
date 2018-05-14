from __future__ import division
import subprocess
import sys
import re
import matplotlib.pyplot as plt
import numpy as np


header = ">Token_header"

#base_input_file = str(sys.argv[1])
#comparison_input_file = str(sys.argv[2])

list_of_files = sys.argv[1:]
print(list_of_files)

#You may need to add colours if you have many datasets to compare.
color_list=[
            ("#dde8b9"),
            ("#bfacc8"),
            ("#ffadc6"),
            ("#188fa7"),
            ("#7d70ba"),
            ]

# list_of_scales = ["hessa.pl", "ww.pl",
#                  "eisenberg.pl", "kd.pl"]

list_of_scales = ["kd.pl"]
halfway_value_for_alignment = 0
longest_tmh = 0

for file_number, file in enumerate(list_of_files):
    results = []
    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))

    halfway_value_for_alignment = 0
    for entry in results:
        if entry == results[0]:
            pass
        else:
            tmh_sequence = entry[6]
            n_flank_sequence = entry[7]
            c_flank_sequence = entry[8]

            if len(n_flank_sequence) + (len(tmh_sequence) / 2) > halfway_value_for_alignment:
                # the 0.5 is added in the case of an odd number.
                # Int(float) rounds the float down. This would
                # cause problems later since the full TMD lengths
                # may go below the 0th vector position
                halfway_value_for_alignment = int(
                    len(n_flank_sequence) + (len(tmh_sequence) / 2) + 0.5)
            if len(n_flank_sequence) + len(tmh_sequence) + len(c_flank_sequence) > longest_tmh:
                longest_tmh = len(n_flank_sequence) + len(tmh_sequence) + len(c_flank_sequence) +1 #For zero base counting

for scale in list_of_scales:
    average_scales = []
    list_of_all_hydrophobicity_vectors =[]
    for file_number, file in enumerate(list_of_files):


        # This generates an empty list for each potential possition. Values of hydrophobicity will be added to this later and will contribute to the average.
        all_hydrophobicity_vectors=[]
        for n in range(longest_tmh):
            all_hydrophobicity_vectors.append([])

        output_filename = "%s_%s_hydrophobicity.csv" % (scale, file)
        print "calculating hydrophobicity for", file
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
                tmh_start_location = entry[2]
                tmh_end_location = entry[3]
                tmh_length = entry[4]
                sequence = entry[5]
                tmh_sequence = entry[6]
                n_flank_sequence = entry[7]
                c_flank_sequence = entry[8]
                correction_number = 0

                # This assumes that the filters have already been applied and that only the n terminal will be the flank.
                tmh_unaltered_sequence = str(
                    n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)

                #The alignment seems off, check this since it could be a mathemtatical or biological artefact.
                correction_number = (
                    (halfway_value_for_alignment - len(tmh_sequence)) / 2 + (len(n_flank_sequence)))
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
                hydrophobicity = result
                hydrophobicity = hydrophobicity.split()
                for i in range(int(correction_number)):
                    hydrophobicity.insert(0, np.nan)
                hydrophobicity_converted = []
                for n, i in enumerate(hydrophobicity):
                    hydrophobicity_converted.append(float(i))
                hydrophobicity = np.array(hydrophobicity_converted)

                #Now we calculate the positions to plot on the X axis in relation to these vectors.
                #positions_for_line=[]
                #for n, item in enumerate(hydrophobicity):
                #    positions_for_line.append(n-halfway_value_for_alignment-1) #Base 0 counting

                #plt.plot(positions_for_line, hydrophobicity, linestyle='', marker='.',
                #         linewidth=0.5, color = color_list[file_number], alpha=(1/len(results)*len(results)/10))

                for position, hydrophobicity_value in enumerate(hydrophobicity):
                    #if str(hydrophobicity_value) != str("nan"):
                    all_hydrophobicity_vectors[position].append(hydrophobicity_value)

        averaged_hydrophobicities=[]
        list_of_all_hydrophobicity_vectors.append(np.array(all_hydrophobicity_vectors))
        for i in all_hydrophobicity_vectors:
            if len(i)>20:
                averaged_hydrophobicities.append(np.mean(i))
            else:
                averaged_hydrophobicities.append(np.mean(np.nan))
        average_scales.append(averaged_hydrophobicities)


    for line_number, line in enumerate(average_scales):
        positions_for_line=[]
        for n, item in enumerate(line):
            positions_for_line.append(n-halfway_value_for_alignment-1)#base 0 counting
        plt.plot(positions_for_line, line, linestyle='-', marker='',
            linewidth=1, color = color_list[line_number], alpha=1)

    for index, number in enumerate(positions_for_line):
        if number == 0:
            central_index=index

    plt.xlim([positions_for_line[central_index-15], positions_for_line[central_index+15]])
    plt.ylabel('Kyte & Doolittle hydrophobicity')
    plt.xlabel('Sequence position')
    plt.show()
