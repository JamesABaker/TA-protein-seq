# Functions to import

from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy
from scipy import stats
import numpy as np
import pylab
import itertools
from time import gmtime, strftime
from random import random, seed
from matplotlib import cm
from matplotlib import rcParams
import operator
import sys


# The total list of files to be processed

list_of_files = sys.argv[1:]
'''
list_of_files = ["TopDB_10_flanklength_flankclashFalse.csv",
                 "TopDB_10_flanklength_flankclashTrue.csv",
                 "TopDB_20_flanklength_flankclashFalse.csv",
                 "TopDB_20_flanklength_flankclashTrue.csv",
                 "TopDB_5_flanklength_flankclashFalse.csv",
                 "TopDB_5_flanklength_flankclashTrue.csv",
                 "UniArch_10_flanklength_flankclashFalse.csv",
                 "UniArch_10_flanklength_flankclashTrue.csv",
                 "UniArch_20_flanklength_flankclashFalse.csv",
                 "UniArch_20_flanklength_flankclashTrue.csv",
                 "UniArch_5_flanklength_flankclashFalse.csv",
                 "UniArch_5_flanklength_flankclashTrue.csv",
                 "UniBacilli_10_flanklength_flankclashFalse.csv",
                 "UniBacilli_10_flanklength_flankclashTrue.csv",
                 "UniBacilli_20_flanklength_flankclashFalse.csv",
                 "UniBacilli_20_flanklength_flankclashTrue.csv",
                 "UniBacilli_5_flanklength_flankclashFalse.csv",
                 "UniBacilli_5_flanklength_flankclashTrue.csv",
                 "UniCress_10_flanklength_flankclashFalse.csv",
                 "UniCress_10_flanklength_flankclashTrue.csv",
                 "UniCress_20_flanklength_flankclashFalse.csv",
                 "UniCress_20_flanklength_flankclashTrue.csv",
                 "UniCress_5_flanklength_flankclashFalse.csv",
                 "UniCress_5_flanklength_flankclashTrue.csv",
                 "UniEcoli_10_flanklength_flankclashFalse.csv",
                 "UniEcoli_10_flanklength_flankclashTrue.csv",
                 "UniEcoli_20_flanklength_flankclashFalse.csv",
                 "UniEcoli_20_flanklength_flankclashTrue.csv",
                 "UniEcoli_5_flanklength_flankclashFalse.csv",
                 "UniEcoli_5_flanklength_flankclashTrue.csv",
                 "UniER_10_flanklength_flankclashFalse.csv",
                 "UniER_10_flanklength_flankclashTrue.csv",
                 "UniER_20_flanklength_flankclashFalse.csv",
                 "UniER_20_flanklength_flankclashTrue.csv",
                 "UniER_5_flanklength_flankclashFalse.csv",
                 "UniER_5_flanklength_flankclashTrue.csv",
                 "UniFungi_10_flanklength_flankclashFalse.csv",
                 "UniFungi_10_flanklength_flankclashTrue.csv",
                 "UniFungi_20_flanklength_flankclashFalse.csv",
                 "UniFungi_20_flanklength_flankclashTrue.csv",
                 "UniFungi_5_flanklength_flankclashFalse.csv",
                 "UniFungi_5_flanklength_flankclashTrue.csv",
                 "UniGolgi_10_flanklength_flankclashFalse.csv",
                 "UniGolgi_10_flanklength_flankclashTrue.csv",
                 "UniGolgi_20_flanklength_flankclashFalse.csv",
                 "UniGolgi_20_flanklength_flankclashTrue.csv",
                 "UniGolgi_5_flanklength_flankclashFalse.csv",
                 "UniGolgi_5_flanklength_flankclashTrue.csv",
                 "UniHuman_10_flanklength_flankclashFalse.csv",
                 "UniHuman_10_flanklength_flankclashTrue.csv",
                 "UniHuman_20_flanklength_flankclashFalse.csv",
                 "UniHuman_20_flanklength_flankclashTrue.csv",
                 "UniHuman_5_flanklength_flankclashFalse.csv",
                 "UniHuman_5_flanklength_flankclashTrue.csv",
                 "UniPM_10_flanklength_flankclashFalse.csv",
                 "UniPM_10_flanklength_flankclashTrue.csv",
                 "UniPM_20_flanklength_flankclashFalse.csv",
                 "UniPM_20_flanklength_flankclashTrue.csv",
                 "UniPM_5_flanklength_flankclashFalse.csv",
                 "UniPM_5_flanklength_flankclashTrue.csv",
                 ]
'''
'''
# For figure 1
list_of_files = [
    "UniHuman_5_flanklength_flankclashTrue.csv",
    "TopDB_5_flanklength_flankclashTrue.csv",
]

# For Figure 2,3,& 4
list_of_files = [
    "UniHuman_20_flanklength_flankclashTrue.csv",
    "TopDB_20_flanklength_flankclashTrue.csv",
    "UniHuman_10_flanklength_flankclashTrue.csv",
    "TopDB_10_flanklength_flankclashTrue.csv",
]

# For Figure 6
list_of_files = [
    "TopDB_20_flanklength_flankclashTrue.csv",
    "UniHuman_20_flanklength_flankclashTrue.csv",
    "UniER_20_flanklength_flankclashTrue.csv",
    "UniGolgi_20_flanklength_flankclashTrue.csv",
    "UniPM_20_flanklength_flankclashTrue.csv",
    "UniFungi_20_flanklength_flankclashTrue.csv",
    "UniCress_20_flanklength_flankclashTrue.csv",
    "UniEcoli_20_flanklength_flankclashTrue.csv",
    "UniBacilli_20_flanklength_flankclashTrue.csv",
    "UniEcoli_20_flanklength_flankclashTrue.csv",
    "UniArch_20_flanklength_flankclashTrue.csv",
]

# For NEW FIGURE
list_of_files = [
    "UniHuman_5_flanklength_flankclashTrue.csv",
    "TopDB_5_flanklength_flankclashTrue.csv",
]
'''
# The type of TM proteins.
type_of_proteins = ["multipass", "simple", "complex"]
type_of_proteins = ["singlepass", "multipass", ]


# To avoid repetition, the residues will be held in a list to be looped
# through.
residues = ["I", "V", "L", "F", "C", "M", "A", "G", "T",
            "S", "W", "Y", "P", "H", "E", "Q", "D", "N", "K", "R"]


# Function for checking if the entry has 1 or more TMHs depending on if a
# single-pass result or multi-pass result is needed.
def single_or_multi(type_of_protein, total_tmd_count, record_id):
    if type_of_protein == "singlepass" and total_tmd_count == 1:
        return True
    elif type_of_protein == "multipass" and total_tmd_count > 1:
        return True
    elif type_of_protein == "simple" and total_tmd_count == 1:
        with open('List_of_Simple_ids.txt') as f:
            simple_id_list = f.read().splitlines()
            if record_id in simple_id_list:
                return True
            else:
                return False
    elif type_of_protein == "complex" and total_tmd_count == 1:
        with open('List_of_Complex_ids.txt') as f:
            complex_id_list = f.read().splitlines()
            if record_id in complex_id_list:
                return True
            else:
                return False
    else:
        return False

# A function to find the compensation size to correctly align the TMHs at
# position 0 in the middle of the TMH.

# We need to know if the N terminal flank or the C terminal
# flank is appropriate for finding the longest half to
# align the vectors to.


def compensation_checker(tmh_sequence, n_flank_sequence, n_terminal_start, halfway_value_for_alignment):
    '''
    Returns how many residues are needed to align the TMH to the central 0th residue.
    '''
    if n_terminal_start == "Inside":
        if len(n_flank_sequence) + len(tmh_sequence) / 2 > halfway_value_for_alignment:
            # the 0.5 is added in the case of an odd number.
            # Int(float) rounds the float down. This would
            # cause problems later since the full TMD lengths
            # may go below the 0th vector position
            halfway_value_for_alignment = int(
                len(n_flank_sequence) + len(tmh_sequence) / 2 + 0.5)
    elif n_terminal_start == "Outside":
        if len(c_flank_sequence) + len(tmh_sequence) / 2 > halfway_value_for_alignment:
            halfway_value_for_alignment = int(
                len(c_flank_sequence) + len(tmh_sequence) / 2 + 0.5)
    else:
        halfway_value_for_alignment = halfway_value_for_alignment
    return halfway_value_for_alignment

# This returns the sequence according to the inside-> outside orientation
# with a string of "Js" that count for nothing when calculating residue
# abundances.


def oriented_sequence(n_terminal_start, n_flank_sequence, tmh_sequence, c_flank_sequence, halfway_value_for_alignment):
    '''
    Returns the sequence oriented from inside to outside rather than N to C.
    '''
    if "Outside" in str(n_terminal_start):
        tmh_unaltered_sequence = str(
            n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
        tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
        correction_number = halfway_value_for_alignment - \
            ((len(tmh_sequence) / 2) + len(c_flank_sequence))
        # +1 used since python is counting from 0, however I'm not sure entirely that it's needed here.
        tmh_segment = "J" * (int(correction_number) + 1)
        tmh_segment = tmh_segment + tmh_reversed_sequence
    elif "Inside" in str(n_terminal_start):
        tmh_unaltered_sequence = str(
            n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
        correction_number = halfway_value_for_alignment - \
            ((len(tmh_sequence) / 2) + len(n_flank_sequence))
        tmh_segment = "J" * (int(correction_number) + 1)
        tmh_segment = tmh_segment + tmh_unaltered_sequence
    return tmh_segment

# cleans a list for csv compatible output.


def clean_list_output(list):
    '''
    Returns a clean string from a list.
    '''
    list = ''.join(str(list))
    list = list.replace("[", "")
    list = list.replace("]", "")
    return list

# raw, absolute vectors foor each residue at each position.


def vectors_for_neg30_to_pos30(amino_acid_tally, halfway_value_for_alignment):
    '''
    Returns the vectors from -30 to +30 and the totals in those slices for each residue. This is no longer used.
    Positions are recounted in case the last time it was called was a previous dataset.
    We also need to check this belongs to a longer flanking region dataset, otherwise we don't need to print this for the bahadur calculations.
    '''
    print "\n\nRAW VECTORS from positions -30 to 30."
    if len(amino_acid_tally[1][4]) > 60:
        sequence_position = []
        # Getting the position coordinates
        for number, item in enumerate(range(0, len(amino_acid_tally[1][4]))):
            value_for_position = number - halfway_value_for_alignment
            sequence_position.append(value_for_position)
        # The values go to +31 because we need to include the 0th position
        positions = sequence_position[
            halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]
        positions = clean_list_output(positions)

        # print "\n", amino_acid_tally[1][0], amino_acid_tally[1][1]
        # print "Residue type, Total residues,", "Position,", positions

        # Printing the absolute vectors. These are useful for Bahadur values to
        # be calculated externally.
        for residue_types in amino_acid_tally:
            # This simply cleans up the output for the terminal
            absolute_values = residue_types[4][
                halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]
            absolute_values = clean_list_output(absolute_values)

            print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",", residue_types[2], ",", sum(residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]), ",,", absolute_values

        # Now we get the absolute vectors for combinations R & K, and D & E.
        for residue_types in amino_acid_tally:
            if str(residue_types[2]) == "R":
                for other_residue_types in amino_acid_tally:
                    if str(other_residue_types[2]) == "K":
                        absolute_values = []
                        for n, value in enumerate(residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]):
                            this_absolute_value = value + \
                                other_residue_types[4][
                                    halfway_value_for_alignment - 30:halfway_value_for_alignment + 31][n]
                            absolute_values.append(this_absolute_value)
                        absolute_values = clean_list_output(absolute_values)
                        print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",", "R & K,", sum(residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]) + sum(other_residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]), ",,", absolute_values
            elif str(residue_types[2]) == "D":
                for other_residue_types in amino_acid_tally:
                    if str(other_residue_types[2]) == "E":
                        absolute_values = []
                        for n, value in enumerate(residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]):
                            this_absolute_value = value + \
                                other_residue_types[4][
                                    halfway_value_for_alignment - 30:halfway_value_for_alignment + 31][n]
                            absolute_values.append(this_absolute_value)
                        absolute_values = clean_list_output(absolute_values)
                        print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",",  "D & E,", sum(residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]) + sum(other_residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]), ",,", absolute_values


def normalised_vectors_for_neg30_to_pos30(amino_acid_tally, halfway_value_for_alignment):
    '''
    Returns the normalised vectors from -30 to +30. This is no longer used since we switched to -20 to +20 and process the data externally for figure 4. The principal is the same.
    Positions are recounted in case the last time it was called was a previous dataset.
    We also need to check this belongs to a longer flanking region dataset, otherwise we don't need to print this for the bahadur calculations.
    '''

    print "\n\nRAW VECTORS from positions -30 to 30."

    if len(amino_acid_tally[1][4]) > 60:
        sequence_position = []
        # Getting the position coordinates
        for number, item in enumerate(range(0, len(amino_acid_tally[1][4]))):
            value_for_position = number - halfway_value_for_alignment
            sequence_position.append(value_for_position)
        # The values go to +31 because we need to include the 0th position
        positions = sequence_position[
            halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]
        positions = clean_list_output(positions)

        # print "\n", amino_acid_tally[1][0], amino_acid_tally[1][1]
        # print "Residue type, Total residues,", "Position,", positions

        # Printing the normalised vectors. These are useful for Bahadur values to
        # be calculated externally.
        for residue_types in amino_acid_tally:
            # This simply cleans up the output for the terminal
            normalised_values = residue_types[5][
                halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]
            normalised_values = clean_list_output(normalised_values)

            print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",", residue_types[2], ",", sum(residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]), ",,", normalised_values

        # Now we get the absolute vectors for combinations R & K, and D & E.
        for residue_types in amino_acid_tally:
            if str(residue_types[2]) == "R":
                for other_residue_types in amino_acid_tally:
                    if str(other_residue_types[2]) == "K":
                        normalised_values = []
                        for n, value in enumerate(residue_types[5][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]):
                            this_normalised_value = value + \
                                other_residue_types[5][
                                    halfway_value_for_alignment - 30:halfway_value_for_alignment + 31][n]
                            normalised_values.append(this_normalised_value)
                        normalised_values = clean_list_output(
                            normalised_values)
                        print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",", "R & K,", sum(residue_types[5][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]) + sum(other_residue_types[5][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]), ",,", normalised_values
            elif str(residue_types[2]) == "D":
                for other_residue_types in amino_acid_tally:
                    if str(other_residue_types[2]) == "E":
                        normalised_values = []
                        for n, value in enumerate(residue_types[5][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]):
                            this_normalised_value = value + \
                                other_residue_types[5][
                                    halfway_value_for_alignment - 30:halfway_value_for_alignment + 31][n]
                            normalised_values.append(this_normalised_value)
                        normalised_values = clean_list_output(
                            normalised_values)
                        print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",",  "D & E,", sum(residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]) + sum(other_residue_types[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]), ",,", normalised_values

        print "\nMax & Mins"
        for residue_types in amino_acid_tally:
            # This simply cleans up the output for the terminal
            normalised_values = residue_types[5][
                halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]
            min_index, min_value = min(
                enumerate(normalised_values), key=operator.itemgetter(1))
            max_index, max_value = max(
                enumerate(normalised_values), key=operator.itemgetter(1))
            print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",", residue_types[2], ", max ,", max(normalised_values), ",position,", max_index + 1
            print amino_acid_tally[1][0], ",", amino_acid_tally[1][1], ",", residue_types[2], ", min ,", min(normalised_values), ",position,", min_index + 1

        # for residue_types

for file in list_of_files:
    print file
    for type_of_protein in type_of_proteins:
        print type_of_protein
        list_to_print = []
        results = []

        with open(file) as inputfile:
            for line in inputfile:
                results.append(line.strip().split(','))
        # The total list of residues is held in a list. This will be turned
        # into a string for counting later.
        list_of_slices = []

        # Values for longest TMD, and everything can be aligned by correcting
        # to this number
        halfway_value_for_alignment = 0

        for entry in results:
            # Don't read header line
            if entry == results[0]:
                pass
            else:
                #We are assuming the database topology is correct.

                total_tmd_count = int(1)
                n_terminal_start = str("Inside")
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

                # Check that it mathces the type
                if single_or_multi(type_of_protein, total_tmd_count, id) == True:
                    halfway_value_for_alignment = compensation_checker(
                        tmh_sequence, n_flank_sequence, n_terminal_start, halfway_value_for_alignment)

        # halfway_value_for_alignment now holds the number needed to align to 0

        # this will be filled with lists of [file, type_of_protein, amino_acid,
        # sum(sequence_residue_count), sequence_residue_count ]
        amino_acid_tally = []

        for amino_acid in residues:
            record_count = 0
            sequence_residue_count = [0]
            #+1 avoids issues with maximum lengths being odd resulting in shorter list lengths than positions in the sequence.
            for n in range(0, 2 * halfway_value_for_alignment + 1):
                sequence_residue_count.append(0)
            for entry in results:
                # Don't read header line
                if entry == results[0]:
                    pass
                else:
                    total_tmd_count = int(1)
                    n_terminal_start = str("Inside")
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

                    if single_or_multi(type_of_protein, total_tmd_count, id) == True:
                        record_count = record_count + 1
                        sequence_to_measure = oriented_sequence(
                            n_terminal_start, n_flank_sequence, tmh_sequence, c_flank_sequence, halfway_value_for_alignment)

                        for position, an_amino_acid in enumerate(sequence_to_measure):
                            if an_amino_acid == amino_acid:
                                sequence_residue_count[position] = int(
                                    sequence_residue_count[position]) + 1

            absolute_vectors = [amino_acid, sequence_residue_count]

            # This is the normalisation formula.
            normalised_vectors = []
            for positional_frequency in absolute_vectors[1]:
                if sum(sequence_residue_count) > 0:
                    normalised_value = (100 * positional_frequency) / \
                        sum(sequence_residue_count)
                    normalised_vectors.append(normalised_value)
                else:
                    normalised_vectors.append(0)

            # This is calculating the averages across different sections.
            # Like the background level...
            values_for_noise_inside = normalised_vectors[
                int(halfway_value_for_alignment - 30):int(halfway_value_for_alignment - 25)]
            values_for_noise_outside = normalised_vectors[
                int(halfway_value_for_alignment + 25):int(halfway_value_for_alignment + 30)]

            # And the flank averages.
            values_for_flank_inside = normalised_vectors[
                int(halfway_value_for_alignment - 20):int(halfway_value_for_alignment - 10)]
            values_for_flank_outside = normalised_vectors[
                int(halfway_value_for_alignment + 10):int(halfway_value_for_alignment + 20)]

            mean_flank_inside = np.mean(values_for_flank_inside)
            mean_flank_outside = np.mean(values_for_flank_outside)

            # All the information is added to the tally.
            for_tallying = [file, type_of_protein, amino_acid, sum(
                sequence_residue_count), sequence_residue_count, normalised_vectors, values_for_noise_inside, values_for_noise_outside, mean_flank_inside, mean_flank_outside, record_count]
            amino_acid_tally.append(for_tallying)

        # Amino acid tally now contains the following information for each amino acid type in a given dataset:
        # 0 file name
        # 1 type_of_protein (multipass or singlepass)
        # 2 amino_acid type
        # 3 total number of the amino acid
        # 4 how many of the amino acid at each position
        # 5 the amino acid normalised frequency at each position
        # 6 Background frequency inside,
        # 7 Background frequency outside ,
        # 8 Average values for flanks inside,
        # 9 Average values for flanks outside

        # For Figure 1 Uncomment these lines to get the total values of each
        # amino acid printed in a table.
        print "\nFIGURE 1"
        total_amino_acids = 0
        print amino_acid_tally[1][0], amino_acid_tally[1][1], amino_acid_tally[1][10], "records."
        print "Total amino acids,", total_amino_acids, "\n\n"
        for i in amino_acid_tally:
            print i[2], i[3]
            total_amino_acids = total_amino_acids + i[3]

        # Uncomment these lines to get the values from -30 to +30

        print "\nABSOLUTE VECTORS for -30 to +30"
        vectors_for_neg30_to_pos30(
            amino_acid_tally, halfway_value_for_alignment)

        print "\nNORMALISED VECTORS for -30 to +30"
        normalised_vectors_for_neg30_to_pos30(
            amino_acid_tally, halfway_value_for_alignment)

        print "\nABSOLUTE VECTORS"
        print amino_acid_tally[0][0], amino_acid_tally[0][1]
        for i in amino_acid_tally:
            print i[2], i[4]

        print "\nNORMALISED VECTORS"

        for i in amino_acid_tally:
            print i[2], i[5]

        # Uncomment for figure showing absolute vectors and stats comparison of
        # inside versus outside.
        print "\nABSOLUTE VECTORS", amino_acid_tally[1][0],
        amino_acid_tally[1][1]

        # P-value for a KW test. This may run into errors for UniArch. If you
        # are only interested in replicating the data for tables 2 and 3, see
        # the table 2.py scripts."
        '''
        print "\nTable for global average statistics\n"
        print "type_of_protein, sum inside, sum outside, test stat"

        for i in amino_acid_tally:
            if len(i[4])>61:
                absolute_vector = i[4][halfway_value_for_alignment - 30:halfway_value_for_alignment + 31]
                inside_flank = absolute_vector[10:20]
                outside_flank = absolute_vector[40:50]
                stat_test = scipy.stats.kruskal(inside_flank, outside_flank)

                print i[2], ",", sum(inside_flank), ",", sum(outside_flank), ",", stat_test[0], ",", stat_test[1]
        '''
        # For Figure 2
        # print "\nFIGURE 2"

        # First we need to stylise the graph
        ax = plt.subplot(111)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        # The next block of code generates the gray boundaries.
        '''
        # These boxes help us see where the membrane is by drawing in the flank
        # areas used to define the averages, as well as showing the background
        # noise areas.
        ax.add_patch(patches.Rectangle(
            (-10, 0),   # (x,y)
            20,          # width
            9999999999999999,          # height
            alpha=0.01
        )
        )
        ax.add_patch(patches.Rectangle(
            (-30, 0),   # (x,y)
            5,          # width
            9999999999999999,          # height
            alpha=0.1,
            facecolor="gray"
        )
        )
        ax.add_patch(patches.Rectangle(
            (25, 0),   # (x,y)
            5,          # width
            9999999999999999,          # height
            alpha=0.1,
            facecolor="gray"
        )
        )

        ax.add_patch(patches.Rectangle(
            (-20, 0),   # (x,y)
            10,          # width
            9999999999999999,          # height
            alpha=0.1,
            facecolor="gray"
        )
        )
        ax.add_patch(patches.Rectangle(
            (10, 0),   # (x,y)
            10,          # width
            9999999999999999,          # height
            alpha=0.1,
            facecolor="gray"
        )
        )
        plt.plot((0, 0), (0, 9999999999999999), linestyle='-', linewidth=1, color='black')
        '''
        # To maintain consistency of the figures throughout the paper, the
        # lines for key residues are specified and coloured accordingly.

        # These functions describe the difference lines that will be drawn.

        #
        def background_level(background_value, colour):
            '''
            Returns the values for the hashed line for the background value.
            '''
            max_sequence_length = len(background_value)
            values_for_noise_inside = background_value[
                int((max_sequence_length / 2) - 30):int((max_sequence_length / 2) - 25)]
            values_for_noise_outside = background_value[
                int((max_sequence_length / 2) + 25):int((max_sequence_length / 2) + 30)]
            noise_inside = np.mean(values_for_noise_inside)
            noise_outside = np.mean(values_for_noise_outside)
            noise_normalised_value = np.mean([noise_inside, noise_outside])

            background_value = noise_normalised_value
            plt.plot((-50, +50), (background_value, background_value),
                     linestyle='--', color=colour)

        # This plots the flank length averages inside and outside.
        def flank_levels(inside, outside, colour):
            plt.plot((-20, -10), (inside,
                                  inside), linestyle='-', linewidth=4,  color=colour)
            plt.plot((10, 20), (outside,
                                outside), linestyle='-', linewidth=4,  color=colour)

        # This plots the frequency line.
        def solid_line_values(vectors, colour):
            # First we need the position co-ordinates
            sequence_position = []
            for number, item in enumerate(range(0, len(vectors))):
                value_for_position = number - halfway_value_for_alignment
                sequence_position.append(value_for_position)
            plt.plot(sequence_position, vectors,
                     linestyle='-', marker='',  linewidth=2,  color=colour)

        def scatter_line_values(vectors, colour):
            # First we need the position co-ordinates
            sequence_position = []
            for number, item in enumerate(range(0, len(vectors))):
                value_for_position = number - halfway_value_for_alignment
                sequence_position.append(value_for_position)
            plt.plot(sequence_position, vectors,
                     linestyle='', marker='.',  linewidth=2,  color=colour)

        # For absolute frequencies change residue_types[5] to residue_types[4]
        # in the following lines and comment out the flank levels.
        for residue_types in amino_acid_tally:
            if residue_types[2] == "L":
                #background_level(residue_types[5], "blue")
                #solid_line_values(residue_types[5], "blue")

                solid_line_values(residue_types[4], "blue")

            if residue_types[2] == "K":
                #background_level(residue_types[5], "peachpuff")
                #flank_levels(residue_types[8], residue_types[9], "peachpuff")
                #solid_line_values(residue_types[5], "peachpuff")
                #scatter_line_values(residue_types[5], "peachpuff")

                solid_line_values(residue_types[4], "peachpuff")
                #scatter_line_values(residue_types[4], "peachpuff")

            if residue_types[2] == "R":
                #background_level(residue_types[5], "orange")
                #flank_levels(residue_types[8], residue_types[9], "orange")
                #solid_line_values(residue_types[5], "orange")
                #scatter_line_values(residue_types[5], "orange")

                solid_line_values(residue_types[4], "orange")
                #scatter_line_values(residue_types[4], "orange")

            if residue_types[2] == "D":
                #background_level(residue_types[5], "purple")
                #flank_levels(residue_types[8], residue_types[9], "purple")
                #solid_line_values(residue_types[5], "purple")
                #scatter_line_values(residue_types[5], "purple")

                solid_line_values(residue_types[4], "purple")

            if residue_types[2] == "E":
                #background_level(residue_types[5], "thistle")
                #flank_levels(residue_types[8], residue_types[9], "thistle")
                #solid_line_values(residue_types[5], "thistle")
                #scatter_line_values(residue_types[5], "thistle")

                solid_line_values(residue_types[4], "thistle")

        font = {'size': 18}

        plt.xlabel('Sequence Position', **font)
        plt.ylabel('Relative Percentage', **font)
        plt.tick_params(labelsize=16)
        pylab.xlim([-30, 30])
        # This sets the height of the graph - useful for the relative percentage figures.
        # pylab.ylim([0, 7])
        timestamp = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        filename = timestamp + \
            file.replace('.csv', '') + amino_acid_tally[1][1] + ".pdf"
        plt.gcf().subplots_adjust(bottom=0.2)
        plt.savefig(filename)
        # plt.show()

        plt.clf()
        plt.cla()

        # need the averages and max values.
