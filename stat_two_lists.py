import numpy as np
import scipy
from scipy import stats
import sys


def bahadur(list1, list2, pvalue):
    '''
    Takes two lists and a p-value to calculate the Bahadur slope.
    The Bahadur value is calculated from the modulus of natural log of the P-value divided by the sample size powering the test.
    '''
    bahadur_vlaue = (abs(0 - np.log(pvalue))) / (len(list1) + len(list2))
    return bahadur_vlaue


def stats(factor, list1, list2):
    '''Prints stats for two lists.'''
    print(factor, ", Kruskal-Wallis,", scipy.stats.kruskal(list1, list2)[0], ",", scipy.stats.kruskal(
        list1, list2)[1], ",", bahadur(list1, list2, scipy.stats.kruskal(list1, list2)[1]))
    print(factor, ", Kolmogorov-Smirnov, ", scipy.stats.ks_2samp(list1, list2)[0], ",", scipy.stats.ks_2samp(
        list1, list2)[1], ",", bahadur(list1, list2, scipy.stats.ks_2samp(list1, list2)[1]))
    print(factor, ", Student's T-test, ", scipy.stats.ttest_ind(list1, list2)[0], ",", scipy.stats.ttest_ind(
        list1, list2)[1], ",", bahadur(list1, list2, scipy.stats.ttest_ind(list1, list2)[1]))
    # print(factor, ", Chi Squre test, ", scipy.stats.chisquare(list1, list2)[0], ",", scipy.stats.chisquare(
    #    list1, list2)[1], ",", bahadur(list1, list2, scipy.stats.chisquare(list1, list2)[1]))


print("Factor,Test, Test-statistic, P value, Bahadur slope")

list_of_files = sys.argv[1:]
disorder_sets = []
flanks_disorder_sets = []

hydrophobicity_sets = []
flanks_hydrophobicity_sets = []

entropy_sets = []
flanks_entropy_sets = []

for file_number, file in enumerate(list_of_files):
      # This generates an empty list for each potential possition. Values of hydrophobicity will be added to this later and will contribute to the average.
    results = []
    disorder_set = []
    hydrophobicity_set = []
    flanks_disorder_set = []
    flanks_hydrophobicity_set = []
    entropy_set = []
    flanks_entropy_set = []

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
            tmh_hydrophobicity = float(entry[9])
            tmh_flank_hydrophobicity = float(entry[10])
            tmh_disorder = float(entry[11])
            tmh_flank_disorder = float(entry[12])
            tmh_entropy = float(entry[13])
            tmh_flank_entropy = float(entry[14])

            hydrophobicity_set.append(tmh_hydrophobicity)
            flanks_hydrophobicity_set.append(tmh_flank_hydrophobicity)
            disorder_set.append(tmh_disorder)
            flanks_disorder_set.append(tmh_flank_disorder)
            entropy_set.append(tmh_entropy)
            flanks_entropy_set.append(tmh_flank_entropy)
    hydrophobicity_sets.append(hydrophobicity_set)
    flanks_hydrophobicity_sets.append(flanks_hydrophobicity_set)
    disorder_sets.append(disorder_set)
    flanks_disorder_sets.append(flanks_disorder_set)
    entropy_sets.append(entropy_set)
    flanks_entropy_sets.append(flanks_entropy_set)

print("Mean Values And Error")
for n, dataset in enumerate(list_of_files):
    print(list_of_files[n])
    print("Factor, TMH average, TMH and flanks average, TMH standard deviation, TMH and flanks standard deviation")
    print("Hydrophobicity,", np.mean(hydrophobicity_sets[n]),",",np.mean(flanks_hydrophobicity_sets[n]), ",", np.std(hydrophobicity_sets[n]),",",np.std(flanks_hydrophobicity_sets[n]))
    print("Disorder,", np.mean(disorder_sets[n]), ",", np.mean(flanks_disorder_sets[n]), ",", np.std(disorder_sets[n]),",",np.std(flanks_disorder_sets[n]))
    print("Sequence Entropy,", np.mean(entropy_sets[n]), ",", np.mean(flanks_entropy_sets[n]), ",", np.std(entropy_sets[n]),",",np.std(flanks_entropy_sets[n]))



print("Statistical table")
#print("\n\nHydrophobicity of TMH\n")
stats(str("Hydrophobicity of TMH"),
      hydrophobicity_sets[0], hydrophobicity_sets[1])

#print("\n\nHydrophobicity of TMH and flanks\n")
stats(str("Hydrophobicity of TMH and flanks"),
      flanks_hydrophobicity_sets[0], flanks_hydrophobicity_sets[1])


#print("\n\nDisorder of TMH\n")
stats(str("Disorder of TMH"), disorder_sets[0], disorder_sets[1])

#print("\n\nDisorder of TMH and flanks\n")
stats(str("Disorder of TMH and flanks"),
      flanks_disorder_sets[0], flanks_disorder_sets[1])

#print("\n\nentropy of TMH\n")
stats(str("Sequence Entropy of TMH"), entropy_sets[0], entropy_sets[1])

#print("\n\nentropy of TMH and flanks\n")
stats(str("Sequence Entropy of TMH and flanks"),
      flanks_entropy_sets[0], flanks_entropy_sets[1])
