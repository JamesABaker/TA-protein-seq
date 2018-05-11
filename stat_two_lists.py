import numpy as np
import scipy
from scipy import stats
import sys

def stats(list1, list2):
    '''Stats on two lists'''
    print(scipy.stats.kruskal(list1, list2))
    print(scipy.stats.ks_2samp(list1, list2))
    print(scipy.stats.ttest_ind(list1, list2))

list_of_files = sys.argv[1:]
disorder_sets = []
hydrophobicity_sets = []
flanks_disorder_sets = []
flanks_hydrophobicity_sets = []

for file_number, file in enumerate(list_of_files):
      # This generates an empty list for each potential possition. Values of hydrophobicity will be added to this later and will contribute to the average.
    results = []
    disorder_set = []
    hydrophobicity_set = []
    flanks_disorder_set = []
    flanks_hydrophobicity_set = []
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

          disorder_set.append(tmh_disorder)
          hydrophobicity_set.append(tmh_hydrophobicity)
          flanks_disorder_set.append(tmh_flank_disorder)
          flanks_hydrophobicity_set.append(tmh_flank_hydrophobicity)
    disorder_sets.append(disorder_set)
    hydrophobicity_sets.append(hydrophobicity_set)
    flanks_disorder_sets.append(flanks_disorder_set)
    flanks_hydrophobicity_sets.append(flanks_hydrophobicity_set)


print("\n\nDisorder of TMH\n")
stats(disorder_sets[0], disorder_sets[1])

print("\n\nHydrophobicity of TMH\n")
stats(hydrophobicity_sets[0], hydrophobicity_sets[1])

print("\n\nDisorder of TMH and flanks\n")
stats(flanks_disorder_sets[0], hydrophobicity_sets[1])

print("\n\nHydrophobicity of TMH and flanks\n")
stats(flanks_hydrophobicity_sets[0], hydrophobicity_sets[1])
