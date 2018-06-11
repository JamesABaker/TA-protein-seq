import sys

list_of_files = sys.argv[1:]


def listify(file):
    with open(file) as f:
        lines = f.read().splitlines()
    return(lines)

for number_of_list, list in enumerate(list_of_files):
    for number_of_comparison_list, comparison_list in enumerate(list_of_files):
        if str(list) == str(comparison_list):
            pass
        else:
            number_of_matches = 0
            number_of_misses = 0

            for id in listify(list):
                for comparison_id in listify(comparison_list):
                    if str(id) == str(comparison_id):
                        number_of_matches = number_of_matches + 1
                    else:
                        number_of_misses = number_of_misses + 1

            print(list, "against", comparison_list)
            print("list_size, ", len(listify(list)))
            print("matches, ", number_of_matches)
