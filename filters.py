import sys
from Bio import SeqIO
import numpy as np


print("This script was developed by James A Baker under the supervision of Dr Jim Warwicker.\nIt requires an active internet connection and an up to date version of biopython.\nSee readme.md for more information on installation and visit www.github.com/jbkr/TA_predict to report any errors.")

input_file = str(sys.argv[1])
flank_length = 5
# This works with uniprot filetype. From the seqIO biopython wiki: Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also the UniProt XML format if something goes wrong.
input_format = "swiss"
# For future modification, this can be used to look for any annotation in the file.
feature_type = "TRANSMEM"
# Simply the output name, can be anything as it is written in binary (not file-type specific language).
output_filename_fasta = "TMD.fasta"
other_feature_type = "NON-TER"
signal_feature = "SIGNAL"
subcellular_location = "TOPO_DOM"
minimum_tmd_length = 16
maximum_tmd_length = 38
flank_clash_amendment_status = True
length_excluded_tmds = []


# For each file, a table is generated for each of the flank lengths set.


# Single-pass and mult-pass are treated differently because at the end
# we report on the numbers so we avoid re-iterating through the larger
# datasets, and throughout the process they are handled differently for
# example when looking for clashes of flanking regions.
number_of_records = 0
number_of_records_correct_length = 0
number_of_records_single = 0
number_of_records_correct_length_single = 0

# File output name
output_filename = input_file.replace(".txt", ".csv")

# The header row in the file.
with open(output_filename, 'w') as my_file:
    my_file.write("Name and description, ID, N terminal inside/outside, tmh start location, tmh end location, full protein sequence, tmh sequence, N flank sequence, C flank sequence, transmembrane helix sequential number, number of transmembrane helices in protein\n")
my_file.closed

# We need to check against nearby features to prevent overlapping
# flanking regions. Note here we want to avoid clashing with INTRAMEM
# regions, since their flanking regions may be similar, however INTRAMEM regions
# should not be included in the logged transmembrane features since
# they may have extremely short transmembrane sequences which would
# disrupt very much so the alignment of the flanking regions.
avoid_features = ["TRANSMEM", "INTRAMEM"]
unknown = 0
reliable_flank_length = flank_length

# We iterate through each record, parsed by biopython.
for record in SeqIO.parse(input_file, input_format):
    new_record = True
    tmd_count = 0
    for i, f in enumerate(record.features):
        if f.type == feature_type:

            n_terminal_start = "None"
            tmh_record = []
            tmd_count = tmd_count + 1

            # Calls the parsing functions as strings for later use.
            id_of_record = record.id
            name_of_record = str(
                record.description).replace(",", "")

            # The first separation of single or multi-pass is based on
            # annotation itself, however this is prone to error, so is
            # checked later.
            if "; Single-pass" in str(record):
                single_or_multi_topology = "Single-pass"
            if "; Multi-pass" in str(record):
                single_or_multi_topology = "Multi-pass"

            # Some transmembrane annotations have unknown sequence
            # positions where it is ambiguous. These features are
            # discounted.
            if "UnknownPosition" in str(f.location):
                pass
                unknown = unknown + 1
                print id_of_record, "had an unknown position TMH."
            else:

                full_sequence = str(record.seq)
                tmh_start = int(f.location.start)
                tmh_stop = int(f.location.end)
                tmh_sequence = str(
                    record.seq[(f.location.start):(f.location.end)])

                # If a clash is detected in either flank, the
                # script will compensate for the clash and stop
                # checking to see if there are additional
                # clashes at that flank. This is reset per TMH.

                flank1_length = reliable_flank_length
                flank2_length = reliable_flank_length

                # We iterate through each feature in the record, to
                # assertain if there are any clashes.
                for each_features in record.features:
                    # Another check for unknown positions to rule
                    # out.
                    if "UnknownPosition" in str(each_features.location):
                        pass
                    else:
                        #  We identify if any of these features are within the flanking residue distance.
                        # If they are clashes, then the flank length is half the distance between the features.
                        # Note now, we only need to worry about the flank length of the logged feature, since the other
                        # flank length will be dealt with when
                        # iterating through the features looking for TM
                        # regions.

                        # Let's first check if the N terminus of this
                        # feature comes after the C terminus of the
                        # feature we're logging.
                        # if C_terminal_flank_clash is not True:
                        if each_features.location.start >= f.location.end:

                            # Is the feature and the N terminal flank starting within the flank
                            # of the C terminus of the feature that
                            # we're logging?

                            if (each_features.location.start - reliable_flank_length) <= (f.location.end + reliable_flank_length)and flank_clash_amendment_status == True:
                                if str(each_features.type) in avoid_features:

                                    # There will be a clash/overlapping
                                    # of transmembrane flank lengths at
                                    # the C flank.

                                    max_flank_size = each_features.location.start - f.location.end
                                    new_flank2_length = max_flank_size / 2

                                    if new_flank2_length < flank2_length:
                                        flank2_length = new_flank2_length

                                    if each_features.type == "INTRAMEM":
                                        print_intramem_details = True

                                else:
                                    # There is a feature, but it's not
                                    # a transmembrane region, so
                                    # flanking regions won't be double
                                    # counted
                                    pass
                            else:
                                pass
                        else:
                            pass

                        # Now we will check the same as above, however
                        # for features that the C terminus falls before
                        # the feature being logged N-terminus
                        if each_features.location.end <= f.location.start:
                            if (each_features.location.end + reliable_flank_length) >= (f.location.start - reliable_flank_length) and flank_clash_amendment_status == True:
                                if str(each_features.type) in avoid_features:

                                    max_flank_size = f.location.start - each_features.location.end
                                    new_flank1_length = max_flank_size / 2

                                    if new_flank1_length < flank1_length:
                                        flank1_length = new_flank1_length

                                    if each_features.type == "INTRAMEM":
                                        print_intramem_details = True
                                else:
                                    pass
                            else:
                                pass
                        else:
                            pass

                        # Now we will amend the flank length incase it
                        # exceeds the maximum protein length, or the
                        # 0th sequence item (the start of the
                        # sequence)
                        if (f.location.start - int(flank1_length)) >= 0:
                            n_terminal_flank = str(
                                record.seq[(f.location.start - int(flank1_length)):(f.location.start)])
                        elif (f.location.start - int(flank1_length)) < 0:
                            n_terminal_flank = str(
                                record.seq[0:(f.location.start)])

                        if (f.location.end + int(flank2_length)) > len(record.seq):
                            c_terminal_flank = str(
                                record.seq[(f.location.end):(int(len(record.seq)))])
                        else:
                            c_terminal_flank = str(
                                record.seq[(f.location.end):(f.location.end + int(flank2_length))])

                # Now, the orientation is determined for all features
                # according to the cytoplasmic annotation. The total
                # number of cytoplasmic annotations, and non
                # cytoplasmic annotations are counted, and depending on
                # the starting locations falling "inside" or "outside" the cytoplasm,
                # the orientation is determined.

                list_of_cyto_starts = []
                list_of_non_cyto_starts = []

                for feature_number, other_features in enumerate(record.features):

                    # Some features contain unknown locations. These are
                    # discarded.
                    if "UnknownPosition" in str(other_features.location):
                        pass
                    else:
                        # If the feature matches the cytoplasm...
                        if other_features.type == subcellular_location:
                            if "Cytoplasmic" in str(other_features.qualifiers):

                                # The feature start locations are
                                # recorded as being cytoplasmic or not.
                                # This, combined with knowing the start
                                # location, allows us to go through the
                                # TMHs and figure out which way each is
                                # oriented quickly.
                                list_of_cyto_starts.append(
                                    int(other_features.location.start))

                                total_tmd_count = 0
                                for number_of_features, a_feature in enumerate(record.features):

                                    if a_feature.type == feature_type:
                                        total_tmd_count = total_tmd_count + 1

                                new_record = False

                            else:
                                list_of_non_cyto_starts.append(
                                    int(other_features.location.start))
                        else:
                            pass

                # We will now check if the residues preceding the
                # TMD are intra, or extra cytoplasmic. Because this
                # calls locations below the feature location start,
                # it must be called only when above 1. If this
                # cannot be done, the topology can be infered in
                # the next script.
                if n_terminal_start == "none" and f.location.start > 1:
                    previous_feautre_location = f.location.start - 1
                    for index, a_features in enumerate(record.features):
                        if a_features.type == subcellular_location and a_features.location.start < previous_feautre_location and a_features.location.end > previous_feautre_location:
                            if "Cytoplasmic" in str(a_features.qualifiers):
                                n_terminal_start == "Inside"
                            else:
                                n_terminal_start = "Outside"

                # If the previous method did not work to identify topology, we can still
                # imply the topology.
                if n_terminal_start == "None":
                    # This checks that the subcellular starting
                    # locations exist.
                    if len(list_of_cyto_starts) > 0 and len(list_of_non_cyto_starts) > 0:

                        # We need to furthermore check whether the
                        # starting location was inside or outside.
                        if min(list_of_non_cyto_starts) < min(list_of_cyto_starts):

                            # The modulo is used to check if the
                            # current tmd number is even or odd. From
                            # this and the first TMD we know whether
                            # this TMD starts inside or outside.
                            if tmd_count % 2 == 0:
                                n_terminal_start = "Inside"
                            elif tmd_count % 2 == 1:
                                n_terminal_start = "Outside"

                        elif min(list_of_cyto_starts) < min(list_of_non_cyto_starts):
                            if tmd_count % 2 == 0:
                                n_terminal_start = "Outside"
                            elif tmd_count % 2 == 1:
                                n_terminal_start = "Inside"
                        else:
                            pass

                # Checks if the single_or_multi_topology variable
                # exists (it would have been created earlier)
                if 'single_or_multi_topology' not in locals() or 'single_or_multi_topology' not in globals():
                    if total_tmd_count == 1:
                        single_or_multi_topology = "Single-pass"
                    if total_tmd_count > 1:
                        single_or_multi_topology = "Multi-pass"

                if "Inside" in n_terminal_start or "Outside" in n_terminal_start:
                    # This is the information that will be written for the record.
                    #+/-1s are used since slices originally call how many steps to iterate rather than the sequence postion. This matches the Uniprot sequence numbering
                    tmh_record = [name_of_record, id_of_record, tmh_start + 1, tmh_stop, abs(tmh_start - tmh_stop) - 1,
                                  full_sequence, tmh_sequence, n_terminal_flank, c_terminal_flank]

                    number_of_records = number_of_records + 1

                    if total_tmd_count == 1:
                        number_of_records_single = number_of_records_single + 1

                    # Pythonic language would usually concatenate all conditions onto one line, but for the sake of clarity, here we do it statement by statement.
                    if len(tmh_sequence) >= minimum_tmd_length:
                        if len(tmh_sequence) <= maximum_tmd_length:
                            number_of_records_correct_length = number_of_records_correct_length + 1

                            # Is this a single-pass protein?
                            if total_tmd_count == 1:
                                number_of_records_correct_length_single = number_of_records_correct_length_single + 1

                                # Is the N-terminal cytoplasmic?
                                if "Inside" in n_terminal_start:
                                    print str(record.seq)
                                    # Is the C-terminal within 25 residues of the final residue?
                                    if abs(tmh_stop - len(str(record.seq))) < 26:
                                        with open(output_filename, 'a') as my_file:
                                            for i in tmh_record:
                                                my_file.write(str(i))
                                                my_file.write(",")
                                            my_file.write("\n")

                    else:
                        length_exclusion_info = str(
                            id_of_record) + "_" + str(tmd_count)
                        length_excluded_tmds.append(
                            length_exclusion_info)

                # No records should be here with none. This is for
                # debugging only.
                elif "None" in n_terminal_start:
                    pass
                else:
                    pass
# General information regarding the output file useful as a log.
print input_file, ", flank_clash_amendment_status:", flank_clash_amendment_status
print "Number of TMHs in dataset:", number_of_records
print "Number of TMHs after dumping incorrect lengths:", number_of_records_correct_length
print "Records"
print number_of_records
print number_of_records_correct_length
print "Single-pass"
print "Total:", number_of_records_single
print "After length exclusion:", number_of_records_correct_length_single

exclusion_ids_output_filename = input_file.replace(
    ".txt", "_%s_flanklength_flankclash%s_logged_lengthexclusionIDs.txt" % (flank_length, str(flank_clash_amendment_status)))
with open(exclusion_ids_output_filename, 'a') as my_file:
    my_file.write("# ID_HelixNumber\n")
my_file.closed
for item in length_excluded_tmds:
    with open(exclusion_ids_output_filename, 'a') as my_file:
        my_file.write(item)
        my_file.write("\n")
my_file.closed