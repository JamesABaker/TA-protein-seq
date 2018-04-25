from __future__ import division
import sys
from Bio import SeqIO
import numpy as np
import subprocess
import re

input_file = str(sys.argv[1])
flank_length = 5
# This works with uniprot filetype. From the seqIO biopython wiki:
# Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also
# the UniProt XML format if something goes wrong.
input_format = "swiss"
# For future modification, this can be used to look for any annotation in
# the file.
feature_type = "TRANSMEM"
# Simply the output name, can be anything as it is written in binary (not
# file-type specific language).
output_filename_fasta = str("TMD"+str(input_file)+".fasta")

# For each file, a table is generated for each of the flank lengths set.

# File output name
output_filename = input_file.replace(".txt", ".csv")

# The header row in the file.
with open(output_filename, 'w') as my_file:
    my_file.write("Name and description, ID, tmh start location, tmh end location, tmh length, full protein sequence, tmh sequence, N flank sequence, C flank sequence\n")
my_file.closed

# We need to check against nearby features to prevent overlapping
# flanking regions. Note here we want to avoid clashing with INTRAMEM
# regions, since their flanking regions may be similar, however INTRAMEM regions
# should not be included in the logged transmembrane features since
# they may have extremely short transmembrane sequences which would
# disrupt very much so the alignment of the flanking regions.
unknown = 0
reliable_flank_length = flank_length


# We iterate through each record, parsed by biopython.
for record in SeqIO.parse(input_file, input_format):
    for i, f in enumerate(record.features):
        if f.type == feature_type:
            # Calls the parsing functions as strings for later use.
            id_of_record = record.id
            name_of_record = str(
                record.description).replace(",", "")

            # Some transmembrane annotations have unknown sequence
            # positions where it is ambiguous. These features are
            # discounted.
            if "UnknownPosition" in str(f.location):
                pass
                unknown = unknown + 1
                print(id_of_record, "had an unknown position TMH.")
            else:

                full_sequence = str(record.seq)
                tmh_start = int(f.location.start)
                tmh_stop = int(f.location.end)
                tmh_sequence = str(record.seq[(f.location.start):(f.location.end)])
                #These are not c or n terminal for sure. This is just an assumption we make since we are not formally filtering anything in this list.
                n_terminal_flank = record.seq[(f.location.start+1-5):(f.location.start)]
                c_terminal_flank = record.seq[(f.location.end):(f.location.end+5)]

                # This is the information that will be written for the record.
                # +/-1s are used since slices originally call how many steps to iterate rather than the sequence postion. This matches the Uniprot sequence numbering
                tmh_record = [name_of_record, id_of_record, tmh_start + 1, tmh_stop, abs(tmh_start - tmh_stop) - 1,
                            full_sequence, tmh_sequence, n_terminal_flank, c_terminal_flank]

                with open(output_filename, 'a') as my_file:
                    for i in tmh_record:
                        my_file.write(str(i))
                        my_file.write(",")
                    my_file.write("\n")
                with open(output_filename_fasta, 'a') as filtered_fasta_file:
                    #This prevents several Fasta entries for the same record if splice isoforms exist.
                    fasta_written=True
                    fasta_record=str(">"+str(record.id)+"\n"+str(record.seq)+"\n")
                    filtered_fasta_file.write(fasta_record)
