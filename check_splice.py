from Bio import SeqIO
import os
import os.path
import re
import fileinput
import sys
import urllib,urllib2
import subprocess
import shutil


with open('C_terminal_single_TRANSMEM.txt') as custom_list:
    custom_list = custom_list.readlines()
    interaction = list(custom_list)

for i in interaction: #Begins the iteration
    i = i.replace(" \n", "")
    i = i.replace("\n", "")
    i = i.replace(" ", "") #Removes the spaces between lines. This was causing some really weird bugs and cutting the url below in half.
    urllib.urlretrieve("http://www.uniprot.org/uniprot/%s.txt" % i, filename='uniget.dat') #This uses the ID (saved as i) in a file called uniget.dat.
    #print "Getting %s data from the uniprot." % i

    #The individual result held in uniget.dat is saved as lildat, and then added to bigdat.
    with open ("uniget.dat", "r") as lildatfile:
        lildat=lildatfile.read()
        #Here we use biopython to open the data downloaded from uniprot and check the distance from the C terminus to the TMD.
        filenames = ["uniget.dat"]
        input_format = "swiss" #This SHOULD work with uniprot filetype. From the seqIO biopython wiki: Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also the UniProt XML format.
        feature_type = "NON-TER" #For future modification, this can be used to look for any annotation in the .dat file.
        output_filename = "transmembranes.fasta" #Simply the output name, can be anything as it is written in binary (not file-type specific language).

        output = open(output_filename, "w")
        for filename in filenames:
            # Using SeqIO.parse will cope with multi-record files
            for record in SeqIO.parse(filename, input_format): #A biopython module that can automatically parse uniprot .txt files
                for f in record.features:
            #feature type refers to TRANSMEM
                    with open ("splice_variant.txt", "ab") as single_fasta_file:
                        if f.type == feature_type:
                            print "Splice varient detected in ",record.id, ". ", record.id, " will be removed from the final list."


                #This prints the title and then any TRANSMEM regions associated with that entry along with flanking regions.

                            single_fasta_file.write(record.id)
                            single_fasta_file.write("\n")
                        else:
                            pass
                            #print "No splice isoform detected in this feature of ", record.id

#First file open

#print "These are the total splice varients:\n", f.read()
with open('splice_variant.txt') as file1list:
    spliced = file1list.readlines()

#opens second file to compare

with open('C_terminal_single_TRANSMEM.txt') as file2list:
    C_terminal = file2list.readlines()

#Compares the intersections of the lists i.e finds the matches
if len(spliced)>1:
    output = set(C_terminal).remove(set(spliced))
else:
    pass

output = str(output)
for char in "'][,\)(nset":
    output = output.replace(char,'')



#reports hits
#print "The following protein IDs are hits: \n \n \n", output, "\n \n"

#writes hits to file
#print 'Printing hits to output file! ...in progress... \n'
file = open("not_spliced_C_terminal_single_TRANSMEM.txt", "w")
file.write(str(output))
file.close()
