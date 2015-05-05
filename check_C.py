from Bio import SeqIO
import os
import os.path
import re
import fileinput
import sys
import urllib,urllib2
import subprocess
import shutil


with open('single_pass_list.txt') as custom_list:
    custom_list = custom_list.readlines()
    protein_list = list(custom_list)

for i in protein_list: #Begins the iteration
    i = i.replace(" \n", "")
    i = i.replace("\n", "")#Removes the spaces between lines. This was causing some really weird bugs and cutting the url below in half.
    urllib.urlretrieve("http://www.uniprot.org/uniprot/%s.txt" % i, filename='uniget.dat') #This uses the ID (saved as i) in a file called uniget.dat.
    #print "Getting %s data from the uniprot." % i
    #print "http://www.uniprot.org/uniprot/%s.txt" % i


    #The individual result held in uniget.dat is saved as lildat, and then added to bigdat.
    with open ("uniget.dat", "r") as lildatfile:
        lildat=lildatfile.read()
        #Here we use biopython to open the data downloaded from uniprot and check the distance from the C terminus to the TMD.
        filenames = ["uniget.dat"]
        input_format = "swiss" #This SHOULD work with uniprot filetype. From the seqIO biopython wiki: Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also the UniProt XML format.
        feature_type = "TRANSMEM" #For future modification, this can be used to look for any annotation in the .dat file.
        output_filename = "transmembranes.fasta" #Simply the output name, can be anything as it is written in binary (not file-type specific language).

        output = open(output_filename, "w")
        for filename in filenames:
    # Using SeqIO.parse will cope with multi-record files
            for record in SeqIO.parse(filename, input_format): #A biopython module that can automatically parse uniprot .txt files
                for f in record.features:
            #feature type refers to TRANSMEM
                    if f.type == feature_type:
                        C_terminal_distance = len(record.seq) - f.location.end
                        #print "The TMD of", record.id, " is ", C_terminal_distance, " residues from the C terminal."
                        if C_terminal_distance < 26 :
                #These are the flanking regions and the TMD using record.seq. They are identified by simply counting 5 spaces before and after the TMD is annotated as being.
                            #print("Since", record.id, "has a TRANSMEM of less than 25 residues, it will be taken to the next stage.")
                            TMD = f.extract(record.seq)

                #Here we define the ID code and then the start and end positions.
                            title_line = "'%s' | TMH: '%i-%i'" % (record.id, f.location.start+1, f.location.end)

                #This prints the title and then any TRANSMEM regions associated with that entry along with flanking regions.
                            name = str(record.id)
                            with open ("near_c_terminal_single_pass_list.txt", "ab") as single_fasta_file:
                                single_fasta_file.write(record.id)
                                single_fasta_file.write("\n")

                #if there are any TRANSMEM entries, it prints them. Note that if there is no transmem domain then the ID is skipped, you won't get flanking regions without an identified TMD.
                #Prints the output in the terminal

    output.close()
