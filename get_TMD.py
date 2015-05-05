
#Generates overall TMH_and_flanking


#This requires a working version of Biopython.
from Bio import SeqIO





#These are the variables that are repeatedly used throughout the script. The only one that changes is the output_filename that is changed as the separate fasta sequences are generated.
filenames = ["input.dat"]
input_format = "swiss" #This SHOULD work with uniprot filetype. From the seqIO biopython wiki: Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also the UniProt XML format.
feature_type = "TRANSMEM" #For future modification, this can be used to look for any annotation in the .dat file.
output_filename = "TMD.fasta" #Simply the output name, can be anything as it is written in binary (not file-type specific language).

output = open(output_filename, "w")
for filename in filenames:
    # Using SeqIO.parse will cope with multi-record files
    for record in SeqIO.parse(filename, input_format): #A biopython module that can automatically parse uniprot .txt files
        for f in record.features:
            #feature type refers to TRANSMEM
            if f.type == feature_type:
                #These are the flanking regions and the TMD using record.seq. They are identified by simply counting 5 spaces before and after the TMD is annotated as being.

                TMD = f.extract(record.seq)

                #Here we define the ID code and then the start and end positions.
                title_line = "'%s' | TMH: '%i-%i'" % (record.id, f.location.start+1, f.location.end)

                #This prints the title and then any TRANSMEM regions associated with that entry along with flanking regions.
                output.write(">%s\n%s\n" % (title_line, TMD))

                #if there are any TRANSMEM entries, it prints them. Note that if there is no transmem domain then the ID is skipped, you won't get flanking regions without an identified TMD.
                #Prints the output in the terminal
                '''
                print "ID: ", title_line

                print "Trans Membrane Domain: ", f.extract(record.seq)

                print "fasta output to file: \n>%s\n%s\n" % (title_line, TMD)
                '''




output.close()
