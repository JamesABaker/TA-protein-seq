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

        # This is set to False after the first features are counted so the
        # final TMD count isn't calculated twice.
        new_record = True
        # Count through each feature
        tmd_count = 0

        #This sees if any annotation claims the multi or single topology. There are several examples where this is wrong. We still need to check the actual number of TMDs.
        if "; Single-pass" in str(record):
            single_or_multi_topology = "Single-pass"
        if "; Multi-pass" in str(record):
            single_or_multi_topology = "Multi-pass"

        #Finds the total number of TMDs
        for i, f in enumerate(record.features):
            #feature type refers to TRANSMEM
            if f.type == feature_type:

                if "UnknownPosition" in str(f.location):
                    print(record.id, "contained a TMD without sequence number information.")

                else:
                    max_tmd_count = tmd_count + 1

        tmd_count = 0

        #Recounting the number of TMDs to see if the annotation for multi or single.
        # Checking if single-pass.
        if tmd_count > 1:
            single_or_multi_topology = "Multi-pass"
        if tmd_count == 1:
            single_or_multi_topology = "Single-pass"

        if single_or_multi_topology == "Single-pass":

            for i, f in enumerate(record.features):
                #feature type refers to TRANSMEM
                if f.type == feature_type:

                    if "UnknownPosition" in str(f.location):
                        print(record.id, "contained a TMD without sequence number information.")

                    else:
                        tmd_count = tmd_count + 1

                        C_terminal_distance = len(record.seq) - f.location.end+1 #python counts from 0


                        if C_terminal_distance < 15 :
                #These are the flanking regions and the TMD using record.seq. They are identified by simply counting 5 spaces before and after the TMD is annotated as being.
                            #print("Since", record.id, "has a TRANSMEM of less than 25 residues, it will be taken to the next stage.")
                            TMD = f.extract(record.seq)


                        #Here we define the ID code and then the start and end positions.
                        title_line = "'%s' | TMH: '%i-%i'" % (record.id, f.location.start+1, f.location.end)

                        #This prints the title and then any TRANSMEM regions associated with that entry along with flanking regions.
                        output.write(">%s\n%s\n" % (title_line, TMD))






output.close()
