from Bio import SeqIO


filenames = ["input.dat"]
input_format = "swiss" #This works with uniprot filetype. From the seqIO biopython wiki: Swiss-Prot aka UniProt format. Uses Bio.SwissProt internally. See also the UniProt XML format if something goes wrong.
feature_type = "TRANSMEM" #For future modification, this can be used to look for any annotation in the file.
output_filename_fasta = "TMD.fasta" #Simply the output name, can be anything as it is written in binary (not file-type specific language).
other_feature_type = "NON-TER"
signal_feature = "SIGNAL"


output = open(output_filename_fasta, "w")
for filename in filenames:
    multipass_count = 0
    singlepass_count = 0
    within25 = 0
    nosignal = 0
    print "Parsing records. This may take some time..."
    for record in SeqIO.parse(filename, input_format): #A biopython module that can automatically parse uniprot .txt files


        # This is set to False after the first features are counted so the
        # final TMD count isn't calculated twice.
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




        #Recounting the number of TMDs to see if the annotation for multi or single.
        # Checking if single-pass.
        if tmd_count > 1:
            single_or_multi_topology = "Multi-pass"
        if tmd_count == 1:
            single_or_multi_topology = "Single-pass"

        if tmd_count == "Multi-pass":
            multipass_count = multipass_count + 1

        if tmd_count == "Single-pass":
            singlepass_count = singlepass_count + 1



        tmd_count = 0
        nonterminalrecordids = 0
        unknownexceptions = 0

        if single_or_multi_topology == "Single-pass":

            for i, f in enumerate(record.features):




                #The non_ter annotation is used in the case of splice isoforms.
                if f.type == other_feature_type:
                    print "Non terminal annotation found:", record.id
                    nonterminalrecordids = nonterminalrecordids + 1


                elif f.type == feature_type:
                    if "UnknownPosition" in str(f.location):
                        print "contained a TMD without sequence number information:", record.id

                    else:
                        C_terminal_distance = len(record.seq) - f.location.end+1 #python counts from 0

                        #Checking distance from the C terminal
                        if C_terminal_distance < 25 :
                            #These are the flanking regions and the TMD using record.seq. They are identified by simply counting 5 spaces before and after the TMD is annotated as being.
                            #print("Since", record.id, "has a TRANSMEM of less than 25 residues, it will be taken to the next stage.")
                            TMD = f.extract(record.seq)
                            full_sequence = str(record.seq)
                            tmh_start = int(f.location.start)
                            tmh_stop = int(f.location.end)
                            tmh_sequence = str(record.seq[(f.location.start):(f.location.end)])
                            within25=within25 + 1

                            #Checking for Signal peptide annotation
                            signal_present = False
                            for x in (record.features):
                                if x.type == signal_feature:
                                    signal_present = True
                            if signal_present == True:
                                print "Signal detected in:", record.id

                            if signal_present == False:
                                nosignal = nosignal + 1

                                #checking for cytoplasmic facing regions
                                if "Cyto" in
                                    #Printing fasta for hits
                                    print "Potential tail anchored protein:", record.id
                                    #Here we define the ID code and then the start and end positions.
                                    title_line = "'%s' | TMH: '%i-%i'" % (record.id, f.location.start+1, f.location.end)
                                    #This prints the title and then any TRANSMEM regions associated with that entry along with flanking regions.
                                    output.write(">%s\n%s\n" % (title_line, TMD))


    print "Multipass records:", multipass_count
    print "Singlepass records:", singlepass_count
    print "Singlepass within 15 of C terminal records:", within15
    print "Records that also didn't contain a signal peptide:", nosignal

output.close()
