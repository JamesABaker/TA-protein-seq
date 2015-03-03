#####
##### The list of lists holds the [[title1, seq1], [title2, seq2]]
#####


with open ("TMD.fasta", "rb") as fastafile:
    data=fastafile.read().replace("'\n", "':::") #this ::: acts as a marker for makign the list of lists
    stringA = data
    #print stringA
    stringB = [] #Create an empty list
    for line in stringA.split('\n'): # Loop through each line in the original text.  Assumes these are newline line breaks which are represented by \n.
        splitLine = line.split(':::') #Split the pairs on the , into a 2 item list according to the marker
        stringB.append(splitLine) #Append that list to our empty list so thus it will end up with a list of lists.'''
seqs = stringB
#print seqs

#Now we have our list in a lit, we can replace the sequence with whatever values we want.

#However we want to create some more lists based on how many duplicates there are.

#Here we add a new entry to each list so it contains the header, the fasta, and a number of how many TMDs are in that protein.
for i, item in enumerate(seqs): #This is disgusting code.
    this_list = seqs[i]
    last_list = seqs[i - 1]
    #next_list = seqs[i + 1]

    this_header = this_list[0]
    last_header = last_list[0]
    #next_header = next_list[0]

    length_of_list = len(this_list)

    if length_of_list == 2:
        this_list.extend([1])
        if this_header[:10] == last_header[:10]: #if the IDs match (the first 10 characters), then we presume a multipass protein.
    #This isn't a great way of counting since it maxes out at 2. Perhaps looping the script 10 times.
            this_list[2] += 1
            TMDcount = this_list[2]
            last_list[2] = TMDcount
        else:
            pass


#Now that each fasta entry is tagged with a TMD count we can sort them into different lists.


for i in seqs:
    #This goes through the list and splits the items into different text files based on their TMD count.

    if i[2] == 1:
        ID_for_writing = i[0][:10]
        ID_for_writing = ID_for_writing.replace(">", "")
        ID_for_writing = ID_for_writing.replace("'", "")
        print ID_for_writing, " looks like a single pass protein.\nAdding to list.\n"
        with open ("single_pass_list.txt", "ab") as single_fasta_file:
            single_fasta_file.write(ID_for_writing)
            single_fasta_file.write("\n")

    elif i[2] == 2:
        print ID_for_writing, " looks like a multi pass protein.\nNot adding to list.\n"
    else:
        print "Something went wrong, no IDs could be counted."


########################


