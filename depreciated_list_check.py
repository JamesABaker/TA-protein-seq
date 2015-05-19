###############################################################
########################## BioLC ##############################
##################### BIOLIST-CHECKER #########################
########################## V 1.1 ##############################
###############################################################

# Authors: James Baker

# UniPROT query notes:
# keyword:transmembrane

# Getting all the system stuff

import sys
import os
import os.path
import re
import fileinput
import time
import sys

print "\n", "Imported system information. \n"

#First file open
task1 = raw_input('What is the first database list called? \n If your database is called "example.txt" then enter "example".')
print 'Opening',task1, '.txt... \n'
f = open('%s.txt' % task1,'r')
print "The database text file contains the following items:\n", (f.read(200)), "... \n", "et cetera... \n"
with open('%s.txt' % task1,) as file1list:
    interaction = file1list.readlines()

#opens second file to compare
task2 = raw_input('What is the second database list called? \n If your database is called "example.txt" then enter "example".')
print 'Opening',task2, '.txt... \n'
f = open('%s.txt' % task2,'r')
print "The database text file contains the following items:\n", (f.read(200)), "... \n", "et cetera... \n"
with open('%s.txt' % task2,) as file2list:
    tms = file2list.readlines()

print "Opened database files successfully. \n"

#Compares the intersections of the lists i.e finds the matches
print "setting sets"
output = set(interaction) & set(tms)
output = str(output)
for char in "'][,\)(nset":
    output = output.replace(char,'')

#reports hits
print "The following proteins are hits: \n \n \n", output, "\n \n"

#writes hits to file
print 'Printing hits to output file! ...in progress... \n'
outputname = raw_input('What do you want to call this new list? \n Include extensions. \n Example: "Transmembrane_MOM_human.txt"')
file = open("output.txt", "w")
file.write(str(output))
file.close()
print "Saving output as", outputname
os.rename("output.txt", outputname)
print 'BANG! \n And the printing is done. The output file has been made! \n \n If your output file is empty the script ran and found no matches. \n \n \n Copy the printed output or upload the ouput.txt to www.uniprot.org/uploadlists for more information on these proteins.'

'''
import urllib,urllib2

url = 'http://www.uniprot.org/mapping/'

params = {
'from':'ACC',
'to':'P_REFSEQ_AC',
'format':'tab',
'query':'transmembrane'
}

data = urllib.urlencode(params)
request = urllib2.Request(url, data)
contact = "james.baker-5@postgrad.manchester.ac.uk" # Please set your email address here to help us debug in case of problems.
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib2.urlopen(request)
page = response.read(20000000)

print page
'''
