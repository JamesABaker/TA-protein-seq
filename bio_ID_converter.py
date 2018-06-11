import urllib
import urllib2
import os
import os.path
import re
import fileinput
import sys
import urllib
import urllib2
import subprocess
import shutil
import re


# This opens a file (list.txt), replaces the formatting for the uniprot
# url handler
f = open('list.txt', 'r')
with open('list.txt') as file1list:
    fetch_TMs = file1list.readlines()
# In order to handle each entry of the list text file it needs to be
# converted to a pythonic list.
fetch_TMs = list(fetch_TMs)
# removes new line artefacts of string.
fetch_TMs = [fetch_TMs.replace('\n', '') for fetch_TMs in fetch_TMs]

prot_list = str(fetch_TMs)
prot_list = prot_list.replace("[", "")
prot_list = prot_list.replace("]", "")
prot_list = prot_list.replace("'", "")
prot_list = prot_list.replace(",", "")
prot_list = prot_list.replace("\n", "")

# This is the uniprot URL mapping part.

url = 'http://www.uniprot.org/mapping/'

params = {
    # This is the format the list is in. see
    # http://www.uniprot.org/help/uploadlists for help with this and the
    # abbreviations.
    'from': 'P_REFSEQ_AC',
    'to': 'ACC',  # This is what the list will be converted too.
    'format': 'list',
    'query': prot_list,
}

data = urllib.urlencode(params)
request = urllib2.Request(url, data)
# Please set your email address here to help us debug in case of problems.
contact = ""
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib2.urlopen(request)
page = response.read(99999999999999)


printer = {
    # This is the format the list is in. see
    # http://www.uniprot.org/help/uploadlists for help with this and the
    # abbreviations.
    'from': 'P_REFSEQ_AC',
    'to': 'ACC',  # This is what the list will be converted too.
    'format': 'tab',
    'query': prot_list,
}

data = urllib.urlencode(printer)
request = urllib2.Request(url, data)
# Please set your email address here to help us debug in case of problems.
contact = ""
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib2.urlopen(request)
page_for_reading = response.read(9999999999999)


print page_for_reading

text_file = open("Output.txt", "w")
text_file.write("%s" % page)
text_file.close()

'''
foogoo = list(page)

foogoo = [s for s in foogoo if not re.search(r'\NP',s)]

print foogoo
'''
