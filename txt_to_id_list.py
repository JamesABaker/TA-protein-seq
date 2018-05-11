from Bio import SeqIO
import sys

input_file = str(sys.argv[1])
input_format = "swiss"

for record in SeqIO.parse(input_file, input_format):
    print record.id
