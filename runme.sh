#!/bin/bash

#
# -- SGE options :
#

#$ -S /bin/bash
#$ -cwd
#$ -q

################################################
##########       Version 0.1       #############
##########       James Baker       #############
################################################

#
# -- the commands to be executed (programs to be run) :
#
# remember to run `qrsh -l inter -l short` before running on the cluster!!!! Running on the cluster is not advised since this is a somewhat interactive sript and pulls information from different modules and applications.

DATE=$(date +"date:%Y.%m.%d_time:%H:%M:%S")

echo
echo "This script was developed by James A Baker under the supervision of Dr Jim Warwicker."
echo
echo "It requires an active internet connection and an up to date version of biopython."
echo
echo "See readme.md for more information on installation and visit www.github.com/jbkr/TA_predict to report any errors."
sleep 2

echo
echo
echo
echo "What is the name of your input file? (include extension i.e file.dat or file.txt). Do not run this script using a file named 'input.dat'"
echo
echo

read filename

cp $filename input.dat
echo

echo "Extracting TMDs from .dat file marked by TRANSMEM. This includes confirmed TMDs and predicted TMDs according to a consensus of TMHMM, Memsat, Phobius and the hydrophobic moment plot method of Eisenberg and coworkers."
echo

python get_TMD.py

echo
echo "Getting a list of the single pass proteins by filtering out proteins with more than one TRANSMEM annotated TMD associated with their ID."
echo

python single_pass.py

echo
echo "Filtering each entry according to C terminal distance with a cutoff of the TRANSMEM annotated TMD being  less than 15 residues from the C-terminal."
echo

python check_C.py

echo
echo "Checking for any non terminal annotations to remove the possibility of splice varients of signal anchors or multipass proteins."
echo

python check_splice.py

echo
echo "Analysis complete. A text file containing all the single pass transmembrane proteins with their TMD near the C terminal are contained in a text file."
echo


mkdir ./$filename_date_$DATE
mv TMD.fasta ./$filename_date_$DATE/TRANSMEM.fasta
mv single_pass_list.txt ./$filename_date_$DATE/single_TRANSMEM.txt
mv near_c_terminal_single_pass_list.txt ./$filename_date_$DATE/C_terminal_single_TRANSMEM.txt
mv splice_variant.txt ./$filename_date_$DATE/splice_variant.txt

rm uniget.dat
rm transmembranes.fasta
rm input.dat

echo
echo "Finished moving files."
