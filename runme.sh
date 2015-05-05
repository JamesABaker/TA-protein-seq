#!/bin/bash

#
# -- SGE options :
#

#$ -S /bin/bash
#$ -cwd
#$ -q

#
# -- the commands to be executed (programs to be run) :
#
# remember to run `qrsh -l inter -l short` before running on the cluster!!!!

DATE=$(date +"date:%Y.%m.%d_time:%H:%M:%S")

echo
echo "This script was developed by James A Baker under the supervision of Dr Jim Warwicker."
echo
echo "It requires an active internet connection and an up to date version of biopython."
echo
echo "See readme.md for more information on installation and visit www.github.com/jbkr/TA_predict to report any errors."
sleep 2

echo
echo "Extracting TMDs from .dat file marked by TRANSMEM. This includes confirmed TMDs and predicted TMDs according to a consensus of TMHMM, Memsat, Phobius and the hydrophobic moment plot method of Eisenberg and coworkers."
echo

python get_TMD.py

echo
echo "Getting a list of the single pass proteins by filtering out proteins with more than one TRANSMEM annotated TMD associated with their ID."
echo

python single_pass.py

echo
echo "Filtering each entry according to C terminal distance with a cutoff of the TRANSMEM annotated TMD being  less than 25 residues from the C-terminal."
echo

python check_C.py

echo
echo "Analysis complete. A text file containing all the single pass transmembrane proteins with their TMD near the C terminal are contained in a text file."
echo


mkdir ./date_$DATE
mv TMD.fasta ./date_$DATE/TRANSMEM.fasta
mv single_pass_list.txt ./date_$DATE/single_TRANSMEM.txt
mv near_c_terminal_single_pass_list.txt ./date_$DATE/C_terminal_single_TRANSMEM.txt


echo
echo "Finished moving files."
