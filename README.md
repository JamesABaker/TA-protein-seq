####Current Issues
These are issues that have not yet been logged on github, but need addressing.
 1. Paralogues are included and this skews the count. These should be removed at every stage.
 2. It is unclear if we have signal splice isoforms in our results.
 3. The results are comprehensive and irrespective of target membrane.
 4. Are *only* swissprot entries held at the end?

#Tail Anchor Prediction

##Introduction

This bundle of scripts aims to filter TRANSMEM annotated sequences from a uniprot text file into a list of what are potentially Tail Anchors.

##Filter Method

There are a series of filters to deduce tail anchored proteins.

- The script filters proteins by those with TRANSMEM regions. 
- Then by those with only a single TRANSMEM region.
- The script then gets an up to date copy of the file of the protein ID directly from uniprot and filters those by any that contain TRANSMEM regions within 25 residues of the C terminal residue.

This list can be directly uploaded for batch retrieval from [***uniprot***](http://www.uniprot.org/uploadlists).

##System Requirements

This script reguires python 2.7, biopython, and an active internet connection.

###Installing Biopython:
 
 In **OSX** or **Linux**:
 
 - Open a terminal.
 - Run the following commands:
  
 	`sudo easy_install pip`
 	
 	`pip install numpy`
 	
	`pip install Biopython`
	
*If you come across any errors it is probably because python is not installed in the default locations, or a package has already been installed before you did these commands. Type in the above commands one after the other regardless and then try running the script.*

##Running the Prediction

Save your uniprot text downloaded dataset as `input.dat` in the folder containing the scripts.

Then:

 - Open a terminal
 - Navigate to the folder. *`cd Downloads/Tail_Anchor_Prediction` for example*.
 - From within the directory containing the scripts run `bash runme.sh`.
 - The script will ask you for a folder name where your results will be saved.
 
*Please contact baker.james.jb@gmail.com if you run into any problems. If you need your issue addressing **quickly** and thoroughly, log an issue on GitHub.*