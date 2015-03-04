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

####Installing Biopython:
 
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
 - Navigate to the folder containing the scripts. (*`cd Downloads/TAPredict` for example*.)
 - Run `bash runme.sh`.
 - The results will be saved in a folder with the current time and date.
 
*If you run into any problems or have any suggestions and corrections, please*, [log an issue on 
GitHub](https://github.com/jbkr/TApredict/issues/new).