# Tail Anchor Prediction

## Introduction

This bundle of scripts aims to filter a uniprot text file into a list of potential Tail Anchors. Basically the list filters a uniprot text file to a list of uniprot IDs that are likely tail anchors.

### Academic abstract originally accepted at FEBS Congress 2015.
#### A bioinformatic method to identify potential SNARE proteins.

Tail anchored proteins are a topologically distinct class of intracellular proteins defined by their single carboxy-terminal transmembrane domain with a cytosolic facing amino-terminus. Tail anchored proteins are involved in a range of key cellular functions including protein translocation and apoptosis. Additionally, within the tail anchored class of proteins are a set of vesicle fusion proteins called SNARE proteins. There is biomedical interest in SNARE drug delivery mechanisms. SNAREs can fuse liposomes containing various drug payloads into the membrane. This study aims to identify SNARE proteins in eukaryotic proteomes by filtering through large datasets using automatically predicted TrEMBL consensus, and manually annotated SWISS-PROT transmembrane regions. The pipeline generates a list of singlepass proteins with a transmembrane domain close to the C terminal, that are not splice isoforms. A previous study by Kalbfleisch _et al._ published in Traffic 2007 (**8**: _1687-1694_) predicted 411 tail anchor proteins. This study uses more stringent filtering methods, and a larger dataset, to identify 351 novel predicted tail anchored proteins from a comprehensive human dataset. The tools developed herein are openly available for re-application to other datasets. Notably, known SNARE transmembrane helices are highly hydrophobic even compared to other tail anchored transmembrane helices. We compare Kyte and Doolittle hydrophobicity profiles of our filtered human protein list against the profiles of previously known SNARE and tail anchored proteins. This provided a list of potential SNARE proteins in addition to potential spontaneously inserting tail anchored proteins similar to cytochrome b5 which have the least hydrophobic transmembrane helices.

## Running the Prediction

Save your uniprot text downloaded dataset in the folder containing the scripts.

Then:

 - Open a terminal
 - Navigate to the folder containing the scripts. (*`cd Downloads/TAPredict` for example*.)
 - Type `bash runme.sh` and hit enter.
 - The script will prompt you to enter the name of your input file.
 - The results will be saved in a folder with the current time and date.

 It is important that you do not add files to the folder, or tamper with input files, or script files whilst the script is running as the script may not report an error and blunder on recklessly with erroneous I/O files.


## Filter Method

There are a series of filters to deduce tail anchored proteins.

<img src="TA_filter_flow.png" width="400">

1. The script filters proteins by those with `TRANSMEM` regions. `TRANSMEM` annotation includes experimentally confirmed TMDs and predicted TMDs. Predictions of TMDs are according to a consensus of TMHMM, Memsat, Phobius and the hydrophobic moment plot method of Eisenberg and coworkers and is calculated by uniprot itself.
2. The TMDs of each protein are counted. If a protein has more than 1 `TRANSMEM` region, it is not added to the list. The list now contains those protein IDs with only a single `TRANSMEM` region.
3. The script counts the distance between the `TRANSMEM` region and the C terminus. If the final residue of the `TRANSMEM` annotated region is within 15 residues of the C terminal residue, the ID is added to the list.
4. Finally the script removes any proteins that contain `NON_TER` annotation. This removes potential signal anchor protein and multipass protein splice isoforms from the list.

This final list, and any intermediate list, can be directly uploaded for batch retrieval from [***uniprot***](http://www.uniprot.org/uploadlists) for more information.

## System Requirements

This script requires python 2.7, biopython, and an active internet connection.

#### Installing Biopython:

 In **OSX** or **Linux**:

 - Open a terminal.
 - Run the following commands, entering your password where necessary:

 	`sudo easy_install pip; sudo pip install numpy; sudo pip install Biopython`

## Reporting errors

If you run into any problems or have any suggestions and corrections, please, [log an issue on
GitHub](https://github.com/jbkr/TApredict/issues/new).

**Common problems.**

If you come across any errors during the installation it is probably because python is not installed in the default locations (and you probably know how to fix it), or the package has already been installed before you did these commands. Type in the above commands one after the other regardless and then try running the script.

If errors come up during running the script, it is probably because network connection has dropped, or files have been moved or added prematurely.
