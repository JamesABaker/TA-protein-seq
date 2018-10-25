# Tail Anchor Prediction

## Introduction

This bundle of scripts aims to filter a uniprot text file into a list of potential Tail Anchors. Basically the list filters a uniprot text file to a list of uniprot IDs that are likely tail anchors.

### Academic abstract originally accepted at FEBS Congress 2015.

Initially the project aimed to identify highly hydrophobic SNARE TMHs from a list of TA proteins.

>#### A bioinformatic method to identify potential SNARE proteins.
>
>Tail anchored proteins are a topologically distinct class of intracellular proteins defined by their single carboxy-terminal transmembrane domain with a cytosolic facing amino-terminus. Tail anchored proteins are involved in a range of key cellular functions including protein translocation and apoptosis. Additionally, within the tail anchored class of proteins are a set of vesicle fusion proteins called SNARE proteins. There is biomedical interest in SNARE drug delivery mechanisms. SNAREs can fuse liposomes containing various drug payloads into the membrane. This study aims to identify SNARE proteins in eukaryotic proteomes by filtering through large datasets using automatically predicted TrEMBL consensus, and manually annotated SWISS-PROT transmembrane regions. The pipeline generates a list of singlepass proteins with a transmembrane domain close to the C terminal, that are not splice isoforms. A previous study by Kalbfleisch _et al._ published in Traffic 2007 (**8**: _1687-1694_) predicted 411 tail anchor proteins. This study uses more stringent filtering methods, and a larger dataset, to identify 351 novel predicted tail anchored proteins from a comprehensive human dataset. The tools developed herein are openly available for re-application to other datasets. Notably, known SNARE transmembrane helices are highly hydrophobic even compared to other tail anchored transmembrane helices. We compare Kyte and Doolittle hydrophobicity profiles of our filtered human protein list against the profiles of previously known SNARE and tail anchored proteins. This provided a list of potential SNARE proteins in addition to potential spontaneously inserting tail anchored proteins similar to cytochrome b5 which have the least hydrophobic transmembrane helices.

Since then the project has been developed to investigate subtle, but biologically distinctive, compositional differences between various groups of TA TMPs. The aim of these scripts is to generate a list of TA proteins from various sources and analyse sequence distributions of the residues.
There are a series of filters to deduce tail anchored proteins.
Then, the lists can be stratified into distinct biological groups, for example different species or different sub cellular organelles.

<img src="TA_filter_flow.png" width="800">


## Running the Prediction

Save your uniprot text downloaded dataset in the folder containing the scripts.

Then:

 - Open a terminal
 - Navigate to the folder containing the scripts. (*`cd Downloads/TAPredict` for example*.)
 - Run `python filters.py YOURFILE.txt`

 `YOURFILE` can be .txt or .dat downloaded from Uniprot as a text file.

## Filter Method



## System Requirements

This script requires python and biopython.

#### Installing Biopython:

 In **OSX** or **Linux**:

 - Open a terminal.
 - If you've not used pip before, run the following commands, entering your password where necessary:

 	`sudo easy_install pip; sudo pip install numpy; sudo pip install Biopython`

## Reporting errors

If you run into any problems or have any suggestions and corrections, please, [log an issue on
GitHub](https://github.com/jbkr/TApredict/issues/new).

**Common problems.**

If you come across any errors during the installation it is probably because python is not installed in the default locations (and you probably know how to fix it), or the package has already been installed before you did these commands. Type in the above commands one after the other regardless and then try running the script.

If errors come up during running the script, it is probably because network connection has dropped, or files have been moved or added prematurely.
