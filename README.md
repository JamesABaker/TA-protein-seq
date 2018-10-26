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

The file structure is outlined below:

        ../TA-protein-seq
        ├── Interaction
        │   ├── get3-A+B-interactors-uniprot-11-6-2018.txt
        │   ├── get3-BIOGRID-GENE-31962-3.4.161.tab2-11-6-2018.txt
        │   ├── get3-interactors-A+B-Entrez.txt
        │   ├── hsc70-A+B-interactors-uniprot-11-6-2018.txt
        │   ├── hsc70-BIOGRID-GENE-109544-3.4.161.tab2-11-6-2018.txt
        │   ├── hsc70-interactor-A+B-Entrez.txt
        │   ├── hsp-interactors-A+B-Entrez.txt
        │   ├── hsp40-A+B-interactors-uniprot-11-6-2018.txt
        │   ├── hsp40-BIOGRID-GENE-119699-3.4.161.tab2-11-6-2018.txt
        │   ├── interaction\ command.txt
        │   ├── interaction-chaperone-profiles-log.txt
        │   ├── sgt2-A+B-interactors-uniprot-11-6-2018.txt
        │   ├── sgt2-BIOGRID-GENE-34410-3.4.161.tab2-11-6-2018.txt
        │   ├── sgt2-interactors-A+B-Entrez.txt
        │   ├── sgta-A+B-interactors-uniprot-11-6-2018.txt
        │   ├── sgta-BIOGRID-GENE-112347-3.4.161.tab2-11-6-2018.txt
        │   ├── sgta-interactors-A+B-Entrez.txt
        │   ├── sl-9908-7-6-2018.csv
        │   ├── snd1-A+B-interactors-uniprot-11-6-2018.txt
        │   ├── snd1-BIOGRID-GENE-32240-3.4.161.tab2-11-6-2018.txt
        │   ├── snd1-interactors-A+B-Entrez.txt
        │   ├── trc40-A+B-interactors-uniprot-11-6-2018.txt
        │   ├── trc40-BIOGRID-GENE-106931-3.4.161.tab2-11-6-2018.txt
        │   └── trc40-interactors-A+B-Entrez.txt
        ├── SL-9908_24_4_2018
        │   ├── Unfiltered
        │   │   ├── SL-9908_24-4-2018.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_er.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_golgi.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_human-and-mouse.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_human.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_mito.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_mouse.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_plant.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_pm.csv
        │   │   ├── SL-9908_24_4_2018-after_cdhit_yeast.csv
        │   │   ├── SL-9908_24_4_2018_afterCDHIT.csv
        │   │   └── cd-hit
        │   │       ├── LOG
        │   │       ├── README
        │   │       ├── output
        │   │       ├── output-sorted.clstr
        │   │       ├── output.clstr
        │   │       └── output.log
        │   ├── cd-hit
        │   │   └── sl-9908-7-6-2018-cd-hit.txt.clstr
        │   ├── sl-9908-7-6-2018-cd-hit-Mouse.csv
        │   ├── sl-9908-7-6-2018-cd-hit-er.csv
        │   ├── sl-9908-7-6-2018-cd-hit-golgi.csv
        │   ├── sl-9908-7-6-2018-cd-hit-human.csv
        │   ├── sl-9908-7-6-2018-cd-hit-mammal.csv
        │   ├── sl-9908-7-6-2018-cd-hit-mito.csv
        │   ├── sl-9908-7-6-2018-cd-hit-plant.csv
        │   ├── sl-9908-7-6-2018-cd-hit-pm.csv
        │   ├── sl-9908-7-6-2018-cd-hit-yeast.csv
        │   └── sl-9908-7-6-2018-cd-hit.csv
        ├── TA_filter_flow.pdf
        ├── TA_filter_flow.png
        ├── TrEMBL_25_4_2018
        │   ├── TrEMBL-TRANSMEM_25-4-2018.csv
        │   ├── TrEMBL-TRANSMEM_25-4-2018_aftercdhit.csv
        │   └── cd-hit
        │       ├── LOG
        │       ├── README
        │       ├── output
        │       ├── output-sorted.clstr
        │       ├── output.clstr
        │       └── output.log
        ├── baseline-hydrophobicity
        │   ├── UniER_5_flanklength_flankclashTrue.csv
        │   ├── UniGolgi_5_flanklength_flankclashTrue.csv
        │   ├── UniPM_5_flanklength_flankclashTrue.csv
        │   ├── mean-results.txt
        │   └── nni-csv-hydrophobicity.py
        ├── filter_by_list_results.txt
        ├── heatmaps
        │   ├── Figure_1_2_3_4_6(nni-clone).py
        │   ├── heatmaps\ only.xlsx
        │   └── heatmaps.log
        ├── readme.md
        ├── reference_list.csv
        ├── scripts
        │   ├── bio_ID_converter.py
        │   ├── eisenberg.pl
        │   ├── filter-by-list.py
        │   ├── filter_with_topology.py
        │   ├── graphs.py
        │   ├── hessa.pl
        │   ├── hydrophobicity_profile.py
        │   ├── kd.pl
        │   ├── list_checker.py
        │   ├── no_topology_filter.py
        │   ├── stat_two_lists.py
        │   ├── tail_lengths.py
        │   └── txt_to_id_list.py
        ├── spontaneous_TA_case_study
        │   ├── Modeller_9375144(1).pdb
        │   ├── consurf2cytb5.png
        │   ├── consurf2cytb5180.png
        │   ├── consurf_new.py
        │   ├── cytb5-180.png
        │   ├── cytb5-clean.pdb
        │   ├── cytb5.pdb
        │   ├── cytb5.png
        │   ├── cytb5_clean_pdb_ATOMS_section_With_ConSurf.pdb
        │   ├── cytb5_clean_pdb_With_Conservation_Scores.pdb
        │   ├── cytb5electro-180.png
        │   ├── cytb5hydro-180.png
        │   ├── cytb5tmhelectro.png
        │   ├── cytb5tmhhydro.png
        │   ├── cytochrome\ b5\ figure.pdf
        │   ├── hydrophobicity.py
        │   ├── ptp1bhelix.png
        │   ├── ptp1bhelix180.png
        │   ├── ptp1bhelixelectro.png
        │   ├── ptp1bhelixelectro180.png
        │   ├── ptp1bhelixhydrophobicity.png
        │   ├── ptp1bhelixhydrophobicity180.png
        │   ├── ptphelixclean.pdb
        │   └── view\ for\ cyt\ b5\ pymol.py
        └── swissprot_24_4_2018
            ├── SwissProt_24-4-2018_AfterCDHIT.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_er.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_golgi.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_human.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_human_and_mouse.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_mito.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_mouse.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_plants.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_pm.csv
            ├── SwissProt_24-4-2018_AfterCDHIT_yeast.csv
            ├── SwissProt_24_4_2018.csv
            └── cd-hit
                ├── LOG
                ├── README
                ├── output
                ├── output-sorted.clstr
                ├── output.clstr
                └── output.log



## System Requirements

This script requires python and biopython.

#### Installing Biopython:

 In **OSX** or **Linux**:

 - Open a terminal.
 - If you've not used pip before, run the following commands, entering your password where necessary:

 	`sudo easy_install pip; pip install numpy; pip install Biopython`

## Reporting errors

If you run into any problems or have any suggestions and corrections, please, [log an issue on
GitHub](https://github.com/jbkr/TApredict/issues/new).

**Common problems.**

If you come across any errors during the installation it is probably because python is not installed in the default locations (and you probably know how to fix it), or the package has already been installed before you did these commands. Type in the above commands one after the other regardless and then try running the script.

If errors come up during running the script, it is probably because network connection has dropped, or files have been moved or added prematurely.
