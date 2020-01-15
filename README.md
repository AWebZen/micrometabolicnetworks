# micrometabolicnetworks
Building and analysing prokaryotic metabolic networks in regard to their optimal growth temperature. 

Code for paper "Robust structure measures of metabolic networks that predict prokaryotic optimal growth temperature", from Weber Zendrera, Sokolovska and Soula (2019, BMC Bioinformatics)

For details on the algorithm and results, see the paper.

As of 2020, added code for metabolic network scope analysis, for our second paper on the subject.

## Dependencies
Python 2.7
The package requirements can be found in **requirements.txt**

## Installation
Download repertory.
Unzip **database\_species\_list\_assembly\_files.zip**. It is then ready to use.

-----------------
## To replicate our paper results from 2019, run :

### To build our metabolic networks
- graphs\_takemoto\_species.py
- graphs\_Bacdivedb\_species.py
- graphs\_HPMCdb\_species.py

### To analyse our metabolic networks
- analysis\_graph\_properties.py
- analysis\_laplacian\_spectrum.py


## To replicate our paper results from 2020, after building the networks, run:
- scope_analysis.py

-----------------
An example of use of our toolbox **utils.py** can be found in **example.py**


## File overview

**requirements.txt** - Package requirements

**12859\_2006\_1675\_MOESM1\_ESM.xls** - Takemoto et al 2007 supplementary figure (Takemoto, K., Nacher, J.C., and Akutsu, T. (2007). Correlation between structure and temperature in prokaryotic metabolic networks. *BMC Bioinformatics* 8, 303.)
In **database\_species\_list\_assembly\_files.zip** :
  **assembly\_summary\_prokaryotes.txt** - Genbank's assembly summary (archaea, bacteria), March 2018
  **bacteriaClassification\_hpmcd\_microbiota.txt** - HPMC database bacteria classification, May 2018
  **species.txt** - Ensembl species summary file, March 2018

**utils.py** and **utils\_general.py** - toolbox with functions to build and analyse metabolic networks

**graphs\_takemoto\_species.py** - to build networks for the species from Takemoto et al 2007 paper
**graphs\_Bacdivedb\_species.py** - to find and build networks for BacDive species
**graphs\_HPMCdb\_species.py** - to find species and build networks from HPMC database

**analysis\_graph\_properties.py** - analysing and plotting different network structural properties vs optimal growth temperature
**analysis\_laplacian\_spectrum.py** - evaluating and plotting Laplacian spectrum
**scope\_analysis.py** - evaluating and analysing/plotting metabolic network scopes, and using a spectral clustering on the species

**example.py** - toy example of usage of our toolbox

**get\_cDNA1.0.sh** - finding cDNA fasta file through species name
**get\_cDNA2.0.sh** - finding cDNA fasta file through taxonomy ID
**get\_cDNA3.0.sh** - finding cDNA fasta file through assembly ID

