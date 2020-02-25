# Genome Mining Script
Script for genome mining in search of lasso peptides

## Installation instructions

### Pre-requisites
This software relies on the MEME suite for identifying motifs using MAST, and it relies on NCBI's Entrez system for downloading genomes based on accession numbers.

#### MEME Suite:
To download MEME, follow these instructions: http://meme-suite.org/doc/install.html?man_type=web#quick

#### E-Direct:
To download Entrez, type in the following into a terminal:
```
sudo apt-get install entrez-direct
```
#### GitHub Codebase:
Open a terminal and type the following commands in a terminal:
```
git clone https://github.com/nathanalam/genomemining.git
pip3 install virtualenv
python3 -m venv genomemining
cd genomemining
```
The command line should now have a (genomemining) tag behind it.

#### Python libraries:
Then, run the following script to install the required libraries:
```
source shellscripts/setup.sh
```

## Run the script

The following is an example command:
```
esearch -db nucleotide -query "MN695290.1" | efetch -format fasta | python3 lassoMine.py
```