# Genome Mining Script
Script for genome mining in search of lasso peptides

## Installation instructions

### Pre-requisites
This software relies on the MEME suite for identifying motifs using MAST, and it relies on NCBI's Entrez system for downloading genomes based on accession numbers.

#### MEME Suite:
To download MEME, follow these instructions: http://meme-suite.org/doc/install.html?man_type=web#quick

#### GitHub Codebase:
Open a terminal and type the following commands in a terminal:
```
git clone https://github.com/nathanalam/genomemining.git
pip3 install virtualenv
python3 -m venv lassomining
cd genomemining
```
The command line should now have a (lassomining) tag behind it.

#### Python libraries
To install the python libraries in the virtual environment, type the following:
```
git clone https://github.com/nathanalam/genomemining.git
pip3 install virtualenv
python3 -m venv lassomining
cd genomemining
```

#### File Directory Setup

The default directories are listed in Launch.py as the default arguments for the command line

All of these can be changed manually, but also during runtime


### Example Usage
Using the default directories, the program can be run as follows:
```
python3 Launch.py
```
Which can also be achieved, alternatively, with the shellscript that also ensures the python environment is set up:
```
source run.sh
```

The following lists all of the example arguments (all are optional, and defaults will be used instead)"
```
python3 Launch.py -name testRun -pattern M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\* -cutuff -100 -genome /tigress/nalam/genomeMining/genomes/* -o /tigress/nalam/genomeMining/matches.db -rundata /tigress/nalam/genomeMining/runs/ -model /tigress/nalam/genomeMining/models/b.faa:3:25,/tigress/nalam/genomeMining/models/c.faa:4:25
```

Note that the flag "-model" has the format such that individual models are separated by ",", and each model requires three sections, separated by ":", corresonding to the file location, the number of motifs to generate, and the maximum width per motif.