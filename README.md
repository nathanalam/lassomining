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
git clone https://github.com/nathanalam/lassomining.git
pip3 install virtualenv
python3 -m venv lassomining
cd lassomining
```

Then, once inside the directory, activate the virtual environment:
```
source bin/activate
```
You will see a (lassomining) identifier in the command line now, indicating that you are in the virtual environment and are safe to freely download libraries without making global changes to your python installation. Now, you can install the necessary libraries (you can copy and paste the following into the command prompt):
```
pip3 install biopython
pip3 install pandas
pip3 install jupyter
```
## Operation
You can run this on a local computer, or you can run this using a [SLURM script](https://researchcomputing.princeton.edu/education/online-tutorials/getting-started/introducing-slurm)

### Configuration
The default directories are listed in *config.yaml* file within the repository. Change them as you see fit for your use case - the file is commented with more details about what each of the parameters mean.

This file holds parameters which will be specific to each run, and also other variables that you may want to alter.

### Example Usage
After editing the *config.yaml* file, the program can be run as follows:
```
python3 Mine.py
```

## Viewing Results
You can view the results using the jupyter notebook inside the directory. To run Jupyter notebook, enter the following into a command prompt (NOTE: you must be in the virtual environment, or in other words have the "(lassomining)" tag in the command line):
```
jupyter notebook
```
Then, you can enter the notebook and analyze the results as you wish. More instructions and example analyses are available on the notebook.
