# Genome Mining Script
Script for genome mining in search of lasso peptides

## Installation instructions
#### On a Linux/MacOS system:
Open a terminal and type the following commands in a terminal:
```
git clone https://github.com/nathanalam/genomemining.git
pip3 install virtualenv
python3 -m venv genomemining
cd genomemining
```
The command line should now have a (genomemining) tag behind it.

Then, run the following script to install the required libraries:
```
source shellscripts/setup.sh
```

### To Generate the Lasso Peptides from a set of genomes
Copy and paste all of the FASTA files (.fna or .faa) of the genomes into a folder called 'genomes' (CASE-SENSITIVE) in your working directory.
Then, open up a terminal and type in the following:
```
python Translator.py
python Miner.py
```

This script will then read the genomes, convert any .fna files into .faa files by translating the DNA, and then store the results in a JSON file "matches.json" and a CSV file "matches.csv", under the output folder.

### Viewing results in a browser
To make the results visible from a website, there is code to set up a small server that allows access to the lasso peptides via RESTful requests. To start the server, open up a terminal and type:
```
source shellscripts/run.sh
```

Then, navigate to the static directory and open up the "index.html" folder there on your local machine. This well let you search through your array of lasso peptides.
