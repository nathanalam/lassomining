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
python3 -m venv lassomining
cd genomemining
```
The command line should now have a (lassomining) tag behind it.

#### Python libraries:
Then, run the following script to install the required libraries:
```
source shellscripts/setup.sh
```

#### Web Server Setup:
Finally, you must update the URL for the site to correspond to the URL that you have access to (i.e. one that can receive and respond to HTTP requests and run this process). To do this, change the instances of **50.116.48.39:8080/** in [Server.py](https://github.com/nathanalam/lassomining/blob/master/Server.py) and [shellscripts/view.sh](https://github.com/nathanalam/lassomining/blob/master/shellscripts/view.sh). 

To start the server, open up a terminal and type:
```
source shellscripts/view.sh
```

Then, go to the URL that you set BASE_URL to - it should be running web server. Otherwise, you may need to work on having the correct server administration credentials.
