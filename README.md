# Genome Mining Script
Script for genome mining in search of lasso peptides

## Installation instructions

### Pre-requisites
This software relies on the MEME suite for identifying motifs using MAST, and it relies on NCBI's Entrez system for downloading genomes based on accession numbers.

#### MEME Suite:
To download MEME, follow these instructions: http://meme-suite.org/doc/install.html?man_type=web#quick

The following is *roughly* what you will do to install, but it's honestly a little convoluted and mileage may vary:
```
wget http://meme-suite.org/meme-software/5.1.1/meme-5.1.1.tar.gz
tar zxf meme-5.1.1.tar.gz
cd meme-5.1.1
./configure --prefix=$HOME/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
```
That's what they say on the website. Most of the time, things don't work though, because there is some pre-requisite software that isn't installed. Here are things that I found to be useful steps to resolve this (the goal is such that after running make test, the commands for "meme" and "mast" pass (these are the two suites that we will use)

**PERL**
Hopefully you have Perl! If not, this is another pre-requisite, although usually it comes pre-installed. Here is the reference: https://www.perl.org/get.html
You will also need to install the Perl dependencies. To check which ones are missing, type the following in the meme directory:
```
cd scripts
perl dependencies.pl
```
You will see a lot of things missing most likely. This is what I did to fix the issues, might vary depending on your machine:
```
sudo apt-get install libexpat-dev
sudo cpan Math::CDF
sudo cpan File::Which
sudo cpan HTML::Template
sudo cpan JSON
sudo cpan XML::Simple
sudo cpan XML::Parser::Expat
```

**zlib**
This one is a little weird to install - try the following commands, if neither work you may need to look up how to install zlib separately (on Della, had to install it manually from a download link):
```
sudo port install zlib
```
```
sudo apt-get install zlib1g-dev
```
**Other things that helped**
This was on UNIX, I got them by going to the scripts directory and running `perl dependencies.pl`

**Trying again**
After downloading necessary pre-requisites, try the installation again by running
```
make test
make install
```

Finally, keep track of where the meme installation ultimately got installed. You will enter this into the *config.yaml* file so that the miner can run MEME automatically.

#### GitHub Codebase:
Open a terminal and type the following commands in a terminal:
```
git clone https://github.com/nathanalam/lassomining.git
```

#### Python libraries
To install the python libraries in the virtual environment, type the following:
```
cd lassomining
pip3 install virtualenv
python3 -m venv venv
```
The command line should now have a (lassomining) tag behind it.

Then, once inside the directory, activate the virtual environment:
```
source venv/bin/activate
```
You will see a (lassomining) identifier in the command line now, indicating that you are in the virtual environment and are safe to freely download libraries without making global changes to your python installation. Now, you can install the necessary libraries (you can copy and paste the following into the command prompt):
```
pip install pip-tools
pip-compile
pip install -r requirements.txt
```
## Operation
You can run this on a local computer, or you can run this using a [SLURM script](https://researchcomputing.princeton.edu/education/online-tutorials/getting-started/introducing-slurm)

### Configuration
The default directories are listed in *config.yml.template* file within the repository. Copy this file and rename it to *config.yml* and change them as you see fit for your use case - the file is commented with more details about what each of the parameters mean.

This file holds parameters which wizkll be specific to each run, and also other variables that you may want to alter.

Alternatively, you can run this with more fine control in the Jupyter notebook entitled "View.ipynb", where you can again specify configuration schemes within the notebook itself.

### Example Usage
After editing the *config.yaml* file, the program can be run as follows:
```
python3 Mine.py
```

## Viewing Results
You can also run and view the results using the jupyter notebook inside the directory. To run Jupyter notebook, enter the following into a command prompt (NOTE: you must be in the virtual environment, or in other words have the "(lassomining)" tag in the command line):
```
jupyter notebook
```
Then, you can enter the notebook and analyze the results as you wish. More instructions and example analyses are available on the notebook.

Sometimes, jupyter can also give issues, and you might get an issue such as "cannot find notebook" when running this command. The following worked for me:
```
sudo apt-get remove ipython
sudo apt-get purge ipython
sudo apt-get autoremove
pip install jupyter
```

And then to run the visualizer, use:
```
python -m IPython notebook
```

From there, the notebook can also guide the mining process - note that the MEME suite is also required for the Jupyter notebook method, and the installation directory of that suite must be recorded.
