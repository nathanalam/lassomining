#!/bin/bash

source bin/activate
pip3 install django
pip3 install biopython
pip3 install pandas
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.1.0:$PATH
source ~/.bash_profile
source ~/.profile

echo FINISHED INSTALLING PYTHON LIBRARIES