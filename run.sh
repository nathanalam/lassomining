#!/bin/bash

source /home/nalam/lassomining/bin/activate
pip3 install django
pip3 install biopython
pip3 install pandas
source /home/nalam/.bash_profile

echo FINISHED INSTALLING PYTHON LIBRARIES

echo ATTEMPTING TO RUN Launch.py

python3 /home/nalam/lassomining/Launch.py
