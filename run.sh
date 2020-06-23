#!/bin/bash
source /home/nalam/.bash_profile
source /home/nalam/.bashrc
source /home/nalam/lassomining/bin/activate
pip3 install django
pip3 install biopython
pip3 install pandas
source /home/nalam/.bash_profile

echo FINISHED INSTALLING PYTHON LIBRARIES

echo ATTEMPTING TO RUN Mine.py

/home/nalam/lassomining/bin/python3 /home/nalam/lassomining/Mine.py -name refSeqMine -cutoff -20
