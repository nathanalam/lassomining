#!/bin/bash

source bin/activate
pip3 install flask

export FLASK_APP=Genome_Access_Server.py
flask run