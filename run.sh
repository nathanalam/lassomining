#!/bin/bash

source bin/activate
pip3 install flask
pip3 install flask_cors

export FLASK_APP=Genome_Access_Server.py
flask run