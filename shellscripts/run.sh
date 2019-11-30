#!/bin/bash

source shellscripts/setup.sh

python3 Downloader.py
python3 Translator.py
python3 Miner.py

source shellscripts/view.sh