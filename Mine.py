import os
import time
import json
import shutil
import sys
import yaml
from pathlib import Path
import re
import traceback
import sys
import math
import pandas as pd
import sqlite3
from mining import mine, export_to_csv, generate_motifs

# Reading parameters from the config file
config = None
with open("config.yml") as file:
    config = yaml.load(file)

if (config is None):
    print("Error reading config.yml")
    sys.exit(1)

runName = config["runName"]
pattern = config["pattern"]
cutoffRank = config["cutoffRank"]
genomeDir = config["genomeDir"]
databaseDir = config["database"]
motifs = config["motifs"]
memeDir = config["memeDir"]
print("Beginning run " + runName)
print("cutting off hits below " + str(cutoffRank))
print("searching for pattern " + pattern)
print("Using these motifs:")
print(motifs)
print("Genomes being read from " + str(genomeDir))
print("writing output to " + databaseDir)

try:
    # start a timer
    t0 = time.time()

    # store meta data about the particular run
    runStatus = {
        "name": runName,
        "pattern": pattern,
        "input": [],
        "progress": 0.0,
        "peptides": 0,
        "cutoff": cutoffRank
    }

    ## create filepaths
    # create genome folder if not already there
    if not os.path.exists(genomeDir):
        print("could not find " + genomeDir + ", attempting to make it")
        os.makedirs(genomeDir)
    # create output database if not already there
    path = databaseDir.split("/")
    databaseFolder = "/".join(path[0:len(path) - 1])
    if not os.path.exists(databaseFolder):
        print("creating database directory " + databaseFolder)
        os.makedirs(databaseFolder)
    if not os.path.exists(databaseDir):
        print("Could not find " + databaseDir + ", attempting to create...")
        Path(databaseDir).touch()

    mine(genomeDir, runName, pattern, cutoffRank, databaseDir, memeDir, motifs)
    print("finished all the runs for " + runName)

except Exception as error:
    print("An error occured while mining")
    traceback.print_tb(sys.exc_info()[2])
    print(str(error))

# Export the information into CSVs
export_to_csv(config["runName"], config["database"], os.path.join(config["database"][:config["database"].rfind('/')], 'csvs/'))
