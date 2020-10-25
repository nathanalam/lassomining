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
from mining import mine, export_to_csv

# Reading parameters from the config file
config = None
with open("config.yml") as file:
    config = yaml.load(file)

if (config is None):
    print("Error reading config.yml")
    sys.exit(1)

print("read the following from config.yml")
print(config)

runName = config["runName"]
pattern = config["pattern"]
cutoffRank = config["cutoffRank"]
genomeDir = config["genomeDir"]
databaseDir = config["database"]
runDir = config["outputLogs"]
memeJobs = []
for model in config["models"]:
    memeJobs.append(
        [model["location"], model["numOfMotifs"], model["maxWidth"]])

localMotifDir = config["localMotifDir"]
memeDir = config["memeDir"]

print("Beginning run " + runName)
print("cutting off hits below " + str(cutoffRank))
print("searching for pattern " + pattern)
print("Meme jobs to be run:")
print(memeJobs)
print("Genomes being read from " + str(genomeDir))
print("writing output to " + databaseDir)

try:
    # check if the localMotifDir exists
    if not os.path.exists(localMotifDir):
        print("creating a folder " + localMotifDir + " for temporary files")
        os.makedirs(localMotifDir)
    # Generate motifs and store them in localMotifDir
    for memeJob in memeJobs:
        memeName = memeJob[0].split("/")
        modelDir = "/".join(memeName[0:len(memeName) - 1]) + "/"
        memeName = memeName[len(memeName) - 1]
        model = memeName
        memeName = memeName[0:len(memeName) - 4]
        nmotifs = memeJob[1]
        width = memeJob[2]
        print("creating meme motifs for " + memeName)
        print("reading from " + modelDir + model)
        command = memeDir + "/bin/meme -nmotifs " + str(
            nmotifs) + " -maxw " + str(
                width
            ) + " " + modelDir + model + " -o " + localMotifDir + memeName
        print(command)
        os.system(command)

        os.rename(localMotifDir + memeName + "/meme.txt",
                  localMotifDir + memeName + "Results.txt")
        shutil.rmtree(localMotifDir + memeName)
except Exception as error:
    print("An error occured while running MEME")
    traceback.print_tb(sys.exc_info()[2])
    print(str(error))

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

    # create a run file to log the progress of this run
    if not os.path.exists(runDir):
        print("creating output directory " + runDir)
        os.makedirs(runDir)

    if not os.path.exists(runDir + runName + ".json"):
        print("writing output logs to " + runDir + runName + ".json")
        Path(runDir + runName + ".json").touch()

    # create genome folder if not already there
    if not os.path.exists(genomeDir):
        print("could not find " + genomeDir + ", attempting to make it")
        os.makedirs(genomeDir)

    # create output database if not already there
    # create a database file if one doesn't already exist
    path = databaseDir.split("/")
    databaseFolder = "/".join(path[0:len(path) - 1])

    if not os.path.exists(databaseFolder):
        print("creating database directory " + databaseFolder)
        os.makedirs(databaseFolder)
    if not os.path.exists(databaseDir):
        print("Could not find " + databaseDir + ", attempting to create...")
        Path(databaseDir).touch()

    # get the list of queries
    queries = os.listdir(genomeDir)

    # define a function to progressively update the current status of the run


    def updateRun(message, number, count, accession):
        runStatus["phase"] = message
        runStatus["totalTime"] = time.time() - t0
        runStatus["input"].append(accession)
        runStatus["progress"] = str((number * 1.0) / len(queries))
        runStatus["peptides"] = count
        with open(runDir + runName + '.json', 'w+') as outputFile:
            outputFile.write(json.dumps(runStatus))

    count = 0
    peptideCount = 0

    for query in queries:
        peptideCount += mine(genomeDir, runName, pattern, cutoffRank,
                             databaseDir, memeDir, localMotifDir)
        count += 1
        updateRun("processing" + query, count, peptideCount, query)

    print("finished all the runs for " + runName)
except Exception as error:
    print("An error occured while mining")
    traceback.print_tb(sys.exc_info()[2])
    print(str(error))

# Delete all of the temporary MEME files
if os.path.exists(localMotifDir):
    shutil.rmtree(localMotifDir)

# Export the information into CSVs
export_to_csv(config["runName"], config["database"], config["outputLogs"])
