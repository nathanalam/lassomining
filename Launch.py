import argparse
import os
from Miner import mine
import time
import json
import os
import shutil
import sys


## Locations for dependencies - memeDir points to installation location of meme suite, localmotifDir is a temp folder
memeDir = "/home/nalam/meme"
localmotifDir = "/home/nalam/lassomining/motifs/"

# Method to convert a string into a list of meme jobs
def meme_job(string):
    jobs = string.split(",")
    jobList = []
    for j in jobs:
        s = j.split(":")
        file = s[0]
        nmotif = s[1]
        width = s[2]
        jobList.append([file, int(nmotif), int(width)])
    return jobList

# Optional command line arguments with default values
parser = argparse.ArgumentParser()
parser.add_argument('-name', action="store", dest="runName", default="testRun")
parser.add_argument('-pattern', action="store", dest="pattern", default="M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*")
parser.add_argument('-cutoff', action="store", dest="cutoff", default="-100")
parser.add_argument('-genome', action="store", dest="genome", default="/tigress/nalam/genomeMining/genomes/*")
parser.add_argument('-o', action="store", dest="outputDir", default="/tigress/nalam/genomeMining/matches.db")
parser.add_argument('-rundata', action="store", dest="rundata", default="/tigress/nalam/genomeMining/runs/")
parser.add_argument('-model', action="store", dest="model", 
    default="/tigress/nalam/genomeMining/models/b.faa:3:25,/tigress/nalam/genomeMining/models/c.faa:4:25")

# Reading command line arguments
args = parser.parse_args()
runName = args.runName
pattern = args.pattern
cutoffRank = float(args.cutoff)
genome = args.genome
genomeDir = ""
if (genome[len(genome) - 1] == "*"):
    genomeDir = genome[:len(genome) - 1]
databaseDir = args.outputDir
runDir = args.rundata
memeJobs = meme_job(args.model)



print("Meme jobs to be run:")
print(memeJobs)

# Generate motifs and store them in localmotifDir
for memeJob in memeJobs:
    memeName = memeJob[0].split("/")
    modelDir = "/".join(memeName[0: len(memeName) - 1]) + "/"
    memeName = memeName[len(memeName) - 1]
    model = memeName
    memeName = memeName[0: len(memeName) - 4]
    nmotifs = memeJob[1]
    width = memeJob[2]
    print("creating meme motifs for " + memeName)
    print("reading from " + modelDir + model)
    command = memeDir + "/bin/meme -nmotifs " + str(nmotifs) + " -maxw " + str(width) + " " + modelDir + model + " -o " + localmotifDir + memeName
    print(command)
    os.system(command)

    os.rename(localmotifDir + memeName + "/meme.txt", localmotifDir + memeName + "Results.txt")
    shutil.rmtree(localmotifDir + memeName)

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
if not os.path.exists(runDir + runName + ".json"):
    os.mknod(runDir + runName + ".json")

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
    peptideCount += mine(genomeDir, runName, pattern, cutoffRank, databaseDir, memeDir)
    count += 1
    updateRun("processing" + query, count, peptideCount, query)

print("finished all the runs for " + runName)

# Delete all of the temporary files
for memeJob in memeJobs:
    memeName = memeJob[0].split("/")
    memeName = memeName[len(memeName) - 1]
    memeName = memeName[0: len(memeName) - 4]
    if os.path.exists(localmotifDir + memeName + "Results.txt"):
        os.remove(localmotifDir + memeName + "Results.txt")

