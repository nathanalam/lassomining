import argparse
import os
from Miner import mine
import time
import json
import os
import shutil
import sys

runName = "testRun"
pattern = "M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*"
cutoffRank = -100
genomeDir = "/tigress/nalam/genomeMining/genomes/"
databaseDir = "/tigress/nalam/genomeMining/matches.db"
memeDir = "/home/nalam/meme"
motifDir = "/home/nalam/lassomining/motifs/"
runDir = "/tigress/nalam/genomeMining/runs/"
memeJobs = [
    ["/tigress/nalam/genomeMining/models/b.faa", 3, 25],
    ["/tigress/nalam/genomeMining/models/c.faa", 4, 25]
]

print("Meme jobs to be run:")
print(memeJobs)
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
    command = memeDir + "/bin/meme -nmotifs " + str(nmotifs) + " -maxw " + str(width) + " " + modelDir + model + " -o " + motifDir + memeName
    print(command)
    os.system(command)

    os.rename(motifDir + memeName + "/meme.txt", motifDir + memeName + "Results.txt")
    shutil.rmtree(motifDir + memeName)

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
    if os.path.exists(motifDir + memeName + "Results.txt"):
        os.remove(motifDir + memeName + "Results.txt")

