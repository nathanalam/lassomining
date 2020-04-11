import argparse
import os
from Motifmaker import makeMeme
from Miner import mine
import time
import json

runName = "testRun"
pattern = "M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*"
queries = ["JF784339.1"]
cutoffRank = -100
numMemes = 0
memeJobs = []

# for i in range(1, numMemes + 1):
#     motifName = request.POST.get("motifName" + str(i))
#     motifSeqs = request.POST.get("sequences" + str(i)).split(";")
#     nmotifs = request.POST.get("nmotifs" + str(i))
#     memeJobs.append({
#         "name": motifName,
#         "seqs": motifSeqs,
#         "nmotifs": nmotifs
#     })

print("Meme jobs to be run:")
print(memeJobs)

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
if not os.path.exists("runs/" + runName + ".json"):
    os.mknod("runs/" + runName + ".json")

# define a function to progressively update the current status of the run

def updateRun(message, number, count, accession):
    runStatus["phase"] = message
    runStatus["totalTime"] = time.time() - t0
    runStatus["input"].append(accession)
    runStatus["progress"] = str((number * 1.0) / len(queries))
    runStatus["peptides"] = count
    with open('runs/' + runName + '.json', 'w+') as outputFile:
        outputFile.write(json.dumps(runStatus))

count = 0
peptideCount = 0

# create the new motifs with meme
for memeJob in memeJobs:
    seeq = []
    for p in memeJob["seqs"]:
        pair = p.split(",")
        seeq.append((pair[0], pair[1]))

    makeMeme(seeq, memeJob["name"], memeJob["nmotifs"])

for query in queries:
    peptideCount += mine(query, runName, pattern, cutoffRank, "matches.db", "/home/blucheez/Projects/meme")
    count += 1
    updateRun("processing" + query, count, peptideCount, query)

print("finished all the runs for " + runName)

# Delete all of the temporary files
for memeJob in memeJobs:
    if os.path.exists("motifs/" + memeJob["name"] + "Results.txt"):
        os.remove("motifs/" + memeJob["name"] + "Results.txt")

