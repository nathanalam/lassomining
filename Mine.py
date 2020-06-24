import os
import time
import json
import shutil
import sys
import yaml

import re
import traceback
import sys
import math
import pandas as pd
import sqlite3
import Bio
import time
from Bio.Seq import Seq, reverse_complement, translate
from Bio.Alphabet import IUPAC

# Reading parameters from the config file
config = None
with open("config.yaml") as file:
    config = yaml.load(file)
    
if (config is None):
    print("Error reading config.yaml")
    sys.exit(1)

print("read the following from config.yaml")
print(config)

runName = config["runName"]
pattern = config["pattern"]
cutoffRank = config["cutoffRank"]
genomeDir = config["genomeDir"]
databaseDir = config["database"]
runDir = config["outputLogs"]
memeJobs = []
for model in config["models"]:
    memeJobs.append([
            model["location"],
            model["numOfMotifs"],
            model["maxWidth"]
        ])

localMotifDir = config["localMotifDir"]
memeDir = config["memeDir"]

'''
Define a function that takes as input the relative path of a FASTA formatted text file, return 
an object that contains a list of sequence objects. Each sequence object has a description field 
["description"] and a sequence field ["sequence"].

From http://www.csbio.sjtu.edu.cn/bioinf/virus-multi/example.htm, specification of a FASTA 
formatted file:
- The first line of each query protein input format must begin with a greater-than (">") symbol 
  in the first column. The word following the ">" symbol is the identifier and description of the 
  sequence, but both are optional.
- The sequence (in single-character code) begins in a different line and ends if another line 
  starting with a ">" appears, which indicates the start of another query protein.
'''
def readFASTA(name, cleanspace = 0):
    descriptions = []
    sequences = []
    sequenceList = []
    tempSequences = []     
        
    with open(name) as file:
        count = -1
        for line in file:
            
            if(line[0] == '>'):
                # if begins with a >, then a description
                descriptions.append(line[1:].replace('\n', ''))
                count += 1
                # skip the first time
                if count > 0 :
                    # combine the tempSequences into a single string and
                    # add it to sequences
                    newSequence = ' '.join(tempSequences)
                    # now remove all of the whitespaces
                    newSequence = newSequence.replace(' ', '')
                    newSequence = newSequence.replace('\n', '')
                    
                    sequences.append(newSequence)
                    # refresh the tempSequence list
                    tempSequences = []
                    
                    sequenceList.append({
                        "description": descriptions[count - 1],
                        "sequence": sequences[count - 1]
                    })
            else:
                tempSequences.append(line)
                
        # combine the tempSequences into a single string and
        # add it to sequences
        newSequence = ' '.join(tempSequences)
        # now remove all of the whitespaces
        newSequence = newSequence.replace(' ', '')
        newSequence = newSequence.replace('\n', '')

        sequences.append(newSequence)
        # refresh the tempSequence list
        tempSequences = []
        
        sequenceList.append({
            "description": descriptions[count],
            "sequence": sequences[count]
        })
                
                
    if len(descriptions) != len(sequences):
        print("ERROR: Number of descriptions does not match number of sequences")
        print("Number of descriptions: " + str(len(descriptions)))
        print("Number of sequences: " + str(len(sequences)))
        sys.exit(1);
        
    print("Read " + str(count + 1) + " objects from FASTA file " + name)
        
    return sequenceList


'''
Define a method that lets us take a sequence and identify motif matches from it. 
Requires the sequence and the directory of the pre-generated MEME motif files as input, 
and returns a tuple containing the B matches and the C matches.
'''
def mastSearch(sequence, memeDir, memeInstall):
    motifMatches = []
    with open("tempseq.txt", "w") as file:
        file.write("> " + "temporary" + "\n")
        file.write(sequence)
        file.close()

    for dir in os.listdir(memeDir):
        command = memeInstall + '/bin/mast -hit_list ' + memeDir + "/" + dir + ' tempseq.txt > tempout' + dir
        # print(command)
        os.system(command)

    for dir in os.listdir(memeDir):
        matchedProts = []
        with open("tempout" + dir, "r") as file:
            inlines = file.readlines()
            inlines = inlines[2:len(inlines) - 1]

        
            for line in inlines:
                # remove ending newline character
                line = line[:len(line) - 1]
                params = line.split(' ')
                while('' in params) : 
                    params.remove('') 
                try:
                    newProt = {
                        "strand" : int(params[1]),
                        "motif" : params[2],
                        "start" : int(params[4]),
                        "end" : int(params[5]),
                        "score" : float(params[6]),
                        "p-value" : float(params[7]),
                        "memeDir": dir[0:len(dir) - 11]
                    }
                    matchedProts.append(newProt)
                except:
                    print("error in parsing line - " + line)
                    print("params: " + str(params))
            file.close()
            os.remove("tempout" + dir)
            motifMatches.append(matchedProts)

    os.remove("tempseq.txt")

    return motifMatches

### Some helper functions
def isOverlapping(start1, end1, start2, end2):
    if (start1 <= start2) and (end1 >= start2):
        return True
    if(start2 <= start1) and (end2 >= start1):
        return True
    
    
    return False

def adjustRangeByORF(ORF, length, start, end):
    if ORF == 2:
        start += 1
        end += 1
    elif ORF == 3:
        start += 1
        end += 1
    elif ORF == -1:
        temp = start
        start = length - end
        end = length - temp
    elif ORF == -2:
        temp = start
        start = length - end - 1
        end = length - temp - 1
    elif ORF == -3:
        temp = start
        start = length - end - 2
        end = length - temp - 2
    
    return [start, end]

'''
Takes in 6 sequences corresponding to 6 ORFs of a given sequence, and then identifies 
B/C gene clusters within them. Then, use the python regular expression library to identify 
precursor peptides matching the regular expression pattern at the top of the script. 
The function returns a list of matched proteins, which have a specific sequence, ORF, 
nearest B/C cluster, and range within the overall sequence.
'''
def patternMatch(sequenceORFs, pattern, filenam, runName, cutoffRank, memeInstall, motifDir):
    Aproteins = []
    AuxProteins = []

    for i in range(0, len(os.listdir(motifDir))):
        AuxProteins.append([])
    
    
    ## generate all motif match sets to go into AuxProteins
    ORFs = [1, 2, 3, -1, -2, -3]
    ORF = 0
    for pair in sequenceORFs:
        overallSequence = pair["sequence"]
        motifMatches = mastSearch(overallSequence, motifDir, memeInstall)
        
        index = 0
        for matchSet in motifMatches:
            # adjust the matchSet for the current ORF
            for b in matchSet:
                b["ORF"] = ORFs[ORF]
                prange = adjustRangeByORF(ORFs[ORF], len(overallSequence) * 3, b["start"] * 3, b["end"] * 3)
                b["start"] = prange[0]
                b["end"] = prange[1]

            AuxProteins[index].extend(matchSet)
            index += 1
    
        ORF += 1
        
    ## create all of the cluster points, and give them a score
    ORF = 0
    for pair in sequenceORFs:
        overallSequence = pair["sequence"]
        description = pair["description"]
        # find all matches in protein that match
        matchIter = re.finditer(pattern, overallSequence)
        done_looping = False
        while not done_looping:
            try:
                match = next(matchIter)
            except StopIteration:
                done_looping = True
            else:
                # get the correct range based on span
                indices = list(match.span())
                indices = adjustRangeByORF(ORFs[ORF], len(overallSequence) * 3, indices[0] * 3, indices[1] * 3)

                # make the ranking calculation and find closest B and C
                start = indices[0]
                end = indices[1]
                def sortFunct(prot):
                    return (prot["start"] - start) ** 2;
                
                rank = 1

                closestProts = []
                closestProtLists = []

                
                for IProteins in AuxProteins:
                    if (rank is 0):
                        continue
                    termE = 1
                    closest = float("inf")
                    closestI = None
                    closestIs = []

                    # create a symbol table linking motifs to arrays of proteins
                    motifTable = {}
                    for prot in IProteins:
                        if isOverlapping(start, end, prot["start"], prot["end"]):
                            # print(str(start) + "-" + str(end) + " |B at " + str(prot["start"]) + "-" + str(prot["end"]))
                            # print(prot)
                            continue
                        if not prot["motif"] in motifTable:
                            motifTable[prot["motif"]] = []
                        motifTable[prot["motif"]].append(prot)
                        closestIs.append(prot)

                    # iterate over each symbol table value to get a summation, multiply those together
                    for motif in motifTable:
                            
                        termi = 0
                        for prot in motifTable[motif]:
                            
                            # skip if not in same direction
                            if not ((prot["ORF"] > 0 and ORFs[ORF] > 0) or (prot["ORF"] < 0 and ORFs[ORF] < 0)):
                                continue
                            
                            diffsquared = (prot["start"] - start) ** 2
                            diffsquared = diffsquared * prot["p-value"]
                            if diffsquared < closest:
                                closestI = prot
                                closest = diffsquared
                            termi += (1.0 / float(diffsquared))
                        
                        termE = termE * termi
                    rank = rank * termE
                    closestIs.sort(key=sortFunct)
                    closestProts.append(closestI)
                    closestProtLists.append(closestIs[0:10])

                if rank <= 0:
                    continue
                else:
                   rank = math.log(rank, 10)

                descriptors = description.split()
                # append the protein to the list of proteins
                if (rank >= cutoffRank):
                    print("Found peptide " + match.group(0) + ", rank " + str(rank) +", in " + filenam)
                    Aproteins.append({
                        "description": description,
                        "sequence": match.group(0),
                        "searchPattern": match.re.pattern,
                        "searchRange": indices,
                        "overallLength": len(overallSequence) * 3,
                        "rank": rank,
                        "closestProts": closestProts,
                        "closestProtLists": closestProtLists,
                        "ORF": ORFs[ORF],
                        "genome": descriptors[1] + " " + descriptors[2],
                        "index": descriptors[0],
                        "runName": runName
                        ## "overallString": match.string
                    })
                
        ORF += 1
    return Aproteins


def scanGenome(runName, pattern, cutoffRank, databaseDir, memeInstall, genomeDir, motifDir):
    # create a database file if one doesn't already exist
    databaseFolder = databaseDir.split("/")
    databaseFolder = "/".join(databaseFolder[0:len(memeName) - 1])
    
    if not os.path.exists(databaseFolder):
        print("creating database directory " + databaseFolder)
        os.makedirs(databaseFolder)
    if not os.path.exists(databaseDir):
        print("Could not find " + databaseDir + ", attempting to create...")
        os.mknod(databaseDir)
    
    conn = sqlite3.connect(databaseDir)

    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS lassopeptides
             (sequence text, start integer, end integer, overallLength integer, rank real, orf integer, genome text, accession text, runName text, closestProts text, closestProtLists text)''')

    matchedProteins = []
    
    readSequences = readFASTA(genomeDir)
    if(not (len(readSequences) % 6) == 0):
        print("Error: sequence in file " + genomeDir + " does not have multiple of 6 sequences")
        print("Instead, it has " + str(len(readSequences)))
        raise RuntimeError
    for i in range(0, len(readSequences), 6):
        buffer = patternMatch(readSequences[i: i + 6], pattern, genomeDir, runName, cutoffRank, memeInstall, motifDir)
        print("Found " + str(len(buffer)) + " peptides in this set of ORFs")
        for peptide in buffer:
            print("Inserting " + peptide['sequence'] + " into sqlite database")
            c.execute("INSERT INTO lassopeptides VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", 
                [peptide['sequence'],
                peptide['searchRange'][0],
                peptide['searchRange'][1],
                peptide['overallLength'],
                peptide['rank'],
                peptide['ORF'],
                peptide['genome'],
                peptide['index'],
                peptide['runName'],
                json.dumps(str(peptide['closestProts'])),
                json.dumps(str(peptide['closestProtLists']))]
            )
            if(toFirestore):
                print("Inserting " + peptide['sequence'] + " into Firestore")
                data = {
                    "sequence": peptide['sequence'], 
                    "start": peptide['searchRange'][0], 
                    "end": peptide['searchRange'][1], 
                    "overallLength": peptide['overallLength'], 
                    "rank": peptide['rank'], 
                    "orf": peptide['ORF'], 
                    "genome": peptide['genome'], 
                    "accession": peptide['index'], 
                    "runName": peptide['runName'], 
                    "closestProts": json.dumps(str(peptide['closestProts'])), 
                    "closestProtLists": json.dumps(str(peptide['closestProtLists']))
                }
        
        matchedProteins.extend(buffer)

    submitted = False
    while(not submitted):
        try:
            conn.commit()
            submitted = True
        except:
            print(databaseDir + " is busy, waiting 5 seconds")
            time.sleep(5)
    conn.close()
    
    return matchedProteins

# takes a three letter codon, and returns the corresponding amino acid or 'X' if unknown
def translate_codon(codonString):
    if(not len(codonString) == 3):
        raise InputError()
        
    try:
        return translate(codonString, to_stop=False, table = 11)
    except:
        return 'X'

# An adapter function for the biopython's translate, takes in a DNA sequence and 
# returns a list of protein sequences. If reading fails, reads each character manually
# and insert X if unable to translate properly
def get_orfs(DNAseq):
    AAList = []
    
    codonArr = []
    seqLen = len(DNAseq) - (len(DNAseq) % 3)
    seq = ''
    try:
        seq = translate(DNAseq[0:seqLen])
    except:
        for i in range(0, seqLen, 3):
            codonArr.append(translate_codon(DNAseq[i:i + 3]))
        seq = ''.join(codonArr)
    AAList.append({
        "ORF": 1,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(DNAseq) - ((len(DNAseq) - 1) % 3)
    seq = ''
    try:
        seq = translate(DNAseq[1:seqLen])
    except:
        for i in range(1, seqLen, 3):
            codonArr.append(translate_codon(DNAseq[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({
        "ORF": 2,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(DNAseq) - ((len(DNAseq) - 2) % 3)
    seq = ''
    try:
        seq = translate(DNAseq[2:seqLen])
    except:
        for i in range(2, seqLen, 3):
            codonArr.append(translate_codon(DNAseq[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({
        "ORF": 3,
        "sequence": seq
    })
    
    backwards_dna = reverse_complement(DNAseq)
    codonArr = []
    seqLen = len(backwards_dna) - (len(backwards_dna) % 3)
    seq = ''
    try:
        seq = translate(backwards_dna[0:seqLen])
    except:
        for i in range(0, seqLen, 3):
            codonArr.append(translate_codon(backwards_dna[i:i + 3]))
        seq = ''.join(codonArr)
    AAList.append({
        "ORF": -1,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(backwards_dna) - ((len(backwards_dna) - 1) % 3)
    seq = ''
    try:
        seq = translate(backwards_dna[1:seqLen])
    except:
        for i in range(1, seqLen, 3):
            codonArr.append(translate_codon(backwards_dna[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({
        "ORF": -2,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(backwards_dna) - ((len(backwards_dna) - 2) % 3)
    seq = ''
    try:
        seq = translate(backwards_dna[2:seqLen])
    except:
        for i in range(2, seqLen, 3):
            codonArr.append(translate_codon(backwards_dna[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({
        "ORF": -3,
        "sequence": seq
    })
    
    return AAList   

def mine(genomeFolder, runName, pattern, cutoffRank, databaseDir, memeInstall, motifDir):

    ## translate the downloaded file into amino acid sequence
    count = 0
    
    if not os.path.exists(genomeFolder):
        print("could not find " + genomeFolder + ", attempting to make it")
        os.makedirs(genomeFolder)

    print("translating fna files in directory folder " + genomeFolder)
    ALLDIRNAMES = os.listdir(genomeFolder)
    for dirname in ALLDIRNAMES:
        translatedDirectory = ""
        if(((dirname[len(dirname) - 3:] == "fna") or (dirname[len(dirname) - 5:] == "fasta")) and not (dirname[:len(dirname) - 3] + "faa") in ALLDIRNAMES):
            print("Opening up " + dirname + " and converting into peptide sequences...")
            DNAseqs = []
            seqDescriptions = []
            try:
                for fastaobj in readFASTA(genomeFolder + dirname):
                    DNAseqs.append(fastaobj["sequence"])
                    seqDescriptions.append(fastaobj["description"])
            except:
		
                continue

            try:
                os.remove(genomeFolder + dirname)
            except:
                continue
                
            entries = []
            for i in range(0, len(DNAseqs)):
                print("converting " + str(len(DNAseqs[i])) + " base pairs from " + seqDescriptions[i])
                aalist = get_orfs(DNAseqs[i])
                print("created " + str(len(aalist)) + " peptide sequences from " + seqDescriptions[i])
                for e in range(0, len(aalist)):
                    entries.append({
                        "sequence": aalist[e]["sequence"],
                        "description": str(seqDescriptions[i] + " - ORF " + str(aalist[e]["ORF"])) 
                    })
            suffixNum = 5
            if(dirname[len(dirname) - 3:] == "fna"):
                suffixNum = 3

            translatedDirectory = genomeFolder + dirname[:len(dirname) - suffixNum] + "faa"
            
            print("writing read peptides into '" + translatedDirectory + "'")
            with open(translatedDirectory, 'w') as outfile:
                for ent in entries:
                    outfile.write("> " + ent["description"] + "\n")
                    outfile.write(ent["sequence"] + "\n\n")
        else:
            continue
        # launch the actual mining of the translated genomes
        print("scanning " + dirname + " for lassos")
        results = scanGenome(runName, pattern, cutoffRank, databaseDir, memeInstall, translatedDirectory, motifDir)
        count += len(results)
        print("found " + str(count) + " peptides")

        ## clear the genomes subdirectory
        print("removing " + translatedDirectory)
        os.remove(translatedDirectory)
        

    return count

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
        modelDir = "/".join(memeName[0: len(memeName) - 1]) + "/"
        memeName = memeName[len(memeName) - 1]
        model = memeName
        memeName = memeName[0: len(memeName) - 4]
        nmotifs = memeJob[1]
        width = memeJob[2]
        print("creating meme motifs for " + memeName)
        print("reading from " + modelDir + model)
        command = memeDir + "/bin/meme -nmotifs " + str(nmotifs) + " -maxw " + str(width) + " " + modelDir + model + " -o " + localMotifDir + memeName
        print(command)
        os.system(command)
    
        os.rename(localMotifDir + memeName + "/meme.txt", localMotifDir + memeName + "Results.txt")
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
        peptideCount += mine(genomeDir, runName, pattern, cutoffRank, databaseDir, memeDir, localMotifDir)
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