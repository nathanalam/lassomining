import tkinter as tk
from tkinter import filedialog
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
from pathlib import Path
from Bio.Seq import Seq, reverse_complement, translate
from Bio.Alphabet import IUPAC

# Mining specific functions
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

# GUI specific functions

# read in default values
config = {
    "runName": "nameOfRun",
    "pattern": "M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*",
    "cutoffRank": "0",
    "genomeDir": "C:/Users/natha/Desktop/genomeminingdemo/genomes/",
    "database": "C:/Users/natha/Desktop/genomeminingdemo/output/matches.db",
    "outputLogs": "C:/Users/natha/Desktop/genomeminingdemo/output/runs/",
    "models": '''[{'location': '/Users/nalam/Desktop/genomeminingdemo/models/b.faa',
  'maxWidth': '25',
  'numOfMotifs': '3'}, {'location': '/Users/nalam/Desktop/genomeminingdemo/models/c.faa',
  'maxWidth': '25',
  'numOfMotifs': '4'}]
''',
    "memeDir": "/Users/nalam/meme",
    "localMotifDir": "C:/Users/natha/Desktop/genomeminingdemo/temp",
}
try:
    with open("config.yaml") as file:
        config = yaml.load(file)
except Exception as error:
    print("Could not find an adjacent config.yaml")


# Make the GUI features    
root= tk.Tk()
scrollbar = tk.Scrollbar(root)

listbox = tk.Listbox(root, yscrollcommand=scrollbar.set)
for i in range(1000):
    listbox.insert(tk.END, str(i))

scrollbar.config(command=listbox.yview)

tk.Label(root, text="Run Name").grid(row=0)
tk.Label(root, text="Search Pattern").grid(row=1)
tk.Label(root, text="Cutoff Rank").grid(row=2)
tk.Label(root, text="Genome Directory").grid(row=3)
tk.Label(root, text="Database Path").grid(row=4)
tk.Label(root, text="Output Directory").grid(row=5)
tk.Label(root, text="Models").grid(row=6)
tk.Label(root, text="Motif Temporary Space").grid(row=7)
tk.Label(root, text="MEME Install Path").grid(row=8)

def fillGenome():
    currdir = os.getcwd()
    tempdir = filedialog.askdirectory(parent=root, initialdir=currdir, title='Please select a directory')
    e4.delete(0, tk.END)
    e4.insert(10, tempdir + "/")

def fillDatabase():
    currdir = os.getcwd()
    tempdir = filedialog.askopenfilename(parent=root, initialdir=currdir, title='Please select a directory')
    e5.delete(0, tk.END)
    e5.insert(10, tempdir)

def fillOutput():
    currdir = os.getcwd()
    tempdir = filedialog.askdirectory(parent=root, initialdir=currdir, title='Please select a directory')
    e6.delete(0, tk.END)
    e6.insert(10, tempdir + "/")

def fillLocalMotif():
    currdir = os.getcwd()
    tempdir = filedialog.askdirectory(parent=root, initialdir=currdir, title='Please select a directory')
    e8.delete(0, tk.END)
    e8.insert(10, tempdir)

def fillMeme():
    currdir = os.getcwd()
    tempdir = filedialog.askdirectory(parent=root, initialdir=currdir, title='Please select a directory')
    e9.delete(0, tk.END)
    e9.insert(10, tempdir)
    

e1 = tk.Entry(root, width=60)
e1.insert(10, config["runName"])
e2 = tk.Entry(root, width=60)
e2.insert(10, config["pattern"])
e3 = tk.Entry(root, width=60)
e3.insert(10, config["cutoffRank"])
e4 = tk.Entry(root, width=60)
e4.insert(10, config["genomeDir"])
tk.Button(root, text='Browse...', command=fillGenome).grid(row=3, column=2, sticky=tk.W, pady=4)
e5 = tk.Entry(root, width=60)
e5.insert(10, config["database"])
tk.Button(root, text='Browse...', command=fillDatabase).grid(row=4, column=2, sticky=tk.W, pady=4)
e6 = tk.Entry(root, width=60)
e6.insert(10, config["outputLogs"])
tk.Button(root, text='Browse...', command=fillOutput).grid(row=5, column=2, sticky=tk.W, pady=4)
e7 = tk.Text(root, width=60)
e7.insert(tk.END, config["models"])
e8 = tk.Entry(root, width=60)
e8.insert(10, config["localMotifDir"])
tk.Button(root, text='Browse...', command=fillLocalMotif).grid(row=7, column=2, sticky=tk.W, pady=4)
e9 = tk.Entry(root, width=60)
e9.insert(10, config["memeDir"])
tk.Button(root, text='Browse...', command=fillMeme).grid(row=8, column=2, sticky=tk.W, pady=4)

e1.grid(row=0, column=1)
e2.grid(row=1, column=1)
e3.grid(row=2, column=1)
e4.grid(row=3, column=1)
e5.grid(row=4, column=1)
e6.grid(row=5, column=1)
e7.grid(row=6, column=1)
e8.grid(row=7, column=1)
e9.grid(row=8, column=1)

def runGui():
    print("read the following from user input")
    runName = e1.get()
    pattern = e2.get()
    cutoffRank = float(e3.get())
    genomeDir = e4.get()
    databaseDir = e5.get()
    runDir = e6.get()
    memeJobs = []
    models = yaml.load(e7.get("1.0", tk.END))
    for model in models:
        print(model)
        memeJobs.append([
                model["location"],
                model["numOfMotifs"],
                model["maxWidth"]
            ])
    localMotifDir = e8.get()
    memeDir = e9.get()

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

    # Now output all of these entries into CSV files for easier viewing
    
    theStuff = ""
    with open(e6.get() + e1.get() + ".json") as file:
        theStuff = json.load(file)

    print("Total peptides found: " + str(theStuff["peptides"]))
    print("Pattern used: " + str(theStuff["pattern"]))
    print("Cutoff rank: " + str(theStuff["cutoff"]))
    print("Time elapsed: " + str(theStuff["totalTime"] / 3600) + " hours")
    print("Genomes mined: " + str(len(theStuff["input"])))
    print("Average time per genome: " + str(theStuff["totalTime"] / len(theStuff["input"])) + " seconds")

    # regular expression function for regular expression search
    def regexp(expr, item):
        reg = re.compile(expr)
        return reg.search(item) is not None

    conn = sqlite3.connect(e5.get())
    conn.create_function("REGEXP", 2, regexp)
    c = conn.cursor()

    selectionStringGenomes = "SELECT DISTINCT genome FROM lassopeptides WHERE runname is '" + config["runName"] + "'"
    distinctGenomes = []
    for row in c.execute(selectionStringGenomes):
        distinctGenomes.append(row[0])

    c.close()

    print("Number of genomes with lasso peptides: " + str(len(distinctGenomes)))

    conn = sqlite3.connect(e5.get())
    conn.create_function("REGEXP", 2, regexp)
    c = conn.cursor()

    selectionStringGenomes = "SELECT DISTINCT sequence, rank, genome, start, end, accession, closestProtLists FROM lassopeptides WHERE runname is '" + config["runName"] + "'"
    peptides = []
    for row in c.execute(selectionStringGenomes):
        peptides.append({
            "sequence": row[0],
            "rank": row[1],
            "genome": row[2],
            "start": row[3],
            "end": row[4],
            "accession": row[5],
            "closests": row[6]
        })
    c.close()
    print("DISTINCT lasso peptide hits: " + str(len(peptides)))

    genomeDict = {}
    genomeArr = []
    for genome in distinctGenomes:
        peptideArr = []
        runningSum = 0
        for peptide in peptides:
            if(peptide["genome"] == genome):
                peptideArr.append(peptide)
                runningSum += peptide["rank"]
        genomeArr.append({
                "genome": genome,
                "average" : (1.0 * runningSum) / len(peptideArr),
                "count": len(peptideArr)
            }
        )
        genomeDict.update({
            genome: peptideArr
        })
    genomeArr.sort(reverse=True, key=lambda i: i['average'])
    for gen in genomeDict.keys():
        genomeDict[gen].sort(reverse=True, key=lambda i: i["rank"])
    for gen in genomeDict.keys():
        rankList = []
        peptideList = []
        startList = []
        endList = []
        accessionList = []
        closestList = []
        for peptide in genomeDict[gen]:
            rankList.append(peptide["rank"])
            peptideList.append(peptide["sequence"])
            startList.append(peptide["start"])
            endList.append(peptide["end"])
            accessionList.append("https://www.ncbi.nlm.nih.gov/nuccore/" + peptide["accession"] + "?report=genbank&from=" + str(peptide["start"]) + "&to=" + str(peptide["end"]))
            closestList.append(peptide["closests"])
        genomeDict.update({
            gen: {
                "ranks": rankList,
                "sequences": peptideList,
                "starts": startList,
                "ends": endList,
                "urls": accessionList,
                "closests": closestList,
            }
        })

    for i in range(0, len(genomeArr)):
        print("Exporting " + e6.get() + genomeArr[i]["genome"] + ".csv")

        sequences = genomeDict[genomeArr[i]["genome"]]["sequences"]
        ranks = genomeDict[genomeArr[i]["genome"]]["ranks"]
        starts = genomeDict[genomeArr[i]["genome"]]["starts"]
        ends = genomeDict[genomeArr[i]["genome"]]["ends"]
        urls = genomeDict[genomeArr[i]["genome"]]["urls"]
        closests = genomeDict[genomeArr[i]["genome"]]["closests"]

        newDict = {}
        newDict["Sequence"] = sequences
        newDict["Rank"] = ranks
        newDict["Start"] = starts
        newDict["End"] = ends
        newDict["URL"] = urls
        newDict["Closest Motifs"] = closests

        precsv = pd.DataFrame.from_dict(newDict)
        precsv.to_csv(e6.get() + genomeArr[i]["genome"] + ".csv")

tk.Button(root, text='Mine', command=runGui).grid(row=9, column=0, sticky=tk.W, pady=4)

root.mainloop()