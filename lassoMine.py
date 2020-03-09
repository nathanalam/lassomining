import sys
import os
import Bio
from Bio.Seq import Seq, reverse_complement, translate
from Bio.Alphabet import IUPAC
import re
import json
import math
import pandas as pd
import sqlite3

# Major parameters:
pattern ="M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*"
cutoffRank = -100


## Read as a FASTA file
descriptions = []
sequences = []
sequenceList = []
tempSequences = []     
    
count = -1
for line in sys.stdin:
        
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
    sys.stderr.write("ERROR: Number of descriptions does not match number of sequences")
    sys.stderr.write("Number of descriptions: " + str(len(descriptions)))
    sys.stderr.write("Number of sequences: " + str(len(sequences)))
    sys.exit(1);
    
sys.stderr.write("Read " + str(count + 1) + " objects from FASTA file.")
    
# Now the sequenceList has the FASTA read objects

# Now, translate the DNA sequences into AAs
        
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

'''
Now, convert each of the FNA files found in ALLDIRNAMES into corresponding FAA files
'''

DNAseqs = []
seqDescriptions = []

for fastaobj in sequenceList:
    DNAseqs.append(fastaobj["sequence"])
    seqDescriptions.append(fastaobj["description"])
            
entries = []
for i in range(0, len(DNAseqs)):
    sys.stderr.write("converting " + str(len(DNAseqs[i])) + " base pairs from " + seqDescriptions[i])
    aalist = get_orfs(DNAseqs[i])
    sys.stderr.write("created " + str(len(aalist)) + " peptide sequences from " + seqDescriptions[i])
    for e in range(0, len(aalist)):
        entries.append({
            "sequence": aalist[e]["sequence"],
            "description": str(seqDescriptions[i] + " - ORF " + str(aalist[e]["ORF"])) 
        })
        



# Now take the translations in entries

'''
Define a method that lets us take a sequence and identify B and C proteins from it. 
Requires the sequence and the directory of the pre-generated MEME motif files as input, 
and returns a tuple containing the B matches and the C matches.
'''
def mastSearch(sequence, memeDirB, memeDirC):
    Bproteins = []
    Cproteins = []
    with open("tempseq.txt", "w") as file:
        file.write("> " + "temporary" + "\n")
        file.write(sequence)
        file.close()

    os.system('mast -hit_list ' + memeDirB + ' tempseq.txt > tempoutB.txt')
    os.system('mast -hit_list ' + memeDirC + ' tempseq.txt > tempoutC.txt')

    with open("tempoutB.txt", "r") as file:
        inlines = file.readlines()
        inlines = inlines[2:len(inlines) - 1]

        
        for line in inlines:
            # remove ending newline character
            line = line[:len(line) - 1]
            params = line.split(' ')
            while('' in params) : 
                params.remove('') 
            try:
                newB = {
                    "strand" : int(params[1]),
                    "motif" : params[2],
                    "start" : int(params[4]),
                    "end" : int(params[5]),
                    "score" : float(params[6]),
                    "p-value" : float(params[7])
                }
                Bproteins.append(newB)
            except:
                sys.stderr.write("error in parsing line - " + line)
                sys.stderr.write("params: " + str(params))
        file.close()
    with open("tempoutC.txt", "r") as file:
        inlines = file.readlines()
        inlines = inlines[2:len(inlines) - 1]

        
        for line in inlines:
            # remove ending newline character
            line = line[:len(line) - 1]
            params = line.split()
            while('' in params) : 
                params.remove('') 
            try:
                newC = {
                    "strand" : int(params[1]),
                    "motif" : params[2],
                    "start" : int(params[4]),
                    "end" : int(params[5]),
                    "score" : float(params[6]),
                    "p-value" : float(params[7])
                }
                Cproteins.append(newC)
            except:
                sys.stderr.write("error in parsing line - " + line)
                sys.stderr.write("params: " + str(params))
        file.close()

    os.remove("tempseq.txt")
    os.remove("tempoutB.txt")
    os.remove("tempoutC.txt")

    return (Bproteins, Cproteins)

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
def patternMatch(sequenceORFs, pattern, cutoffRank):
    Aproteins = []
    Bproteins = []
    Cproteins = []
    
    
    ## generate all of the Bs and Cs, and label them by ORF
    ORFs = [1, 2, 3, -1, -2, -3]
    ORF = 0
    for pair in sequenceORFs:
        overallSequence = pair["sequence"]
        pair = mastSearch(overallSequence, "motifs/memeb.txt", "motifs/memec.txt")
        
        newBs = pair[0]
        newCs = pair[1]
        
        for b in newBs:
            b["ORF"] = ORFs[ORF]
            prange = adjustRangeByORF(ORFs[ORF], len(overallSequence) * 3, b["start"] * 3, b["end"] * 3)
            b["start"] = prange[0]
            b["end"] = prange[1]
            
            
        for c in newCs:
            c["ORF"] = ORFs[ORF]
            prange = adjustRangeByORF(ORFs[ORF], len(overallSequence) * 3, c["start"] * 3, c["end"] * 3)
            c["start"] = prange[0]
            c["end"] = prange[1]
        
        Bproteins.extend(newBs)
        Cproteins.extend(newCs)
        
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
                rank = 0
                start = indices[0]
                end = indices[1]
                def sortFunct(prot):
                    return (prot["start"] - start) ** 2;

                term1 = 1
                closest = float("inf")
                closestB = None
                closestBs = []
                
                if len(Cproteins) != 0:

                    # create a symbol table linking motifs to arrays of proteins
                    motifTable = {}
                    for prot in Bproteins:
                        if isOverlapping(start, end, prot["start"], prot["end"]):
                            # print(str(start) + "-" + str(end) + " |B at " + str(prot["start"]) + "-" + str(prot["end"]))
                            # print(prot)
                            continue
                        if not prot["motif"] in motifTable:
                            motifTable[prot["motif"]] = []
                        motifTable[prot["motif"]].append(prot)
                        closestBs.append(prot)

                    # iterate over each symbol table value to get a summation, multiply those together
                    for motif in motifTable:
                            
                        termi = 0
                        for prot in motifTable[motif]:
                            
                            diffsquared = (prot["start"] - start) ** 2
                            diffsquared = diffsquared * prot["p-value"]
                            if diffsquared < closest:
                                closestB = prot
                                closest = diffsquared
                            termi += (1.0 / float(diffsquared))
                        
                        term1 = term1 * termi
                            
                closestBs.sort(key=sortFunct)

                term2 = 1
                closest = float("inf")
                closestC = None
                closestCs = []
                if len(Bproteins) != 0:

                    # create a symbol table linking motifs to arrays of proteins
                    motifTable = {}
                    for prot in Cproteins:
                        if isOverlapping(start, end, prot["start"], prot["end"]):
                            # print(str(start) + "-" + str(end) + " |B at " + str(prot["start"]) + "-" + str(prot["end"]))
                            # print(prot)
                            continue
                        if not prot["motif"] in motifTable:
                            motifTable[prot["motif"]] = []
                        motifTable[prot["motif"]].append(prot)
                        closestCs.append(prot)

                    # iterate over each symbol table value to get a summation, multiply those together
                    for motif in motifTable:
                            
                        termi = 0
                        for prot in motifTable[motif]:
                            
                            diffsquared = (prot["start"] - start) ** 2
                            diffsquared = diffsquared * prot["p-value"]
                            if diffsquared < closest:
                                closestC = prot
                                closest = diffsquared
                            termi += (1.0 / float(diffsquared))
                        
                        term2 = term2 * termi

                closestCs.sort(key=sortFunct)

                rank = term1 * term2
                if rank <= 0:
                    continue
                else:
                   rank = math.log(rank, 10)

                descriptors = description.split()
                # append the protein to the list of proteins
                if (rank >= cutoffRank):
                    Aproteins.append({
                        "description": description,
                        "sequence": match.group(0),
                        "searchPattern": match.re.pattern,
                        "searchRange": indices,
                        "overallLength": len(overallSequence) * 3,
                        "rank": rank,
                        "closestB": closestB,
                        "closestBs": closestBs[0:10],
                        "closestC": closestC,
                        "closestCs": closestCs[0:10],
                        "ORF": ORFs[ORF],
                        "genome": descriptors[1] + " " + descriptors[2],
                        "index": descriptors[0]
                        ## "overallString": match.string
                    })
                
        ORF += 1
    return Aproteins


matchedProteins = []


for i in range(0, len(entries), 6):
    buffer = patternMatch(entries[i: i + 6], pattern, cutoffRank) 
    matchedProteins.extend(buffer)

## Matched proteins has all of the information associated with the found lasso peptides
# print(json.dumps(matchedProteins))

## Output the matched proteins in a text format
for protein in matchedProteins:
    print("> " + protein["sequence"])
    print("Rank:" + str(protein["rank"]))
    print(protein["genome"])
    print(protein["description"])
    print("ORF: " + str(protein["ORF"]))
    print("Range: " + str(protein["searchRange"]))
    index = protein["description"].split()[0]
    start = protein["searchRange"][0]
    end = protein["searchRange"][1]
    print("https://www.ncbi.nlm.nih.gov/nuccore/" + index + "?report=genbank&from=" + str(start) + "&to=" + str(end))