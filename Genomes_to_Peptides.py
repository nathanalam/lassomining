import re
import sys
import os
import json
import requests
from Bio.Seq import Seq, reverse_complement, translate
from Bio.Alphabet import IUPAC
import pandas as pd

PATTERN = 'M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*'
# PATTERN = 'CC.CGCCC...TGGC.'
# PATTERN = '.*'

# takes a three letter codon, and returns the corresponding amino acid or 'X' if unknown
def translate_codon(codonString):
    if(not len(codonString) == 3):
        raise InputError()
        
    try:
        return translate(codonString, to_stop=False, table = 11)
    except:
        return 'X'

## An adapter function for the biopython's translate, takes in a DNA sequence and returns a list of protein sequences
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
    AAList = []
    
    codonArr = []
    seqLen = len(DNAseq) - (len(DNAseq) % 3)
    for i in range(0, seqLen, 3):
        codonArr.append(translate_codon(DNAseq[i:i + 3]))
    seq = ''.join(codonArr)
    AAList.append({
        "ORF": 1,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(DNAseq) - ((len(DNAseq) - 1) % 3)
    for i in range(1, seqLen, 3):
        codonArr.append(translate_codon(DNAseq[i:i + 3]))

    seq = ''.join(codonArr)
    AAList.append({
        "ORF": 2,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(DNAseq) - ((len(DNAseq) - 2) % 3)
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
    for i in range(0, seqLen, 3):
        codonArr.append(translate_codon(backwards_dna[i:i + 3]))
    seq = ''.join(codonArr)
    AAList.append({
        "ORF": -1,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(backwards_dna) - ((len(backwards_dna) - 1) % 3)
    for i in range(1, seqLen, 3):
        codonArr.append(translate_codon(backwards_dna[i:i + 3]))

    seq = ''.join(codonArr)
    AAList.append({
        "ORF": -2,
        "sequence": seq
    })
    
    codonArr = []
    seqLen = len(backwards_dna) - ((len(backwards_dna) - 2) % 3)
    for i in range(2, seqLen, 3):
        codonArr.append(translate_codon(backwards_dna[i:i + 3]))

    seq = ''.join(codonArr)
    AAList.append({
        "ORF": -3,
        "sequence": seq
    })
    
    return AAList
    
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

def patternMatch(overallSequence, pattern, description):
    matchedProteinList = []
    
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
            ORF = int(description[len(description) - 2:])
            if ORF == 2:
                indices[0] += 1
                indices[1] += 1
            elif ORF == 3:
                indices[0] += 1
                indices[1] += 1
            elif ORF == -1:
                indices[0] = len(overallSequence) - indices[0]
                indices[1] = len(overallSequence) - indices[1]
            elif ORF == -2:
                indices[0] = len(overallSequence) - indices[0] - 1
                indices[1] = len(overallSequence) - indices[1] - 1
            elif ORF == -3:
                indices[0] = len(overallSequence) - indices[0] - 2
                indices[1] = len(overallSequence) - indices[1] - 2
            matchedProteinList.append({
                "description": description,
                "sequence": match.group(0),
                "searchPattern": match.re.pattern,
                "searchRange": indices,
                "overallLength": len(overallSequence),
                "genome": description[:description.index("/")]
                ## "overallString": match.string
            })
    return matchedProteinList

ALLDIRNAMES = []
for dirname in os.listdir("genomes"):
    ## if a regular file, just add to directory
    if (dirname.find(".") != -1):
        ALLDIRNAMES.append("genomes/" + dirname)
    else:
        for filename in os.listdir("genomes/" + dirname):
            ALLDIRNAMES.append("genomes/" + dirname + "/" + filename)

for dirname in ALLDIRNAMES:
    if((dirname[len(dirname) - 3:] == "fna") and not (dirname[:len(dirname) - 3] + "faa") in ALLDIRNAMES):
        print("Opening up " + dirname + " and converting into peptide sequences...")
        DNAseqs = []
        seqDescriptions = []
        for fastaobj in readFASTA(dirname):
            DNAseqs.append(fastaobj["sequence"])
            seqDescriptions.append(fastaobj["description"])
            
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
        
        print("writing read peptides into '" + dirname[len('genomes/'):len(dirname) - 3] + "faa'")
        with open(dirname[:len(dirname) - 3] + "faa", 'w') as outfile:
            for ent in entries:
                outfile.write("> " + ent["description"] + "\n")
                outfile.write(ent["sequence"] + "\n\n")

DIRNAMES = []
for dirname in os.listdir("genomes"):
    if (dirname.find(".") != -1):
        if(dirname[len(dirname) - 3:] == "faa"):
            DIRNAMES.append("genomes/" + dirname)
    else:
        for filename in os.listdir("genomes/" + dirname):
            if(filename[len(filename) - 3:] == "faa"):
                DIRNAMES.append("genomes/" + dirname + "/" + filename)
    
print("\n***********************************\nReading lasso peptides from the following files:\n***********************************");
print(DIRNAMES)

matchedProteins = []
for filename in DIRNAMES:
    readSequences = readFASTA(filename)
    for seq in readSequences:
        matchedProteins.extend(patternMatch(seq["sequence"], PATTERN, filename[8:len(filename) - 4] + " - " + seq["description"]))
    

print("\n*************************************\nFinished reading peptides")
print("Found " + str(len(matchedProteins)) + " that satisfy the pattern: " + PATTERN)
print("*************************************")

print("Writing output to 'matches.json'")
with open('matches.json', 'w') as outfile:
    json.dump(matchedProteins, outfile)

print("Writing output to 'matches.csv'")
pd.read_json("matches.json").to_csv("matches.csv")