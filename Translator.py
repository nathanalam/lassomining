'''
Author: Nathan Alam

Translates nucleic acid sequences in the "genomes" folder into amino acid counterparts.
Reads FNA input and outputs FAA output - each entry in the FNA file gets its own FAA file,
with each entry in the FAA file corresponding to one of six reading frames

'''


import argparse
import os
import Bio
from Bio.Seq import Seq, reverse_complement, translate
from Bio.Alphabet import IUPAC

'''
Add command line arguments
'''
parser = argparse.ArgumentParser()
# add long and short argument
parser.add_argument("--clear", "-c", help="clear all .faa files and translate again", action="store_true")
parser.add_argument('-dir', action="store", dest="dir", default="genomes/")

# read arguments from the command line
args = parser.parse_args()

# check for --clear
clearDir = False
if args.clear:
    clearDir = True

genomeLocation = args.dir

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


# Populate ALLDIRNAMES with all of the files in the genome file
ALLDIRNAMES = []
for dirname in os.listdir(genomeLocation):
    ## if a regular file, just add to directory
    if (dirname.find(".") != -1):
        ALLDIRNAMES.append(genomeLocation + dirname)
    else:
        for filename in os.listdir(genomeLocation + dirname):
            ALLDIRNAMES.append(genomeLocation + dirname + "/" + filename)

if(clearDir):
    print("Clearing all .faa files in genome directory")
    for dirname in ALLDIRNAMES:
        if (dirname[len(dirname) - 3:] == "faa"):
            os.remove(dirname)
        
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

for dirname in ALLDIRNAMES:
    if(((dirname[len(dirname) - 3:] == "fna") or (dirname[len(dirname) - 5:] == "fasta")) and not (dirname[:len(dirname) - 3] + "faa") in ALLDIRNAMES):
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
        
        print("writing read peptides into '" + dirname[len(genomeLocation):len(dirname) - 3] + "faa'")
        with open(dirname[:len(dirname) - 3] + "faa", 'w') as outfile:
            for ent in entries:
                outfile.write("> " + ent["description"] + "\n")
                outfile.write(ent["sequence"] + "\n\n")