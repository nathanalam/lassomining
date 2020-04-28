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