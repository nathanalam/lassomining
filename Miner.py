'''
Author: Nathan Alam

Find the precursor peptides within the .faa files in the genome directory, along with
a ranking and nearest cluster of B and C protein candidates
'''
import re
import sys
import os
import json
import math
import pandas as pd
import sqlite3

PATTERN = 'M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*'
# PATTERN = 'CC.CGCCC...TGGC.'
# PATTERN = '.*'

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
def patternMatch(sequenceORFs, pattern, filenam, runName, cutoffRank, memeInstall):
    Aproteins = []
    AuxProteins = []

    for i in range(0, len(os.listdir("motifs"))):
        AuxProteins.append([])
    
    
    ## generate all motif match sets to go into AuxProteins
    ORFs = [1, 2, 3, -1, -2, -3]
    ORF = 0
    for pair in sequenceORFs:
        overallSequence = pair["sequence"]
        motifMatches = mastSearch(overallSequence, "motifs", memeInstall)
        
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
                    closestProtLists.append(closestIs)

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
                        "closestProts": closestProts,
                        "closestProtLists": closestProtLists[0:10],
                        "ORF": ORFs[ORF],
                        "genome": descriptors[1] + " " + descriptors[2],
                        "index": descriptors[0],
                        "runName": runName
                        ## "overallString": match.string
                    })
                
        ORF += 1
    return Aproteins


def scanGenomes(runName, pattern, cutoffRank, databaseDir, memeInstall):
    conn = sqlite3.connect(databaseDir)

    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS lassopeptides
             (sequence text, start integer, end integer, overallLength integer, rank real, orf integer, genome text, accession text, runName text, closestProts text, closestProtLists text)''')

    DIRNAMES = []
    for dirname in os.listdir("genomes"):
        if (dirname.find(".") != -1):
            if(dirname[len(dirname) - 3:] == "faa"):
                DIRNAMES.append("genomes/" + dirname)
        else:
            for filename in os.listdir("genomes/" + dirname):
                if(filename[len(filename) - 3:] == "faa"):
                    DIRNAMES.append("genomes/" + dirname + "/" + filename)

    matchedProteins = []
    for filename in DIRNAMES:
        readSequences = readFASTA(filename)
        if(not (len(readSequences) % 6) == 0):
            print("Error: sequence in file " + filename + " does not have multiple of 6 sequences")
            print("Instead, it has " + str(len(readSequences)))
            raise RuntimeError
        for i in range(0, len(readSequences), 6):
            buffer = patternMatch(readSequences[i: i + 6], pattern, filename, runName, cutoffRank, memeInstall)
            for peptide in buffer:
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

    conn.commit()
    conn.close()
    
    return matchedProteins

def mine(accession, runName, pattern, cutoffRank, databaseDir, memeInstall):

    accession = accession.strip()
    ## download the genome associated with the accession number
    print("Generating URL File downloads for genomes")
    if(".gz" in accession[len(accession) - 3:]):
        print("found a direct zip file address")
        os.system('echo "' + accession +'" >> fileurls.txt')
    else: 
        print("attempting assembly search")
        str1 = 'esearch -db assembly -query "' + accession + '" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" \'{print $0"/"$NF"_genomic.fna.gz"}\' >> fileurls.txt'
        print(str1)
        os.system(str1)
        print("attempting nucleotide search")
        str1 = 'esearch -db nucleotide -query "' + accession + '" | efetch -format fasta >> genomes/' + accession + '.fna'
        print(str1)
        os.system(str1)
        print("attempting genome search")
        str1 = 'esearch -db genome -query "' + accession + '" | efetch -format fasta >> genomes/' + accession + '.fna'
        print(str1)
        os.system(str1)
        # os.system('esearch -db assembly -query "' + accession + '" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" \'{print $0"/"$NF"_genomic.fna.gz"}\' >> fileurls.txt')

    print("Downloading files using wget into genomes folder")
    os.system("wget --directory-prefix=genomes $( cat fileurls.txt )")

    print("Unzipping downloaded files")
    os.system("gunzip -r genomes")

    os.system("rm fileurls.txt")

    ## translate the downloaded file into amino acid sequence
    print("translating accessions")
    os.system("python3 Translator.py")

    # launch the actual mining of the translated genomes
    print("scanning genomes for lassos")
    results = scanGenomes(runName, pattern, cutoffRank, databaseDir, memeInstall)
    count = len(results)
    print("found " + str(count) + " peptides from " + accession)

     ## clear the genomes subdirectory
    print("clearing the genomes directory...")
    ALLDIRNAMES = []
    for dirname in os.listdir("genomes"):
        ## if a regular file, just add to directory
        if (dirname.find(".") != -1):
            ALLDIRNAMES.append("genomes/" + dirname)
        else:
            for filename in os.listdir("genomes/" + dirname):
                ALLDIRNAMES.append("genomes/" + dirname + "/" + filename)
    for dirname in ALLDIRNAMES:
        os.remove(dirname)

    return count