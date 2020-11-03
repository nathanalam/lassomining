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

import firebase_admin
from firebase_admin import credentials
from firebase_admin import firestore

# some flags for debugging
REMOVE_GENOMES_ON_TRANSLATE = False
PRINT_EACH_FIND = False
PRINT_EACH_WRITE = False
# put -1 to take all above the cutoff
TAKE_TOP_N = 10
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


def readFASTA(name, cleanspace=0):
    descriptions = []
    sequences = []
    sequenceList = []
    tempSequences = []

    with open(name) as file:
        count = -1
        for line in file:

            if (line[0] == '>'):
                # if begins with a >, then a description
                descriptions.append(line[1:].replace('\n', ''))
                count += 1
                # skip the first time
                if count > 0:
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
        print(
            "ERROR: Number of descriptions does not match number of sequences")
        print("Number of descriptions: " + str(len(descriptions)))
        print("Number of sequences: " + str(len(sequences)))
        sys.exit(1)

    print("Read " + str(count + 1) + " objects from FASTA file " + name)

    return sequenceList


'''
Given a protein sequence, meme motifs, and a meme directory, return a list of open reading
frames in that sequence with associated motifs
'''


def mast_orfs(sequence, motifs, memeInstall, readingFrame):
    def find(s, ch):
        return [i for i, ltr in enumerate(s) if ltr == ch]

    end_indices = find(sequence, '*')
    orfs = []
    start_index = 0
    for end_index in end_indices:
        if (end_index - start_index > 200):
            orfs.append({
                "start": start_index,
                "end": end_index,
                "counts": {},
                "motifs": {}
            })
        start_index = end_index

    with open("tempseq.txt", "w") as file:
        file.write("> " + "temporary" + "\n")
        file.write(sequence)
        file.close()

    for motif in motifs:
        (_, name) = os.path.split(motif)
        command = memeInstall + '/bin/mast -hit_list ' + motif + ' tempseq.txt > tempout' + name
        # print(command)
        os.system(command)
        matched_motifs = []
        with open("tempout" + name, "r") as file:
            inlines = file.readlines()
            inlines = inlines[2:len(inlines) - 1]

            for line in inlines:
                # remove ending newline character
                line = line[:len(line) - 1]
                params = line.split(' ')
                while ('' in params):
                    params.remove('')
                try:
                    newProt = {
                        "strand": int(params[1]),
                        "motif": params[2],
                        "start": int(params[4]),
                        "end": int(params[5]),
                        "score": float(params[6]),
                        "p-value": float(params[7]),
                        "memeDir": motif
                    }
                    matched_motifs.append(newProt)
                except:
                    print("error in parsing line - " + line)
                    print("params: " + str(params))
            file.close()
            os.remove("tempout" + name)

            # assign prots to orfs
            for prot in matched_motifs:
                for orf in orfs:
                    if (orf['start'] <= prot['start']
                            and orf['end'] >= prot['end']):
                        if (prot["memeDir"] in orf['counts']):
                            orf['counts'][prot["memeDir"]] += 1
                        else:
                            orf['counts'][prot["memeDir"]] = 1
                        if (prot["memeDir"] in orf['motifs']):
                            orf['motifs'][prot["memeDir"]].append(prot)
                        else:
                            orf['motifs'][prot["memeDir"]] = [prot]

    os.remove("tempseq.txt")

    matched_orfs = []

    for motif in motifs:
        (_, name) = os.path.split(motif)
        these_orfs = []
        for orf in orfs:
            if (motif in orf['counts']):
                these_orfs.append({
                    'start': orf['start'],
                    'end': orf['end'],
                    'count': orf['counts'][motif],
                    'motifs': orf['motifs'][motif],
                    'readingFrame': readingFrame,
                    'motifType': name,
                    'sequence': sequence[orf['start']:orf['end']]
                })
        matched_orfs.append(
            sorted(these_orfs, key=lambda orf: orf['count'],
                   reverse=True)[1:4])
    return matched_orfs


### Some helper functions
def isOverlapping(start1, end1, start2, end2):
    if (start1 <= start2) and (end1 >= start2):
        return True
    if (start2 <= start1) and (end2 >= start1):
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


def patternMatch(sequenceORFs, pattern, filenam, runName, cutoffRank,
                 memeInstall, motifs):
    # the matched peptides themselves
    Aproteins = []
    minRank = cutoffRank
    minRankIndex = 0
    # an array of arrays of auxillary proteins (B, C, any other motif specified proteins)
    AuxProteins = []

    for i in range(0, len(motifs)):
        AuxProteins.append([])

    ## generate all motif match sets to go into AuxProteins
    ReadingFrames = [1, 2, 3, -1, -2, -3]
    RFindex = 0
    for pair in sequenceORFs:
        overallSequence = pair["sequence"]
        motifMatches = mast_orfs(overallSequence, motifs, memeInstall, RFindex)

        index = 0
        for matchSet in motifMatches:
            # adjust the matchSet for the current ORF
            for match in matchSet:
                match["readingFrame"] = ReadingFrames[RFindex]
                prange = adjustRangeByORF(ReadingFrames[RFindex],
                                          len(overallSequence) * 3,
                                          match["start"] * 3, match["end"] * 3)
                match["start"] = prange[0]
                match["end"] = prange[1]

            AuxProteins[index].extend(matchSet)
            index += 1

        RFindex += 1

    ## find and score each of the potential precursor sequences
    RFindex = 0
    for pair in sequenceORFs:
        overallSequence = pair["sequence"]
        description = pair["description"]
        # find all matches in sequence that match the pattern
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
                indices = adjustRangeByORF(ReadingFrames[RFindex],
                                           len(overallSequence) * 3,
                                           indices[0] * 3, indices[1] * 3)

                # make the ranking calculation and find closest B and C
                start = indices[0]
                end = indices[1]

                def sortFunct(prot):
                    center = (start + end) / 2
                    pcenter = (prot['start'] + prot['end']) / 2
                    return abs(pcenter - center)

                rank = 0

                closestOrfs = []
                """
                Rank the peptide according to the identified motifs
                """
                # go through each match set, I ~= B, C, D, etc.
                for IProteins in AuxProteins:
                    closestIs = []
                    for IProtein in IProteins:
                        center = (start + end) / 2
                        Icenter = (IProtein['start'] + IProtein['end']) / 2
                        if(abs(center - Icenter) < 10000 and not isOverlapping(start, end, IProtein['start'], IProtein['end'])):
                            closestIs.append(IProtein)
                            
                    if (len(closestIs) > 0):
                        closestIs.sort(key=sortFunct)
                        rank += closestIs[0]['count']
                        closestOrfs.append(closestIs[0])

                descriptors = description.split()
                # append the protein to the list of proteins
                if (rank >= cutoffRank):

                    if (PRINT_EACH_FIND):
                        print("Found peptide " + match.group(0) + ", rank " +
                              str(rank) + ", in " + filenam)
                    if (TAKE_TOP_N > 0):
                        if (rank >= minRank):
                            newA = {
                                "description": description,
                                "sequence": match.group(0),
                                "searchPattern": match.re.pattern,
                                "searchRange": indices,
                                "overallLength": len(overallSequence) * 3,
                                "rank": rank,
                                "closestOrfs": closestOrfs,
                                "readingFrame": ReadingFrames[RFindex],
                                "genome":
                                descriptors[1] + " " + descriptors[2],
                                "index": descriptors[0],
                                "runName": runName
                                ## "overallString": match.string
                            }
                            if (len(Aproteins) > TAKE_TOP_N):
                                # kick out the lowest ranked
                                minRank = Aproteins[0]["rank"]
                                minRankIndex = 0
                                ind = 0
                                for aProt in Aproteins:
                                    if (aProt["rank"] < minRank):
                                        minRank = aProt["rank"]
                                        minRankIndex = ind
                                    ind += 1
                                Aproteins[minRankIndex] = newA
                                ind = 0
                                for aProt in Aproteins:
                                    ind += 1
                            else:
                                Aproteins.append(newA)

                        continue
                    newA = {
                        "description": description,
                        "sequence": match.group(0),
                        "searchPattern": match.re.pattern,
                        "searchRange": indices,
                        "overallLength": len(overallSequence) * 3,
                        "rank": rank,
                        "closestOrfs": closestOrfs,
                        "readingFrame": ReadingFrames[RFindex],
                        "genome": descriptors[1] + " " + descriptors[2],
                        "index": descriptors[0],
                        "runName": runName
                        ## "overallString": match.string
                    }
                    Aproteins.append(newA)

        RFindex += 1
    return Aproteins


def scanGenome(runName, pattern, cutoffRank, databaseDir, memeInstall,
               genomeDir, motifs):
    conn = sqlite3.connect(databaseDir)

    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS lassopeptides
             (sequence text, start integer, end integer, overallLength integer, rank real, orf integer, genome text, accession text, runName text, closestOrfs text)'''
              )

    matchedProteins = []

    readSequences = readFASTA(genomeDir)
    if (not (len(readSequences) % 6) == 0):
        print("Error: sequence in file " + genomeDir +
              " does not have multiple of 6 sequences")
        print("Instead, it has " + str(len(readSequences)))
        raise RuntimeError
    for i in range(0, len(readSequences), 6):
        buffer = patternMatch(readSequences[i:i + 6], pattern, genomeDir,
                              runName, cutoffRank, memeInstall, motifs)
        print("Found " + str(len(buffer)) + " peptides in this set of ORFs")
        for peptide in buffer:
            if (PRINT_EACH_WRITE):
                print("Inserting " + peptide['sequence'] +
                      " into sqlite database")
            c.execute(
                "INSERT INTO lassopeptides VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    peptide['sequence'], peptide['searchRange'][0],
                    peptide['searchRange'][1], peptide['overallLength'],
                    peptide['rank'], peptide['readingFrame'],
                    peptide['genome'], peptide['index'], peptide['runName'],
                    json.dumps(str(peptide['closestOrfs']))
                ])

        matchedProteins.extend(buffer)

    submitted = False
    while (not submitted):
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
    if (not len(codonString) == 3):
        raise InputError()

    try:
        return translate(codonString, to_stop=False, table=11)
    except:
        return 'X'


# An adapter function for the biopython's translate, takes in a DNA sequence and
# returns a list of protein sequences. If reading fails, reads each character manually
# and insert X if unable to translate properly
def get_reading_frames(DNAseq):
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
    AAList.append({"ORF": 1, "sequence": seq})

    codonArr = []
    seqLen = len(DNAseq) - ((len(DNAseq) - 1) % 3)
    seq = ''
    try:
        seq = translate(DNAseq[1:seqLen])
    except:
        for i in range(1, seqLen, 3):
            codonArr.append(translate_codon(DNAseq[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({"ORF": 2, "sequence": seq})

    codonArr = []
    seqLen = len(DNAseq) - ((len(DNAseq) - 2) % 3)
    seq = ''
    try:
        seq = translate(DNAseq[2:seqLen])
    except:
        for i in range(2, seqLen, 3):
            codonArr.append(translate_codon(DNAseq[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({"ORF": 3, "sequence": seq})

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
    AAList.append({"ORF": -1, "sequence": seq})

    codonArr = []
    seqLen = len(backwards_dna) - ((len(backwards_dna) - 1) % 3)
    seq = ''
    try:
        seq = translate(backwards_dna[1:seqLen])
    except:
        for i in range(1, seqLen, 3):
            codonArr.append(translate_codon(backwards_dna[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({"ORF": -2, "sequence": seq})

    codonArr = []
    seqLen = len(backwards_dna) - ((len(backwards_dna) - 2) % 3)
    seq = ''
    try:
        seq = translate(backwards_dna[2:seqLen])
    except:
        for i in range(2, seqLen, 3):
            codonArr.append(translate_codon(backwards_dna[i:i + 3]))

        seq = ''.join(codonArr)
    AAList.append({"ORF": -3, "sequence": seq})

    return AAList


def mine(genomeFolder, runName, pattern, cutoffRank, databaseDir, memeInstall,
         motifs):

    ## translate the downloaded file into amino acid sequence
    count = 0

    print("translating fna files in directory folder " + genomeFolder)
    ALLDIRNAMES = os.listdir(genomeFolder)
    for dirname in ALLDIRNAMES:
        translatedDirectory = ""
        if (((dirname[len(dirname) - 3:] == "fna") or
             (dirname[len(dirname) - 5:] == "fasta"))
                and not (dirname[:len(dirname) - 3] + "faa") in ALLDIRNAMES):
            print("Opening up " + dirname +
                  " and converting into peptide sequences...")
            DNAseqs = []
            seqDescriptions = []
            try:
                for fastaobj in readFASTA(genomeFolder + dirname):
                    DNAseqs.append(fastaobj["sequence"])
                    seqDescriptions.append(fastaobj["description"])
            except:

                continue

            if (REMOVE_GENOMES_ON_TRANSLATE):
                try:
                    os.remove(genomeFolder + dirname)
                except:
                    continue

            entries = []
            for i in range(0, len(DNAseqs)):
                print("converting " + str(len(DNAseqs[i])) +
                      " base pairs from " + seqDescriptions[i])
                aalist = get_reading_frames(DNAseqs[i])
                print("created " + str(len(aalist)) +
                      " peptide sequences from " + seqDescriptions[i])
                for e in range(0, len(aalist)):
                    entries.append({
                        "sequence":
                        aalist[e]["sequence"],
                        "description":
                        str(seqDescriptions[i] + " - ORF " +
                            str(aalist[e]["ORF"]))
                    })
            suffixNum = 5
            if (dirname[len(dirname) - 3:] == "fna"):
                suffixNum = 3

            translatedDirectory = genomeFolder + dirname[:len(dirname) -
                                                         suffixNum] + "faa"

            print("writing read peptides into '" + translatedDirectory + "'")
            with open(translatedDirectory, 'w') as outfile:
                for ent in entries:
                    outfile.write("> " + ent["description"] + "\n")
                    outfile.write(ent["sequence"] + "\n\n")
        else:
            continue
        # launch the actual mining of the translated genomes
        print("scanning " + dirname + " for lassos")
        results = scanGenome(runName, pattern, cutoffRank, databaseDir,
                             memeInstall, translatedDirectory, motifs)
        count += len(results)
        print("found " + str(count) + " peptides")

        ## clear the genomes subdirectory
        print("removing " + translatedDirectory)
        os.remove(translatedDirectory)

    return count


'''
Export results of a given runName to CSV
'''


def export_to_csv(run_name, database_dir, output_dir):
    def regexp(expr, item):
        reg = re.compile(expr)
        return reg.search(item) is not None

    print(f'Exporting results of run {run_name}')
    conn = sqlite3.connect(database_dir)
    conn.create_function("REGEXP", 2, regexp)
    c = conn.cursor()

    selectionStringGenomes = "SELECT DISTINCT genome FROM lassopeptides WHERE runname is '" + run_name + "'"
    distinctGenomes = []
    for row in c.execute(selectionStringGenomes):
        distinctGenomes.append(row[0])
    print("Number of genomes with lasso peptides: " +
          str(len(distinctGenomes)))

    selectionStringGenomes = "SELECT DISTINCT sequence, rank, genome, start, end, accession, closestOrfs FROM lassopeptides WHERE runname is '" + run_name + "'"
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
            if (peptide["genome"] == genome):
                peptideArr.append(peptide)
                runningSum += peptide["rank"]
        genomeArr.append({
            "genome": genome,
            "average": (1.0 * runningSum) / len(peptideArr),
            "count": len(peptideArr)
        })
        genomeDict.update({genome: peptideArr})
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
            accessionList.append("https://www.ncbi.nlm.nih.gov/nuccore/" +
                                 peptide["accession"] +
                                 "?report=genbank&from=" +
                                 str(peptide["start"]) + "&to=" +
                                 str(peptide["end"]))
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
        print("Exporting " + output_dir + genomeArr[i]["genome"] + ".csv")

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
        precsv.to_csv(os.path.join(output_dir,
                                   genomeArr[i]["genome"] + ".csv"))


def generate_motifs(meme_jobs, output_dir, meme_dir):
    try:
        # check if the output_dir exists
        if not os.path.exists(output_dir):
            print("creating a folder " + output_dir + " for temporary files")
            os.makedirs(output_dir)
        # Generate motifs and store them in output_dir
        for meme_job in meme_jobs:
            memeName = meme_job['fasta'].split("/")
            modelDir = "/".join(memeName[0:len(memeName) - 1]) + "/"
            memeName = memeName[len(memeName) - 1]
            model = memeName
            memeName = memeName[0:len(memeName) - 4]
            nmotifs = meme_job['num_motifs']
            width = meme_job['max_width']
            print("creating meme motifs for " + memeName)
            print("reading from " + modelDir + model)
            command = meme_dir + "/bin/meme -nmotifs " + str(
                nmotifs) + " -maxw " + str(
                    width
                ) + " " + modelDir + model + " -o " + output_dir + memeName
            print(command)
            os.system(command)

            os.rename(output_dir + memeName + "/meme.txt",
                      output_dir + memeName + "Results.txt")
            shutil.rmtree(output_dir + memeName)
    except Exception as error:
        print("An error occured while generating motifs")
        traceback.print_tb(sys.exc_info()[2])
        print(str(error))

def export_to_firebase(db_dir, run_name, cred_file):
    # regular expression function for regular expression search
    def regexp(expr, item):
        reg = re.compile(expr)
        return reg.search(item) is not None

    # read peptide objects from the 'matches.db' file
    def readPeptides(sequence, genome, start, end, runName, maxNum):
        lassopeptides = []

        conn = sqlite3.connect(db_dir)
        conn.create_function("REGEXP", 2, regexp)
        c = conn.cursor()
    
        # get all lasso peptides, sorted by rank
        selectionString = """SELECT * FROM lassopeptides WHERE
        start >= """ + str(start) + """ AND
        end <= """ + str(end) + """ AND
        runName LIKE '%""" + runName + """%'
        ORDER BY 5 DESC 
        LIMIT """ + str(maxNum)
        print(selectionString)
        for row in c.execute(selectionString):
            lassopeptides.append( {
                "sequence": row[0],
                "start": row[1],
                "end": row[2],
                "overallLength": row[3],
                "rank": row[4],
                "ORF": row[5],
                "genome": row[6],
                "index": row[7],
                "runName": row[8],
                "closestOrfs": json.loads(row[9]),
            })
        c.close()
        return lassopeptides

    ## Initialize Firebase access
    # Use a service account
    if(cred_file):
        cred = credentials.Certificate(cred_file)
        firebase_admin.initialize_app(cred)
    db = firestore.client()

    for peptide in readPeptides("", "", -1000, 1000000000000000, run_name, 100000):
        print("uploading " + str(peptide['sequence']))
        data = {
            "sequence": peptide['sequence'], 
            "start": peptide['start'], 
            "end": peptide['end'], 
            "overallLength": peptide['overallLength'], 
            "rank": peptide['rank'], 
            "orf": peptide['ORF'], 
            "genome": peptide['genome'], 
            "accession": peptide['index'], 
            "runName": peptide['runName'], 
            "closestOrfs": json.dumps(str(peptide['closestOrfs']))
        }
        db.collection("genomes").document(data["genome"]).set({})
        db.collection("peptides").document(data["sequence"] + str(data["start"]) + "-" + str(data["end"]) + str(data["genome"])).set(data)

def clear_firebase(cred_file):
    def delete_collection(coll_ref, batch_size):
        docs = coll_ref.limit(batch_size).stream()
        deleted = 0

        for doc in docs:
            print(f'Deleting doc {doc.id} => {doc.to_dict()}')
            doc.reference.delete()
            deleted = deleted + 1

        if deleted >= batch_size:
            return delete_collection(coll_ref, batch_size)
    
    ## Initialize Firebase access
    # Use a service account
    
    if(cred_file):
        cred = credentials.Certificate(cred_file)
        firebase_admin.initialize_app(cred)
    db = firestore.client()

    delete_collection(db.collection("genomes"), 100)
    delete_collection(db.collection("peptides"), 100000)