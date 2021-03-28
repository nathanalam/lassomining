import os
import time
import json
import shutil
import sys
import yaml
from xml.dom import minidom
import re
import traceback
import sys
import math
import pandas as pd
import sqlite3
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import math
from pathlib import Path
from Bio.Seq import Seq, reverse_complement, translate
import subprocess

import firebase_admin
from firebase_admin import credentials
from firebase_admin import firestore

import dill as pickle
import numpy as np

with open(os.path.join(os.path.dirname(__file__), 'NN.pickle'), 'rb') as g:
    clf = pickle.load(g)

# some flags for debugging
REMOVE_GENOMES_ON_TRANSLATE = True
PRINT_EACH_FIND = False
PRINT_EACH_WRITE = False
MAX_MOTIF_NUMS = [3,3,4]
# put -1 to take all above the cutoff per ORF
TAKE_TOP_N = -1
SECONDARY_RANK_CUTOFF = 0
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

    # print("Read " + str(count + 1) + " objects from FASTA file " + name)

    return sequenceList


'''
Given a protein sequence, meme motifs, and a meme directory, return a list of open reading
frames in that sequence with associated motifs
'''


def mast_orfs(sequence, motifs, memeInstall, readingFrame, filenam):

    orfs = []
    m_stack = []
    index = 0
    for amino_acid in sequence:
        if (amino_acid == 'M'):
            m_stack.append(index)
        elif (amino_acid == '*'):
            for m in m_stack:
                orfs.extend([{
                    "start": m,
                    "end": index,
                    "counts": {},
                    "motifs": {}
                } for m in m_stack])
            m_stack = []
        index += 1
    with open(f"{filenam[:len(filenam)-1]}.tempseq.txt", "w") as file:
        file.write("> " + "temporary" + "\n")
        file.write(sequence)
        file.close()

    try:
        for motif in motifs:
            (_, name) = os.path.split(motif)
            # command = memeInstall + '/bin/mast -nostatus -hit_list ' + motif + f' tempseq{filenam.replace("/", "-")}.txt > tempout{filenam.replace("/", "-")}' + name
            # command = memeInstall + '/bin/mast -hit_list ' + motif + f' tempseq{filenam.replace("/", "-")}.txt'
            # -minseqs 3 -remcorr -ev 10.0 -dl http://www.uniprot.org/uniprot/?query=SEQUENCEID&sort=score
            # command = [
            #     memeInstall + '/bin/mast', '-hit_list', motif,
            #     f"{filenam[:len(filenam)-1]}.tempseq.txt"
            # ]
            command = [
                memeInstall + '/bin/mast', '-hit_list', '-minseqs', '3', '-remcorr', '-ev', '10', '-dl', 'http://www.uniprot.org/uniprot/?query=SEQUENCEID&sort=score', motif,
                f"{filenam[:len(filenam)-1]}.tempseq.txt"
            ]
            # print(command)

            matched_motifs = []
            output = subprocess.check_output(command,
                                             stderr=subprocess.DEVNULL)

            inlines = output.decode("utf-8").split('\n')
            inlines = inlines[2:len(inlines) - 2]

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
                    if (newProt["p-value"] < 0.5):
                        matched_motifs.append(newProt)
                except:
                    print("error in parsing line - " + line)
                    print("params: " + str(params))

            # assign prots to orfs
            for prot in matched_motifs:
                for orf in orfs:
                    # add if overlap and higher score than any existing motif match
                    if (orf['start'] <= prot['start']
                            and orf['end'] >= prot['end']):
                        already_exists = False
                        if hasattr(orf['motifs'], prot["memeDir"]):
                            for current_index in range(
                                    0, len(orf['motifs'][prot["memeDir"]])):
                                existing_prot = orf['motifs'][
                                    prot["memeDir"]][current_index]
                                if existing_prot['motif'] == prot['motif']:
                                    already_exists = True
                                    # swap if the new prot has higher score
                                    if prot["score"] > existing_prot["score"]:
                                        orf['motifs'][prot["memeDir"]][
                                            current_index] = prot
                        if (not already_exists):
                            if (prot["memeDir"] in orf['counts']):
                                orf['counts'][prot["memeDir"]] += 1
                            else:
                                orf['counts'][prot["memeDir"]] = 1
                            if (prot["memeDir"] in orf['motifs']):
                                orf['motifs'][prot["memeDir"]].append(prot)
                            else:
                                orf['motifs'][prot["memeDir"]] = [prot]
    except Exception as e:
        if(not "<Signals.SIG" in str(e)):
            print("An error occured with running MAST")
            print(e)
        if (os.path.exists(f"{filenam[:len(filenam)-1]}.tempseq.txt")):
            os.remove(f"{filenam[:len(filenam)-1]}.tempseq.txt")

    if (os.path.exists(f"{filenam[:len(filenam)-1]}.tempseq.txt")):
        os.remove(f"{filenam[:len(filenam)-1]}.tempseq.txt")
    # used 3 motifs for the b motif, 4 motifs for the c motif
    max_motif_nums = [3, 4]
    matched_orfs = []

    max_motif_num = []
    for motif_index, motif in enumerate(motifs):
        (_, name) = os.path.split(motif)
        these_orfs = []
        max_motif_num = MAX_MOTIF_NUMS[motif_index]
        t = 0
        for orf in orfs:
            if (orf["start"] <= t):
                # skip if viewing duplicate orf
                continue
            t = orf["start"]
            if (motif in orf['counts']):
                these_orfs.append({
                    'start': orf['start'],
                    'end': orf['end'],
                    'count': orf['counts'][motif] / MAX_MOTIF_NUMS,
                    'motifs': orf['motifs'][motif],
                    'readingFrame': readingFrame,
                    'motifType': name,
                    'sequence': sequence[orf['start']:orf['end']]
                })
        matched_orfs.append(these_orfs)
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
        motifMatches = mast_orfs(overallSequence, motifs, memeInstall, RFindex,
                                 filenam)

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

                rank = 0

                closestOrfs = []
                """
                Rank the peptide according to the identified motifs
                """
                # go through each match set, I ~= B, C, D, etc.
                for IProteins in AuxProteins:
                    closestIs = []
                    center = (start + end) / 2
                    for IProtein in IProteins:
                        Icenter = (IProtein['start'] + IProtein['end']) / 2
                        if (abs(center - Icenter) < 10000
                                and not isOverlapping(start, end,
                                                      IProtein['start'],
                                                      IProtein['end'])):
                            closestIs.append(IProtein)

                    def sortFunct(prot):
                        return prot["count"] / (prot["start"] - center)

                    if (len(closestIs) > 0):
                        closestIs.sort(key=sortFunct, reverse=True)
                        rank += closestIs[0]['count']
                        closestOrfs.append(closestIs[0])

                rank = rank / len(AuxProteins)

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


def vectorize(sequence_list):
    vector_arr = []
    aa_rep = np.matrix([[
        0.20412415, 0.15567426, 0.03942035, -0.01822488, -0.14819245,
        0.19802112, -0.43177231, 0.3272457, -0.18397928, 0.37255573,
        -0.10809339, 0.2117906, -0.06704938, 0.00071176, -0.22053608,
        0.35334939, -0.0543725, 0.17522578, -0.06132986, 0.22112548,
        -0.18809861, -0.03007737, 0.21668077, 0.02817936
    ],
                        [
                            0.20412415, -0.16717617, 0.04522636, 0.0527608,
                            -0.06514745, 0.10555747, 0.3294438, 0.15377643,
                            -0.14839629, 0.05564137, 0.44142262, 0.26383881,
                            -0.27165656, -0.21415267, -0.03675024, -0.13039611,
                            0.30409053, 0.06596975, -0.37842362, 0.18876561,
                            0.22838388, -0.06952343, -0.03159155, -0.12939871
                        ],
                        [
                            0.20412415, 0.1718496, 0.21887519, 0.1676149,
                            0.0961377, 0.15366935, -0.12880583, 0.16764423,
                            -0.43199321, -0.04598926, -0.06827292, -0.35602773,
                            0.14094867, 0.10767668, 0.18162269, -0.3756231,
                            0.24412855, 0.09623603, -0.10614588, -0.22846791,
                            -0.2298419, 0.02849644, -0.03571025, -0.24963588
                        ]])
    btranslator = np.matrix([[
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0
    ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0
                             ],
                             [
                                 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
                             ],
                             [
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                             ]])

    def get_cores(sequence_list):
        leaders = []
        cores = []
        tails = []
        for seq in sequence_list:
            last_ind = 0
            while ((last_ind < len(seq) - 1) and last_ind != -1):
                try:
                    new_term = seq.index('T', last_ind + 1)
                    if (new_term > 15 and new_term <= 46):
                        leader = seq[0:new_term + 2]
                        remainder = seq[new_term + 2:]
                        core = re.findall("[A-Z]{6,8}[DE]", remainder)[0]
                        remainder = remainder[len(core):]
                        #M[A-Z]{15,45}T[A-Z][A-Z]{6,8}[DE][A-Z]{5,30}\*
                        leaders.append(leader)
                        cores.append(core)
                        tails.append(remainder)
                    last_ind = new_term
                except:
                    break
        return {"leaders": leaders, "cores": cores, "tails": tails}

    def vectorize_swanson(sequence_list):
        maxLen = 9
        arr = np.zeros((26 * maxLen, len(sequence_list)))
        for i in range(0, len(sequence_list)):
            rowIndex = 0
            for char in sequence_list[i]:
                arr[rowIndex + (ord(char) - 65), i] = 1
                rowIndex += 26
            mod_arr = np.zeros((3 * maxLen, len(sequence_list)))
            rowIndex = 0
            for i in range(0, np.shape(arr)[0], 26):
                new_three = np.matmul(
                    aa_rep, np.matmul(btranslator, arr[i:(i + 26), :]))
                mod_arr[rowIndex:(rowIndex + 3), :] = new_three
                rowIndex += 3
        return mod_arr

    s_list = get_cores(sequence_list)["cores"]

    return vectorize_swanson(s_list)


def classify(sequence_list):
    scores = []
    for seq in sequence_list:
        vector_matrix = np.matrix(vectorize([seq]))
        max_score = 0
        for i in range(0, vector_matrix.shape[1]):
            temp_score = np.matrix(
                clf.predict_proba(np.transpose(vector_matrix[:, i])))[:, 1]
            if (temp_score > max_score):
                max_score = temp_score
        scores.append(np.sum(max_score))
        # scores.append(np.maximum(np.matrix(clf.predict_proba(vectorize([seq])))[:,1]))
    return (scores)


def secondary_rank(peptide):

    seq = peptide["sequence"].replace('*', '')
    modifier = np.sum(classify([seq]))
    # multiply modifier by inverse distance rank
    for accessory_orf in peptide["closestOrfs"]:
        dist = ((peptide["searchRange"][0] - accessory_orf["start"])**2)
        modifier = modifier / np.log10(dist)
    return (modifier * peptide["rank"])


def scanGenome(runName, pattern, cutoffRank, databaseDir, memeInstall,
               genomeDir, motifs):
    conn = sqlite3.connect(databaseDir)

    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS lassopeptides
             (sequence text, start integer, end integer, overallLength integer, rank real, orf integer, genome text, accession text, runName text, closestOrfs text, secondaryRank real)'''
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
        # print("Found " + str(len(buffer)) + " peptides in this set of ORFs")
        for peptide in buffer:
            if (PRINT_EACH_WRITE):
                print("Inserting " + peptide['sequence'] +
                      " into sqlite database")
            adjusted_rank = secondary_rank(peptide)
            if (adjusted_rank < SECONDARY_RANK_CUTOFF):
                continue
            submit_buffer = False
            while (not submit_buffer):
                try:
                    c.execute(
                        "INSERT INTO lassopeptides VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        [
                            peptide['sequence'], peptide['searchRange'][0],
                            peptide['searchRange'][1],
                            peptide['overallLength'], peptide["rank"],
                            peptide['readingFrame'], peptide['genome'],
                            peptide['index'], peptide['runName'],
                            json.dumps(str(
                                peptide['closestOrfs'])), adjusted_rank
                        ])
                    submit_buffer = True
                except:
                    print(databaseDir + " is busy, waiting 5 seconds")
                    time.sleep(5)

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


def mine_process(dirname, ALLDIRNAMES, genomeFolder, runName, pattern,
                 cutoffRank, databaseDir, memeInstall, motifs):
    count = 0
    translatedDirectory = ""
    if (((dirname[len(dirname) - 3:] == "fna") or
         (dirname[len(dirname) - 5:] == "fasta"))
            and not os.path.exists(os.path.join(genomeFolder, (dirname[:len(dirname) - 3] + "faa")))):
        # print("Opening up " + dirname +
        #      " and converting into peptide sequences...")
        DNAseqs = []
        seqDescriptions = []
        try:
            for fastaobj in readFASTA(genomeFolder + dirname):
                DNAseqs.append(fastaobj["sequence"])
                seqDescriptions.append(fastaobj["description"])
        except:
            print(f"Error in reading FASTA for {genomeFolder + dirname}")
            return

        if (REMOVE_GENOMES_ON_TRANSLATE):
            try:
                os.remove(genomeFolder + dirname)
            except:
                print(f"Could not remove {genomeFolder + dirname}")
                return

        entries = []
        for i in range(0, len(DNAseqs)):
            # print("converting " + str(len(DNAseqs[i])) + " base pairs from " +
            #      seqDescriptions[i])
            aalist = get_reading_frames(DNAseqs[i])
            # print("created " + str(len(aalist)) + " peptide sequences from " +
            #      seqDescriptions[i])
            for e in range(0, len(aalist)):
                entries.append({
                    "sequence":
                    aalist[e]["sequence"],
                    "description":
                    str(seqDescriptions[i] + " - ORF " + str(aalist[e]["ORF"]))
                })
        suffixNum = 5
        if (dirname[len(dirname) - 3:] == "fna"):
            suffixNum = 3

        translatedDirectory = genomeFolder + dirname[:len(dirname) -
                                                     suffixNum] + "faa"

        # print("writing read peptides into '" + translatedDirectory + "'")
        try:
            with open(translatedDirectory, 'w') as outfile:
                for ent in entries:
                    try:
                        outfile.write("> " + ent["description"] + "\n")
                        outfile.write(ent["sequence"] + "\n\n")
                    except OSError as e:
                        print(e)
        except OSError as e:
            print(e)
    else:
        return
    # launch the actual mining of the translated genomes
    # print("scanning " + dirname + " for lassos")
    try:
        results = scanGenome(runName, pattern, cutoffRank, databaseDir,
                             memeInstall, translatedDirectory, motifs)
        count += len(results)
        # print("found " + str(count) + " peptides")
    except Exception as e:
        print("An error occured while mining " + dirname)
        print(e)

    ## clear the genomes subdirectory
    # print("removing " + translatedDirectory)
    if (os.path.exists(translatedDirectory)):
        os.remove(translatedDirectory)


def mine(genomeFolder, runName, pattern, cutoffRank, databaseDir, memeInstall,
         motifs):

    ## translate the downloaded file into amino acid sequence
    count = 0

    print("translating fna files in directory folder " + genomeFolder)
    ALLDIRNAMES = os.listdir(genomeFolder)

    for dirname in ALLDIRNAMES:
        mine_process(dirname, ALLDIRNAMES, genomeFolder, runName, pattern,
                  cutoffRank, databaseDir, memeInstall, motifs)

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

    selectionStringGenomes = "SELECT DISTINCT sequence, rank, genome, start, end, accession, closestOrfs, secondaryRank FROM lassopeptides WHERE runname is '" + run_name + "'"
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
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
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
            lassopeptides.append({
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
                "secondaryRank": row[10],
            })
        c.close()
        return lassopeptides

    ## Initialize Firebase access
    # Use a service account
    if (cred_file):
        cred = credentials.Certificate(cred_file)
        firebase_admin.initialize_app(cred)
    db = firestore.client()

    for peptide in readPeptides("", "", -1000, 1000000000000000, run_name,
                                100000):
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
            "closestOrfs": json.dumps(str(peptide['closestOrfs'])),
            "secondaryRank": peptide['secondaryRank']
        }
        db.collection("genomes").document(data["genome"]).set({})
        db.collection("peptides").document(data["sequence"] +
                                           str(data["start"]) + "-" +
                                           str(data["end"]) +
                                           str(data["genome"])).set(data)


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

    if (cred_file):
        cred = credentials.Certificate(cred_file)
        firebase_admin.initialize_app(cred)
    db = firestore.client()

    delete_collection(db.collection("genomes"), 100)
    delete_collection(db.collection("peptides"), 100000)
