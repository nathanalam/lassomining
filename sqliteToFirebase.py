import os
import re
import sqlite3
import json

import firebase_admin
from firebase_admin import credentials
from firebase_admin import firestore

# regular expression function for regular expression search
def regexp(expr, item):
    reg = re.compile(expr)
    return reg.search(item) is not None

# read peptide objects from the 'matches.db' file
def readPeptides(sequence, genome, start, end, runName, maxNum):
    lassopeptides = []

    conn = sqlite3.connect('/home/blucheez/Projects/lassomining/output/matches.db')
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
            "closestProts": json.loads(row[9]),
            "closestProtLists": json.loads(row[10]),
        })
    c.close()
    return lassopeptides

## Initialize Firebase access
# Use a service account
cred = credentials.Certificate('/home/blucheez/Projects/lassomining/output/lasso-peptides-51ce2e6250b9.json')
firebase_admin.initialize_app(cred)
    
db = firestore.client()

for peptide in readPeptides("", "", -1000, 1000000000000000, "referencerun", 100000):
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
        "closestProts": json.dumps(str(peptide['closestProts'])), 
        "closestProtLists": json.dumps(str(peptide['closestProtLists']))
    }
    db.collection("genomes").document(data["genome"]).set({})
    db.collection("peptides").document(data["sequence"] + str(data["start"]) + "-" + str(data["end"]) + str(data["genome"])).set(data)
