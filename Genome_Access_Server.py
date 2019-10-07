import json
from flask import Flask
from flask import request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

lassopeptides = []
with open('matches.json', 'r') as storedfile:
    lassopeptides = json.loads(storedfile.read())

# get the set of available genomes
genomeList = {}
for peptide in lassopeptides:
    if not peptide["genome"] in genomeList:
        genomeList[peptide["genome"]] = 1
    else:
        genomeList[peptide["genome"]] += 1


@app.route('/allpeptides')
def give_all():
    return str(lassopeptides)

@app.route('/genomeList')
def get_genomes():
    return str(json.dumps(genomeList))

@app.route('/getpeptides', methods=["GET"])
def give_specific():
    # given the search parameters, return only the peptides that fit the rank
    description = request.args.get("description")
    sequence = request.args.get("sequence")
    searchPattern = request.args.get("searchPattern")
    genome = request.args.get("genome")
    minRange = request.args.get("minRange")
    maxRange = request.args.get("maxRange")

    returnList = []

    for peptide in lassopeptides:
        if (
            (not description or description in peptide["description"]) and 
            (not sequence or sequence in peptide["sequence"]) and 
            (not searchPattern or searchPattern in peptide["searchPattern"]) and 
            (not genome or genome in peptide["genome"]) and 
            (not minRange or int(minRange) <= min(int(peptide["searchRange"][0]), int(peptide["searchRange"][1]))) and
            (not maxRange or int(maxRange) >= max(int(peptide["searchRange"][0]), int(peptide["searchRange"][1])))
        ):
            returnList.append(peptide)

    print("query returned " + str(len(returnList)))
    return str(returnList)