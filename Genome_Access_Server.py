import json
from flask import Flask
from flask import request
app = Flask(__name__)

lassopeptides = []
with open('matches.json', 'r') as storedfile:
    lassopeptides = json.loads(storedfile.read())

@app.route('/allpeptides')
def give_all():
    return str(lassopeptides)

@app.route('/getpeptides', methods=["GET"])
def give_specific():
    # given the search parameters, return only the peptides that fit the rank
    description = request.args.get("description")
    sequence = request.args.get("sequence")
    searchPattern = request.args.get("searchPattern")

    returnList = []

    print(str(description) + " - " + str(sequence) + " - " + str(searchPattern))
    print(not description)
    print(not sequence)
    print(not searchPattern)

    for peptide in lassopeptides:
        if (
            (not description or description in peptide["description"]) and 
            (not sequence or sequence in peptide["sequence"]) and 
            (not searchPattern or searchPattern in peptide["searchPattern"])
        ):
            returnList.append(peptide)

    print("query returned " + str(len(returnList)))
    return str(returnList)