import json
from flask import Flask
from flask import request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

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

    for peptide in lassopeptides:
        if (
            (not description or description in peptide["description"]) and 
            (not sequence or sequence in peptide["sequence"]) and 
            (not searchPattern or searchPattern in peptide["searchPattern"])
        ):
            returnList.append(peptide)

    print("query returned " + str(len(returnList)))
    return str(returnList)