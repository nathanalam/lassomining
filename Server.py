from django.conf.urls import url
from django.http import HttpResponse
import json
import os
from Downloader import downloadGenomes
from Miner import scanGenomes
import time
import sqlite3

BASE_URL = "http://e908db58.ngrok.io"

def readPeptides():
    
    lassopeptides = []

    conn = sqlite3.connect('matches.db')
    c = conn.cursor()
    for row in c.execute("SELECT * FROM lassopeptides"):
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
            "closestB": row[9],
            "closestC": row[10],
        })
    c.close()
    return lassopeptides

# Server code
DEBUG = True
SECRET_KEY = '4l0ngs3cr3tstr1ngw3lln0ts0l0ngw41tn0w1tsl0ng3n0ugh'
ROOT_URLCONF = __name__
ALLOWED_HOSTS = [BASE_URL[7:]]

def home(request):
    html = "Error finding 'index.html'"
    with open("static/index.html", 'r') as file:
        html = file.read()
        html = str(html)
        html = html.replace("$BASE_URL$", BASE_URL)
    
    return HttpResponse(html)

def give_all(request):
    return HttpResponse(str(readPeptides()), content_type="text/plain")

def get_genomes(request):
    genomeList = []
    conn = sqlite3.connect('matches.db')
    c = conn.cursor()
    for genome in c.execute("SELECT DISTINCT genome FROM lassopeptides"):
        genomeList.append(genome)
    c.close()
    return HttpResponse(str(json.dumps(genomeList)), content_type="text/plain")

def specificPeptides(request):
    # given the search parameters, return only the peptides that fit the rank
    description = request.GET.get("description")
    sequence = request.GET.get("sequence")
    searchPattern = request.GET.get("searchPattern")
    genome = request.GET.get("genome")
    minRange = request.GET.get("minRange")
    maxRange = request.GET.get("maxRange")
    runName = request.GET.get("runName")

    returnList = []

    lassopeptides = readPeptides()

    for peptide in lassopeptides:
        if (
            (not description or description in peptide["description"]) and 
            (not sequence or sequence in peptide["sequence"]) and 
            (not runName or runName in peptide["runName"]) and 
            (not searchPattern or searchPattern in peptide["searchPattern"]) and 
            (not genome or genome in peptide["genome"]) and 
            (not minRange or int(minRange) <= min(int(peptide["searchRange"][0]), int(peptide["searchRange"][1]))) and
            (not maxRange or int(maxRange) >= max(int(peptide["searchRange"][0]), int(peptide["searchRange"][1])))
        ):
            returnList.append(peptide)

    print("query returned " + str(len(returnList)))
    return HttpResponse(str(returnList), content_type='text/plain')

def launch(request):
    raw = request.GET.get("accessions")
    pattern = request.GET.get("pattern")
    runName = request.GET.get("runName")

    # start a timer
    t0 = time.time()

    # store meta data about the particular run
    runStatus = {
        "name": runName,
        "pattern": pattern,
        "input": str(raw)
    }
    phases = ["initialized", "downloaded accessions", "translated accessions", "finished mining"]

    # define a function to progressively update the current status of the run
    if not os.path.exists("runs/" + runName + ".json"):
        os.mknod("runs/" + runName + ".json")
    def updateRun(message):
        runStatus["phase"] = message
        runStatus["totalTime"] = time.time() - t0
        with open('runs/' + runName + '.json', 'w+') as outputFile:
            outputFile.write(json.dumps(runStatus))

    # download accession number genomes
    updateRun(phases[0])
    accessions = raw.split(",")
    print("downloading accessions " + str(accessions))
    downloadGenomes(accessions)
    updateRun(phases[1])

    # launch translation of nucleic acids to amino acids
    print("translating accessions")
    os.system("python3 Translator.py")
    updateRun(phases[2])

    # launch the actual mining of the translated genomes
    print("scanning genomes for lassos")
    results = scanGenomes(runName, pattern)
    updateRun(phases[3])
    print("results saved to output/" + "matches" + ".json")

    runStatus["results"] = {
        "quantity": len(results)
    }

    # record the runstatus one final time outside of helper function
    with open('runs/' + runName + '.json', 'w+') as outputFile:
        outputFile.write(json.dumps(runStatus))

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
    print("Clearing all files in genome directory")
    for dirname in ALLDIRNAMES:
        os.remove(dirname)

    return(HttpResponse("Done with run " + runName, content_type="text/plain"))

    
def status(request):
    html = "Error finding 'status.html'"
    with open("static/status.html", 'r') as file:
        html = file.read()
        html = str(html)
        html = html.replace("$BASE_URL$", BASE_URL)
    
    return HttpResponse(html)

def matches(request):
    html = "Error finding 'index.html'"
    with open("static/index.html", 'r') as file:
        html = file.read()
        html = str(html)
        html = html.replace("$BASE_URL$", BASE_URL)
    
    return HttpResponse(html)

def about(request):
    html = "Error finding 'about.html'"
    with open("static/about.html", 'r') as file:
        html = file.read()
        html = str(html)
        html = html.replace("$BASE_URL$", BASE_URL)
    
    return HttpResponse(html)

def getRuns(request):
    allRuns = []
    for dirname in os.listdir("runs"):
        with open("runs/" + dirname, 'r') as file:
            allRuns.append(json.loads(file.read()))

    print(allRuns)
    return HttpResponse(json.dumps(allRuns), content_type="text/plain")
    

urlpatterns = [
    url(r'^$', home),
    url(r'^getpeptides$', specificPeptides),
    url(r'^genomeList$', get_genomes),
    url(r'^allpeptides$', give_all),
    url(r'^launchrun$', launch),
    url(r'^status.html$', status),
    url(r'^matches.html$', matches),
    url(r'^about.html$', about),
    url(r'^getRuns$', getRuns)
]
