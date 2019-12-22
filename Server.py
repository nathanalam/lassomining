from django.conf.urls import url
from django.http import HttpResponse
import json
import os
from Downloader import downloadGenomes
from Miner import scanGenomes
import time
import sqlite3

BASE_URL = "http://bca0fa6d.ngrok.io"

# Server code
DEBUG = True
SECRET_KEY = '4l0ngs3cr3tstr1ngw3lln0ts0l0ngw41tn0w1tsl0ng3n0ugh'
ROOT_URLCONF = __name__
ALLOWED_HOSTS = [BASE_URL[7:]]

def readPeptides(sequence, genome, start, end, runName, maxNum):
    
    lassopeptides = []

    conn = sqlite3.connect('matches.db')
    c = conn.cursor()
    # get all lasso peptides, sorted by rank
    selectionString = """SELECT * FROM lassopeptides WHERE
    sequence LIKE '%""" + sequence + """%' AND
    genome LIKE '%""" + genome + """%' AND
    start >= """ + str(start) + """ AND
    end <= """ + str(end) + """ AND
    runName LIKE '%""" + runName + """%'
    ORDER BY 5 ASC 
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
            "closestBs": json.loads(row[9]),
            "closestCs": json.loads(row[10]),
        })
    c.close()
    return lassopeptides

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
    sequence = request.GET.get("sequence")
    genome = request.GET.get("genome")
    start = request.GET.get("minRange")
    end = request.GET.get("maxRange")
    runName = request.GET.get("runName")
    maxNum = request.GET.get("maxNum")

    if not start:
        start = 0

    if not end:
        end = 1000000000000000

    returnList = readPeptides(sequence, genome, start, end, runName, maxNum)

    print("query returned " + str(len(returnList)))
    return HttpResponse(json.dumps(returnList), content_type='text/json')

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
            particularRun = json.loads(file.read())
            ranks = []

            conn = sqlite3.connect('matches.db')
            c = conn.cursor()
            # get all lasso peptides, sorted by rank
            for row in c.execute("SELECT rank FROM lassopeptides WHERE runName LIKE '" + particularRun["name"] + "%' ORDER BY 1 ASC"):
                ranks.append(row[0])
            c.close()
            particularRun["ranks"] = ranks
            allRuns.append(particularRun)

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
