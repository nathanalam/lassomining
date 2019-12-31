from django.conf.urls import url
from django.conf.urls.static import static
from django.http import HttpResponse
import json
import os
from Miner import scanGenomes, mine
import time
import sqlite3

BASE_URL = "http://50.116.48.39:8080"
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/genomemining/"

# Server code
DEBUG = True
SECRET_KEY = '4l0ngs3cr3tstr1ngw3lln0ts0l0ngw41tn0w1tsl0ng3n0ugh'
ROOT_URLCONF = __name__
ALLOWED_HOSTS = [BASE_URL[7:], BASE_URL[7:len(BASE_URL) - 5], BASE_DIR]
STATIC_URL = '/static/'
STATICFILES_DIRS = [os.path.join(BASE_DIR, 'static'),]
DEBUG=False

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
    accessions = raw.split(",")
    pattern = request.GET.get("pattern")
    runName = request.GET.get("runName")
    cutoffRank = float(request.GET.get("cutoffRank"))

    if os.path.exists("hold.txt"):
        occupant = ""
        with open("hold.txt") as file:
            occupant = file.read()
        return(HttpResponse("Occupied with run " + occupant, content_type="text/plain")) 
    
    with open("hold.txt", "w+") as file:
        file.write(runName)

    # start a timer
    t0 = time.time()

    # store meta data about the particular run
    runStatus = {
        "name": runName,
        "pattern": pattern,
        "input": str(raw),
        "progress": 0.0
    }

    # define a function to progressively update the current status of the run
    if not os.path.exists("runs/" + runName + ".json"):
        os.mknod("runs/" + runName + ".json")
        
    def updateRun(message, number):
        runStatus["phase"] = message
        runStatus["totalTime"] = time.time() - t0
        runStatus["progress"] = str((number * 1.0) / len(accessions))
        with open('runs/' + runName + '.json', 'w+') as outputFile:
            outputFile.write(json.dumps(runStatus))

    count = 0
    for accession in accessions:
        mine(accession, runName, pattern, cutoffRank)
        count += 1
        updateRun("processing" + accession, count)

    print("finished all the runs for " + runName)
    os.remove("hold.txt")

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
            for row in c.execute("SELECT rank FROM lassopeptides WHERE runName LIKE '" + particularRun["name"] + "%' ORDER BY 1 DESC"):
                ranks.append(row[0])
            c.close()
            particularRun["ranks"] = ranks
            allRuns.append(particularRun)

    return HttpResponse(json.dumps(allRuns), content_type="text/plain")

def deleteRun(request):
    runName = request.GET.get("runName")

    os.remove("runs/" + runName + ".json");
    conn = sqlite3.connect('matches.db')
    c = conn.cursor()
    # get all lasso peptides, sorted by rank

    c.execute("DELETE FROM lassopeptides WHERE runName LIKE '" + runName + "%'")

    conn.commit()
    c.close()
    conn.close()

    return HttpResponse("Removed all entries with run name " + runName, content_type="text/plain")

urlpatterns = [
    url(r'^$', home),
    url(r'^getpeptides$', specificPeptides),
    url(r'^genomeList$', get_genomes),
    url(r'^allpeptides$', give_all),
    url(r'^launchrun$', launch),
    url(r'^status.html$', status),
    url(r'^matches.html$', matches),
    url(r'^about.html$', about),
    url(r'^getRuns$', getRuns),
    url(r'^deleteRun$', deleteRun)
]
