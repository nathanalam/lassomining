from django.conf.urls import url
from django.conf.urls.static import static
from django.http import HttpResponse
import json
import os
from Miner import scanGenomes, mine
import time
import sqlite3
import re

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

# regular expression function for regular expression search
def regexp(expr, item):
    reg = re.compile(expr)
    return reg.search(item) is not None

# read peptide objects from the 'matches.db' file
def readPeptides(sequence, genome, start, end, runName, maxNum):
    lassopeptides = []

    conn = sqlite3.connect('matches.db')
    conn.create_function("REGEXP", 2, regexp)
    c = conn.cursor()
    
    # get all lasso peptides, sorted by rank
    selectionString = """SELECT * FROM lassopeptides WHERE
    sequence REGEXP '""" + sequence + """' AND
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

    if (sequence is None) or (len(sequence) is 0):
        sequence = ".*"
    if (maxNum is None) or (len(maxNum) is 0):
        maxNum = 10
    returnList = readPeptides(sequence, genome, start, end, runName, maxNum)

    print("query returned " + str(len(returnList)))
    return HttpResponse(json.dumps(returnList), content_type='text/json')

def launch(request):
    
    try:
        runName = request.POST.get("runName")
        pattern = request.POST.get("pattern")     
        cutoffRank = float(request.POST.get("cutoffRank"))
        raw = request.POST.get("accessions")
        print(len(raw))
        accessions = raw.split(",")
        
    except Exception as e:
        print("Could not read responses due to " + str(e))
        return(HttpResponse("Improper POST input", content_type="text/plain"))

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
        "input": [],
        "progress": 0.0,
        "peptides": 0,
        "cutoff": cutoffRank
    }

    # define a function to progressively update the current status of the run
    if not os.path.exists("runs/" + runName + ".json"):
        os.mknod("runs/" + runName + ".json")
 
    def updateRun(message, number, count, accession):
        runStatus["phase"] = message
        runStatus["totalTime"] = time.time() - t0
        runStatus["input"].append(accession)
        runStatus["progress"] = str((number * 1.0) / len(accessions))
        runStatus["peptides"] = count
        with open('runs/' + runName + '.json', 'w+') as outputFile:
            outputFile.write(json.dumps(runStatus))

    count = 0
    peptideCount = 0
    returnText = "Done with run " + runName

    try:
        # combine accessions with the accessions in the accessions buffer text file
        with open("accessions.txt", "a+") as file:
            for accession in accessions:
                file.write(accession + "\n")
        
        f = open('accessions.txt')
        accession = f.readline()
        while accession:
            peptideCount += mine(accession, runName, pattern, cutoffRank)
            count += 1
            updateRun("processing" + accession, count, peptideCount, accession)
            accession = f.readline()
        f.close()
            

        print("finished all the runs for " + runName)
    except Exception as e: 
        print("error occured while doing run " + runName)
        returnText = "Failed run " + runName + " with error " + str(e)
    
    if os.path.exists("accessions.txt"):
        os.remove("accessions.txt")
    if os.path.exists("hold.txt"):
        os.remove("hold.txt")

    return(HttpResponse(returnText, content_type="text/plain"))

    
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
        try: 
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
        except:
            # return HttpResponse("Failed to read " + dirname, content_type="text/plain")
            print("Failed to read " + dirname)
    currentRun = ""
    busy = False
    if os.path.exists("hold.txt"):
        busy = True
        with open("hold.txt") as file:
            currentRun = file.read()

    status = {
        "meta": {
            "busy": busy,
            "currentRun": currentRun
        },
        "allRuns": allRuns
    }
    return HttpResponse(json.dumps(status), content_type="text/plain")

def deleteRun(request):
    runName = request.GET.get("runName")

    os.remove("runs/" + runName + ".json");
    conn = sqlite3.connect('matches.db')
    c = conn.cursor()
    c.execute("DELETE FROM lassopeptides WHERE runName LIKE '" + runName + "%'")
    conn.commit()
    c.close()
    conn.close()

    return HttpResponse("Removed all entries with run name " + runName, content_type="text/plain")

def cancelRun(request):
    try: 
        if os.path.exists("accessions.txt"):
            os.remove("accessions.txt")
        if os.path.exists("hold.txt"):
            os.remove("hold.txt")
        
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
    except:
        return HttpResponse("Could not cancel run", content_type="text/plain") 

    return HttpResponse("Cancelled the current run", content_type="text/plain")

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
    url(r'^deleteRun$', deleteRun),
    url(r'^cancelRun$', cancelRun)
]
