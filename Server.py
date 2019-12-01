from django.conf.urls import url
from django.http import HttpResponse
import json
import os
from Downloader import downloadGenomes
import time

BASE_URL = "http://127.0.0.1:5000"


# Server code
DEBUG = True
SECRET_KEY = '4l0ngs3cr3tstr1ngw3lln0ts0l0ngw41tn0w1tsl0ng3n0ugh'
ROOT_URLCONF = __name__

def home(request):
    html = "Error finding 'index.html'"
    with open("static/index.html", 'r') as file:
        html = file.read()
        html = str(html)
        html = html.replace("$BASE_URL$", BASE_URL)
    
    return HttpResponse(html)  # don't use user input like that in real projects!

def give_all(request):
    runName = request.GET.get("runName")
    # Read the matched lasso peptides
    print("Reading output/" + runName + ".json...")
    lassopeptides = []
    with open('output/' + runName + '.json', 'r') as storedfile:
        lassopeptides = json.loads(storedfile.read())
    print("Done reading " + runName + ".json!")
    return HttpResponse(str(lassopeptides), content_type="text/plain")

def get_genomes(request):
    runName = request.GET.get("runName")
    # get the set of available genomes
    genomeList = {}
    print("Reading output/" + runName + ".json...")
    lassopeptides = []
    with open('output/' + runName + '.json', 'r') as storedfile:
        lassopeptides = json.loads(storedfile.read())
    print("Done reading " + runName + ".json!")
    for peptide in lassopeptides:
        if "genome" in peptide:
            if not peptide["genome"] in genomeList:
                genomeList[peptide["genome"]] = 1
            else:
                genomeList[peptide["genome"]] += 1
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

    print("Reading output/" + runName + ".json...")
    lassopeptides = []
    with open('output/' + runName + '.json', 'r') as storedfile:
        lassopeptides = json.loads(storedfile.read())
    print("Done reading " + runName + ".json!")

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
    return HttpResponse(str(returnList), content_type='text/plain')

def launch(request):
    raw = request.GET.get("accessions")
    pattern = request.GET.get("pattern")
    runName = request.GET.get("runName")

    t0 = time.time()
    runStatus = {
        "name": runName,
        "pattern": pattern,
        "input": str(raw)
    }

    def updateRun(message):
        runStatus["phase"] = message
        runStatus["totalTime"] = time.time() - t0
        with open('runs/' + runName + '.json', 'w') as outputFile:
            outputFile.write(json.dumps(runStatus))

    updateRun("initialized")
    accessions = raw.split(",")
    print("downloading accessions " + str(accessions))
    downloadGenomes(accessions)
    updateRun("downloaded accessions")

    print("translating accessions")
    os.system("python3 Translator.py")
    updateRun("translated accessions")

    print("scanning genomes for lassos")
    scanGenomes(runName, pattern)
    updateRun("finished mining")
    print("results saved to output/" + runName + ".json")

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
    
urlpatterns = [
    url(r'^$', home),
    url(r'^getpeptides$', specificPeptides),
    url(r'^genomeList$', get_genomes),
    url(r'^allpeptides$', give_all),
    url(r'^launchrun$', launch)
]