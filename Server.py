from django.conf.urls import url
from django.http import HttpResponse
import json

BASE_URL = "http://127.0.0.1:5000"

# Read the matched lasso peptides
print("Reading matches.json...")
lassopeptides = []
with open('output/matches.json', 'r') as storedfile:
    lassopeptides = json.loads(storedfile.read())
# get the set of available genomes
genomeList = {}
for peptide in lassopeptides:
    if "genome" in peptide:
        if not peptide["genome"] in genomeList:
            genomeList[peptide["genome"]] = 1
        else:
            genomeList[peptide["genome"]] += 1
print("Done reading matches.json!")


# Server code
DEBUG = True
SECRET_KEY = '4l0ngs3cr3tstr1ngw3lln0ts0l0ngw41tn0w1tsl0ng3n0ugh'
ROOT_URLCONF = __name__

def home(request):
    html = "Error finding 'index.html'"
    with open("index.html", 'r') as file:
        html = file.read()
        html = str(html)
        html = html.replace("$BASE_URL$", BASE_URL)
    
    return HttpResponse(html)  # don't use user input like that in real projects!

def give_all(request):
    return HttpResponse(str(lassopeptides), content_type="text/plain")

def get_genomes(request):
    return HttpResponse(str(json.dumps(genomeList)), content_type="text/plain")

def specificPeptides(request):
    # given the search parameters, return only the peptides that fit the rank
    description = request.GET.get("description")
    sequence = request.GET.get("sequence")
    searchPattern = request.GET.get("searchPattern")
    genome = request.GET.get("genome")
    minRange = request.GET.get("minRange")
    maxRange = request.GET.get("maxRange")

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
    return HttpResponse(str(returnList), content_type='text/plain')

urlpatterns = [
    url(r'^$', home),
    url(r'^getpeptides$', specificPeptides),
    url(r'^genomeList$', get_genomes),
    url(r'^allpeptides$', give_all),
]