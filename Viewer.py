import sqlite3

# regular expression function for regular expression search in matches.db
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
            "closestProts": json.loads(row[9]),
            "closestProtLists": json.loads(row[10]),
        })
    c.close()
    return lassopeptides

# query for a specific peptide, all entries are optional
def specificPeptides(sequence, genome, start, end, runName, maxNum):
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
    return returnList

def getRuns():
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
    return status

def deleteRun(runName):
    os.remove("runs/" + runName + ".json");
    conn = sqlite3.connect('matches.db')
    c = conn.cursor()
    c.execute("DELETE FROM lassopeptides WHERE runName LIKE '" + runName + "%'")
    conn.commit()
    c.close()
    conn.close()

    print("Removed all entries with run name " + runName)

def cancelRun():
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

genomeList = []
conn = sqlite3.connect('matches.db')
c = conn.cursor()
for genome in c.execute("SELECT DISTINCT genome FROM lassopeptides"):
    genomeList.append(genome)
c.close()
print(genomeList)