import sqlite3 as sqlite3
import json
import re

def regexp(expr, item):
    reg = re.compile(expr)
    return reg.search(item) is not None

sequence = "YG.*YG"
conn = sqlite3.connect('matches.db')
conn.create_function("REGEXP", 2, regexp)
c = conn.cursor()

lassopeptides = []
selectionString = """SELECT * FROM lassopeptides WHERE
    sequence REGEXP '""" + sequence + """'
    ORDER BY 5 DESC 
    LIMIT 10 """

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
print(lassopeptides)