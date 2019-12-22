import sqlite3

conn = sqlite3.connect("matches.db")
c = conn.cursor()

conn = sqlite3.connect('matches.db')
if(True):
    c = conn.cursor()
    lassopeptides = []
    for row in c.execute("SELECT * FROM lassopeptides"):
        # print(row)
        '''
        lassopeptides.extend( {
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
        '''
        lassopeptides.extend(row)

    print(lassopeptides[0])

c.close()
