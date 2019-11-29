'''
Author: Nathan Alam

Programmatic acquisition of genomes
'''

Accessions = [
    "NC_014816",
    "NC_014817",
    "NC_006348",
    "NC_014816"
]

URLs = []
for acc in Accessions:
    URLs.append("https://www.ncbi.nlm.nih.gov/nuccore/" + acc + "?report=fasta&log$=seqview&format=text")

for u in URLs:
    print(u)