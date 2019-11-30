'''
Author: Nathan Alam

Programmatic acquisition of genomes
'''
import os

Accessions = [
    "NC_014816",
    "NC_014817",
    "NC_006348",
"NC_009076",
"NC_007651",
"NZ_ABBG01000168",
"NC_014722",
"NC_015381",
"NC_010335",
"NC_010338",
"NC_014100",
"NC_007777",
"NC_015580",
"NC_015579",
"NC_014006",
"NC_014007",
"NC_015976",
"NC_008048",
"NC_011071",
"NC_003155",
"NC_008278",
"NC_015574",
"NC_008609",
"NC_015957",
"NC_013525",
"NC_016582",
"NC_008702",
"NC_011969",
"NC_016024"

]

URLs = []
print("Generating URL File downloads for genomes")
for acc in Accessions:
    os.system('esearch -db assembly -query "' + acc + '" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" \'{print $0"/"$NF"_genomic.fna.gz"}\' >> fileurls.txt')

print("Downloading files using wget into genomes folder")
os.system("wget --directory-prefix=genomes $( cat fileurls.txt )")

print("Unzipping downloaded files")
os.system("gunzip -r genomes")

os.system("rm fileurls.txt")

for u in URLs:
    print(u)