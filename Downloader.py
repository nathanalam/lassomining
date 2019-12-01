'''
Author: Nathan Alam

Programmatic acquisition of genomes
'''
import os

def downloadGenomes(Accessions):

    print("Generating URL File downloads for genomes")
    for acc in Accessions:
        os.system('esearch -db assembly -query "' + acc + '" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" \'{print $0"/"$NF"_genomic.fna.gz"}\' >> fileurls.txt')

    print("Downloading files using wget into genomes folder")
    os.system("wget --directory-prefix=genomes $( cat fileurls.txt )")

    print("Unzipping downloaded files")
    os.system("gunzip -r genomes")

    os.system("rm fileurls.txt")