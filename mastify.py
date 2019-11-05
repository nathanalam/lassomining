import sys
import os

ALLDIRNAMES = []
for dirname in os.listdir("genomes"):
    ## if a regular file, just add to directory
    if (dirname.find(".") != -1):
        ALLDIRNAMES.append("genomes/" + dirname)
    else:
        for filename in os.listdir("genomes/" + dirname):
            ALLDIRNAMES.append("genomes/" + dirname + "/" + filename)



os.system('export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.1.0:$PATH')

for dirname in ALLDIRNAMES:
    if (dirname[len(dirname) - 3:] == "faa"):
        outdirB = "masts/Bmatches/" + dirname[8:len(dirname) - 4] + "/"
        outdirC = "masts/Cmatches/" + dirname[8:len(dirname) - 4] + "/"
        # print(outdirB)
        # print(outdirC)
        os.makedirs(outdirB)
        os.makedirs(outdirC)
        os.system('mast -oc ' + outdirB + ' motifs/memeb.txt ' + dirname)
        os.system('mast -oc ' + outdirC + ' motifs/memec.txt ' + dirname)