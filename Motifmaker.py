import os
import shutil

# Generate MEME output file from a sequence of proteins
# @ Param - protSeq is a list of tuples, with name sequence pairs
def makeMeme(protSeq, memeName):
    with open("modelProteins/" + memeName + ".faa", "w") as f:
        for seq in protSeq:
            print("writing " + str(seq))
            f.write(">" + seq[0] + "\n" + seq[1] + "\n\n")
    
    command = "/root/meme/bin/meme -nmotifs " + str(len(protSeq)) + " -maxw 25 modelProteins/" + memeName + ".faa -o motifs/" + memeName
    print(command)
    os.system(command)

    os.remove("modelProteins/" + memeName + ".faa")
    os.rename("motifs/" + memeName + "/meme.txt", "motifs/" + memeName + "Results.txt")
    shutil.rmtree("motifs/" + memeName)

# seeq = [
#     ("lmao", "AFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFA"),
#     ("asdf", "BSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSN"),
#     ("qwer", "REGREBERREGREBERREGREBERREGREBERREGREBERREGREBERREGREBERREGREBERREGREBER")
# ]

# makeMeme(seeq, "test")


