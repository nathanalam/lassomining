import os
import shutil

# Generate MEME output file from a list of protein sequences
# @ Param - motifJobs is a list of triplets, with directory - width - number of motif groupings
def makeMemes(motifJobs, memeDir, motifDir):

    for model in os.listdir(modelDir):
        memeName = model[0:len(model) - 4]
        command = memeDir + "/bin/meme -nmotifs " + str(nmotifs) + " -maxw 25" + modelDir + model + "-o " + motifDir + memeName
        print(command)
        os.system(command)

        os.rename(motifDir + memeName + "/meme.txt", motifDir + memeName + "Results.txt")
        shutil.rmtree(motifDir + memeName)

# seeq = [
#     ("lmao", "AFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFAAFAFAFAFAAFA"),
#     ("asdf", "BSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSNBSNSNSNSN"),
#     ("qwer", "REGREBERREGREBERREGREBERREGREBERREGREBERREGREBERREGREBERREGREBERREGREBER")
# ]

# makeMeme(seeq, "test")


