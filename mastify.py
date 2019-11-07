import sys
import os

DIRNAMES = []
for dirname in os.listdir("genomes"):
    if (dirname.find(".") != -1):
        if(dirname[len(dirname) - 3:] == "faa"):
            DIRNAMES.append("genomes/" + dirname)
    else:
        for filename in os.listdir("genomes/" + dirname):
            if(filename[len(filename) - 3:] == "faa"):
                DIRNAMES.append("genomes/" + dirname + "/" + filename)
    
        
print(DIRNAMES)

def readFASTA(name, cleanspace = 0):
    descriptions = []
    sequences = []
    sequenceList = []
    tempSequences = []     
        
    with open(name) as file:
        count = -1
        for line in file:
            
            if(line[0] == '>'):
                # if begins with a >, then a description
                descriptions.append(line[1:].replace('\n', ''))
                count += 1
                # skip the first time
                if count > 0 :
                    # combine the tempSequences into a single string and
                    # add it to sequences
                    newSequence = ' '.join(tempSequences)
                    # now remove all of the whitespaces
                    newSequence = newSequence.replace(' ', '')
                    newSequence = newSequence.replace('\n', '')
                    
                    sequences.append(newSequence)
                    # refresh the tempSequence list
                    tempSequences = []
                    
                    sequenceList.append({
                        "description": descriptions[count - 1],
                        "sequence": sequences[count - 1]
                    })
            else:
                tempSequences.append(line)
                
        # combine the tempSequences into a single string and
        # add it to sequences
        newSequence = ' '.join(tempSequences)
        # now remove all of the whitespaces
        newSequence = newSequence.replace(' ', '')
        newSequence = newSequence.replace('\n', '')

        sequences.append(newSequence)
        # refresh the tempSequence list
        tempSequences = []
        
        sequenceList.append({
            "description": descriptions[count],
            "sequence": sequences[count]
        })
                
                
    if len(descriptions) != len(sequences):
        print("ERROR: Number of descriptions does not match number of sequences")
        print("Number of descriptions: " + str(len(descriptions)))
        print("Number of sequences: " + str(len(sequences)))
        sys.exit(1);
        
    print("Read " + str(count + 1) + " objects from FASTA file " + name)
        
    return sequenceList

os.system('export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.1.0:$PATH')

for fastafile in DIRNAMES:
    for pair in readFASTA(fastafile):
        desc = pair["description"]
        seq = pair["sequence"]

        with open("temp.txt", "w") as file:
            file.write("> " + desc + "\n")
            file.write(seq)

        os.system('mast -hit_list motifs/memeb.txt temp.txt')

        os.remove("temp.txt")



# for dirname in ALLDIRNAMES:
#     if (dirname[len(dirname) - 3:] == "faa"):
#         outdirB = "masts/Bmatches/" + dirname[8:len(dirname) - 4] + "/"
#         outdirC = "masts/Cmatches/" + dirname[8:len(dirname) - 4] + "/"
#         # print(outdirB)
#         # print(outdirC)
#         os.makedirs(outdirB)
#         os.makedirs(outdirC)
#         os.system('mast -oc ' + outdirB + ' motifs/memeb.txt ' + dirname)
#         os.system('mast -oc ' + outdirC + ' motifs/memec.txt ' + dirname)