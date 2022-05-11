#Checks the designatiions for a given set of lineages in an UShER metadata file
#For each designated sequence, checks whether its designated lineage matches its pango_lineage_usher
#If not, writes the sequence to 2 files - 1 containing only sequence names, the other sequence names and details
#Also writes a summary of the designated lineages that don't match
#Only outputs sequences that are in the UShER tree but not in the lineage clade

import argparse
from ssl import VerifyFlags

#Checks headers in the UShER metadata file are in the expected positions
def verifyUShER(uFile):
    with open(uFile) as f:
        h = f.readline()
        uC = 0
        for i, c in enumerate(h.strip().split("\t")):
            if c == "pango_lineage_usher":
                uC = i
        
        if uC != 10:
            raise RuntimeError("pango_lineage_usher column is not in the expected location in the metadata file")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", help = "UShER metadata file containing UShER assigned lineages")
    parser.add_argument("-d", help = "lineages.csv containing current designations")
    parser.add_argument("-l", help = "File containing lineages to be checked, 1 per line with no header")
    parser.add_argument("-o", help = "Prefix of output files, deafult = usher_check", default = "usher_check")
    args = parser.parse_args()

    #Verify UShER headers
    verifyUShER(args.m)

    #Import the lineages and extract to set
    lS = set()
    with open(args.l) as f:
        for l in f:
            lS.add(l.strip())
    
    #Iterate through the UShER metadata file, check if the pango_lineage_usher is in lS, add to
    #lineage dictionary if so
    lD = dict()
    with open(args.m) as f:
        for l in f:
            lineage = l.strip().split("\t")[10]
            if lineage in lS:
                if lineage not in lD:
                    lD[lineage] = set()
                lD[lineage].add(l.strip().split("|")[0])
    
    #Check if any provided lineages are not in the metadata file
    toExclude = list()
    for i in lS:
        if i not in lD.keys():
            toExclude.append(i)
            print("Lineage " + i + " not in metadata file so was not examined")
    if len(toExclude) != 0:
        for i in toExclude:
            lS.remove(i)
    
    outNames = open(args.o + ".txt", "w")
    outDesc = open(args.o + "_lineages.csv", "w")
    outDesc.write("Sequence_name,Designated_lineage\n")

    #Lineage counts that don't match
    lC = dict()

    #Sequences whose designation doesn't match their UShER lineage
    mmS = set()
    mmD = dict()

    #Iterate through the designated lineages, check if their lineage is in the set of interest, if so
    #check if they are in the UShER assignments to that lineage, if not write
    with open(args.d) as f:
        for l in f:
            lineage = l.strip().split(",")[1]
            if lineage in lS:
                sN = l.strip().split(",")[0]
                if sN not in lD[lineage]:
                    mmS.add(sN)
                    mmD[sN] = lineage
                    #outNames.write(sN + "\n")
                    #outDesc.write(sN + "," + lineage + "\n")
                    #if lineage not in lC:
                    #    lC[lineage] = 0
                    #lC[lineage] += 1
    
    #Iterate through the UShER sequences, check if they are in the mismatch sequences and write if so
    with open(args.m) as f:
        for l in f:
            sN = l.strip().split("|")[0] 
            if sN in mmS:
                outNames.write(sN + "\n")
                outDesc.write(sN + "," + mmD[sN] + "\n")

                if mmD[sN] not in lC:
                    lC[mmD[sN]] = 0
                lC[mmD[sN]] += 1
    
    #Write the lineage counts
    outLineages = open(args.o + "_lineage_summary.csv", "w")
    outLineages.write("Lineage,Number_of_sequences\n")

    for eL in lC:
        outLineages.write(eL + "," + str(lC[eL]) + "\n")
    
    outNames.close()
    outDesc.close()
    outLineages.close()