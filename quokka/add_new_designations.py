#Adds a new set of designations to lineages.csv
#Takes lineages.csv with -l
#Takes a file containing all of the sequences to be designated with -s
#If a sequence is already in lineages.csv, it will retain its designation from lineages.csv and
#will not be added again
#Sequences in the file provided with -s will be added to lineages.csv with the designation provided with -d
#The new lineages.csv will be written to the csv file name provided with -o

import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", help = "lineages.csv")
    parser.add_argument("-s", help = "File containing sequences to be designated, 1 sequence name per row, no header")
    parser.add_argument("-d", help = "Name of lineage to which sequences in -s will be designated")
    parser.add_argument("-o", help = "Output csv file containing new designations")
    args = parser.parse_args()
    
    #Sequences already in lineages.csv
    included = list()

    #Import lineages.csv
    lineages = pd.read_csv(args.l)

    outFile = open(args.o, "w")
    outFile.write(lineages.columns[0] + "," + lineages.columns[1] + "\n")

    #Iterate through the designated sequences,  write as they are and add to included
    for i in range(lineages.shape[0]):
        outFile.write(lineages["taxon"][i] + "," + lineages["lineage"][i] + "\n")
        included.append(lineages["taxon"][i])
    
    #Iterate through the sequences to be designated, check if they are already designated and add if not
    with open(args.s) as fileobject:
        for line in fileobject:
            if line.strip() not in included:
                outFile.write(line.strip() + "," + args.d + "\n")

    outFile.close()