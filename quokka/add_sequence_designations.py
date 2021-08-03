#Adds a new set of designations to lineages.csv
#Takes lineages.csv with -l
#Takes a file containing all of the sequences to be designated with -s
#The lineage to which these sequences should be designated is provided with -d
#If a sequence is already in lineages.csv, its designation will be updated to -d
#Sequences in the file provided with -s that are not already in lineages.csv will be added with
#the designation provided with -d
#The new lineages.csv will be written to the csv file name provided with -o
#Differs from add_new_designations.py as that does not update sequences already in lineages.csv but this script does

import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", help = "lineages.csv")
    parser.add_argument("-s", help = "File containing sequences to be designated, 1 sequence name per row, no header")
    parser.add_argument("-d", help = "Name of lineage to which sequences in -s will be designated")
    parser.add_argument("-o", help = "Output csv file containing new designations")
    args = parser.parse_args()
    
    #Sequences to be designated, used to check which sequences in lineages.csv need to be updated
    sequences = list()
    #Iterate through the sequences to be designated and add to sequences
    with open(args.s) as fileobject:
        for line in fileobject:
            sequences.append(line.strip())
    
    #Sequences already in lineages.csv, used to ensure sequences in -s are not written twice
    included = list()

    #Import lineages.csv
    lineages = pd.read_csv(args.l)

    outFile = open(args.o, "w")
    outFile.write(lineages.columns[0] + "," + lineages.columns[1] + "\n")

    #Iterate through the designated sequences, check if they need to be updated, if so update, else write as they are
    #Add all sequences to included
    for i in range(lineages.shape[0]):
        if lineages["taxon"][i] in sequences:
            outFile.write(lineages["taxon"][i] + "," + args.d + "\n")
        else:
            outFile.write(lineages["taxon"][i] + "," + lineages["lineage"][i] + "\n")
        included.append(lineages["taxon"][i])
    
    #Iterate through the sequences to be designated, check if they are already designated and add if not
    with open(args.s) as fileobject:
        for line in fileobject:
            if line.strip() not in included:
                outFile.write(line.strip() + "," + args.d + "\n")

    outFile.close()