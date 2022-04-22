#Adds a set of random designations from a given lineage
#Identifies the current designations to the lineage from the provided lineages.csv file
#Extracts all sequences assigned by UShER to the lineage from the provided UShER metadata file
#Filters the sequences to keep those not currently designated
#Selects a given number of random sequences from these
#Filters these sequences based on the current all samples file
#Adds to lineages.csv

import random
import argparse

#Extracts designations to a given lineage from a lineages.csv file
def getDesignations(lFile, lineage):
    d = list()

    #Iterate through the designations, check if they are to the required lineage, if so add to d
    with open(lFile) as f:
        for l in f:
            if l.strip().split(",")[1] == lineage:
                d.append(l.strip().split(",")[0])
    
    d = set(d)
    return(d)

#Extracts the sequences assigned to a lineage in the UShER tree that are not currently designated
def getNewAssignments(mFile, d, lineage):
    #Verify pango_lineage_usher column in UShER metadata file
    with open(mFile) as f:
        h = f.readline()
        uC = 0
        for i, c in enumerate(h.strip().split("\t")):
            if c == "pango_lineage_usher":
                uC = i
        
        if uC != 10:
            raise RuntimeError("pango_lineage_usher column is not in the expected location in the metadata file")
    
    assignments = list()

    #Iterate through the sequences, check if their UShER lineage is lineage and add to l if they are not in d
    with open(mFile) as f:
        for l in f:
            if l.strip().split("\t")[uC] == lineage:
                if l.strip().split("\t")[0].split("|")[0] not in d:
                    assignments.append(l.strip().split("\t")[0].split("|")[0])
    
    return(assignments)

#Filters a set of sequence names to keep those present in a given list
#From designate_sequences_usher.py
def filter_sequence_names(sequenceNames, metadata):
    #Filtered sequence names
    fn = list()

    sequenceNames = set(sequenceNames)

    with open(metadata) as f:
        for l in f:
            if l.strip().split(",")[0] in sequenceNames:
                fn.append(l.strip().split(",")[0])
    
    #Convert fn to set
    fn = set(fn)
    
    return(fn)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", help = "UShER metadata file containing assigned lineages")
    parser.add_argument("-s", help = "all_samples file from grinch containing filtered sequences that can be designated")
    parser.add_argument("-d", help = "lineages.csv containing current designations")
    parser.add_argument("-l", help = "Lineage to add new designations to")
    parser.add_argument("-n", help = "Number of sequences to extract. This may not be the number designated as this " + 
                        "number will be extracted and then filtered")
    parser.add_argument("-o", help = "Name of new lineages.csv containing updated designations")
    args = parser.parse_args()

    #Extract the current designations to the lineage
    print("Extracting current " + args.l + " designations")
    currentD = getDesignations(args.d, args.l)
    print("Currently " + str(len(currentD)) + " designations to " + args.l)

    #Extract the UShER assignments to the lineage that are not currently designated
    print("Extracting UShER assignments to " + args.l + " that are not currently designated")
    uA = getNewAssignments(args.m, currentD, args.l)
    print(str(len(uA)) + " sequences currently assigned to " + args.l + " that are not designated")

    #Extract a random subsample of these sequences
    ss = random.sample(uA, int(args.n))

    #Filter random subsample
    print("Filtering " + args.n + " sequences")
    filteredNames = filter_sequence_names(ss, args.s)
    print(str(len(filteredNames)) + " sequences post-filter")

    lineage = args.l

    #Add sequences to lineages.csv
    print("Adding designations")
    outFile = open(args.o, "w")
    with open(args.d) as f:
        for l in f:
            outFile.write(l)
    
    for s in filteredNames:
        outFile.write(s + "," + lineage + "\n")
    
    outFile.close()