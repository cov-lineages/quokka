#Extracts the sequences within a given node of the UShER tree, excludes sequences with reversions, filters
#to remove sequences not in grinch metadata and adds the new lineage to lineages.csv

import argparse
from random import sample
import time

#Identify the mutation(s) leading to the given node
def get_root_mutations(samplePaths, node):
    with open(samplePaths) as f:
        for l in f:
            if node in l:
                nm = list(l.strip().split(node)[1].split(" ")[0].split(","))
                break
    
    return(nm)

#Reverts a given list of mutations
def revert_mutations(mutations):
    rm = list()

    for m in mutations:
        rm.append(m[-1] + m[1:-1] + m[0])
    
    return(rm)

#Extract the sequence names in a given clade, exclude sequences with a reversion of one or more of the
#defining mutations
def get_sequence_names(samplePaths, node, nodeMutations):
    #Revert defining mutations
    mutations = revert_mutations(nodeMutations)

    #Sequence names
    sn = list()
    #All sequences in the clade
    allNames = list()
    
    with open(samplePaths) as f:
        for l in f:
            if node in l:
                allNames.append(l.strip().split("|")[0])
                if not any(m in l for m in mutations):
                    sn.append(l.strip().split("|")[0])
    
    #Convert allNames to set
    allNames = set(allNames)
    
    return(sn, allNames)

#Filters a set of sequence names to keep those present in a given list
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
    parser.add_argument("-s", help = "sample-paths file from UShER")
    parser.add_argument("-m", help = "GISAID metadata file on grinch")
    parser.add_argument("-n", help = "Name of root node of new lineage")
    parser.add_argument("-l", help = "Name of new lineage to which sequences will be designated")
    parser.add_argument("-d", help = "lineages.csv containing current designations")
    parser.add_argument("-o", help = "Name of new lineages.csv containing updated designations")
    args = parser.parse_args()

    #Node to be extracted, append : if not already present to not match nodes with the node name within their name
    node = args.n
    if node[-1] != ":":
        node = node + ":"
    
    #Lineage to which sequences will be designated
    lineage = args.l
    
    #Identify the mutation(s) leading to the root node
    print("Identifying mutations on root node")
    nodeMutations = get_root_mutations(args.s, node)
    print("Mutations leading to " + args.n + ": " + ",".join(nodeMutations))

    #Extract the sequence names in the clade, exclude sequences with a reversion of one or more of the
    #defining mutations
    print("Extracting sequences")
    sequenceNames, allNames = get_sequence_names(args.s, node, nodeMutations)
    print(str(len(allNames)) + " sequences in clade " + args.n + " of which " + str(len(sequenceNames)) + " do not have a reversion")

    #Filter the sequence names that are in the metadata
    print("Filtering sequence names")
    filteredNames = filter_sequence_names(sequenceNames, args.m)
    print(str(len(filteredNames)) + " sequences post-filter")

    #Sequences already in lineages.csv, used to ensure sequences in -s are not written twice
    included = list()

    outFile = open(args.o, "w")

    #New designations
    nd = 0
    #Updated designations
    ud = 0
    #Previous designations of updated
    prevD = list()

    #Iterate through the current designations, check if they are in the filtered names, update and write if so, write if not
    print("Adding designations")
    with open(args.d) as f:
        for l in f:
            sName = l.strip().split(",")[0]
            if sName in filteredNames:
                prevD.append(l.strip().split(",")[1])
                outFile.write(sName + "," + lineage + "\n")
                included.append(sName)
                ud += 1
            #Check if the sequence is in the clade but not the filtered sequences, if so do not include
            elif sName in allNames:
                continue
            #The sequence is outside of the lineage clade
            else:
                outFile.write(l)

    #Convert included to set
    included = set(included)
    
    #Iterate through the filtered sequence names and add any that haven't already been added
    for s in filteredNames:
        if s not in included:
            outFile.write(s + "," + lineage + "\n")
            nd += 1

    outFile.close()

    print(str(nd + ud) + " sequences designated to " + lineage + " of which " + str(nd) + " are new designations and " + str(ud) + " were updated")

    print("Summary of updated designations:")

    for eL in set(prevD):
        print(str(prevD.count(eL)) + " designations updated from " + eL)