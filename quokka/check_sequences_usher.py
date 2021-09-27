#Extracts the sequences from the UShER tree in a clade of interest and checks whether sequences designated to a given lineage are outside the
#lineage clade in the UShER tree and whether sequences designated outside the lineage are within the lineage clade in the UShER tree
#Used to check whether sequences designated to a lineage are in the lineage's clade in the UShER tree
#Takes a sample-paths file, created from a protobuf file using matUtils: matUtils extract -i tree.pb -S sample-paths
#Outputs 2 files: cluster_outside.txt containing sequences designated to the lineage that don't cluster within the lineage and cluster_inside.txt containing
#sequences designated outside the lineage that cluster within the lineage

import argparse
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help = "sample-paths containing sequences and their paths in the tree")
    parser.add_argument("-n", nargs = "+", help = "Names of sequences to examine")
    parser.add_argument("-m", help = "Mutation to define the root node of the clade")
    parser.add_argument("-d", help = "lineages.csv containing current sequence designations")
    parser.add_argument("-l", help = "Lineage of interest, sequences designated to this lineage will be checked")
    parser.add_argument("-o", help = "Prefix of output files")
    args = parser.parse_args()

    #Sequence names
    sequences = args.n

    #Mutation of interest
    mutation = args.m

    #Lineage of interest
    lineage = args.l

    outFileOutside = open(args.o + "_cluster_outside.txt", "w")
    outFileInside = open(args.o + "_cluster_inside.txt", "w")

    #Nodes with the mutation of interest leading to the samples of interest
    nodes = list()

    #Iterate through the sample paths to extract the node at which the mutation of interest occurs in each sample of interest
    with open(args.s) as fileobject:
        for line in fileobject:
            if re.sub("\\|.*", "", line.strip().split("\t")[0]) in sequences:
                if mutation not in line:
                    print("Mutation " + mutation + " not in the path to " + re.sub("\\|.*", "", line.strip().split("\t")[0]))
                else:
                    nodes.append(line.strip().split(mutation)[0].split(" ")[-1].split(":")[0])
    
    #Sequences in the clade of interest
    clade = list()
    
    #Check if the mutation occurs at the same node in each case
    if not len(set(nodes)) == 1:
        print("Mutation " + mutation + " occurs on different branches in the provided samples")
    else:
        node = nodes[0]
        #Iterate through the sequences, check if they contain the node in their path and add to the list of sequences in the clade if so
        with open(args.s) as fileobject:
            for line in fileobject:
                if node + ":" in line.strip():
                    clade.append(re.sub("\\|.*", "", line.strip().split("\t")[0]))
    
    #Iterate through the designated sequences, check if they are designated to the lineage of interest and not in the clade or designated to
    #another lineage and in the clade
    with open(args.d) as fileobject:
        next(fileobject)
        for line in fileobject:
            if (line.strip().split(",")[1] == lineage) and (line.strip().split(",")[0] not in clade):
                outFileOutside.write(line.strip().split(",")[0] + "\n")
            if (line.strip().split(",")[1] != lineage) and (line.strip().split(",")[0] in clade):
                outFileInside.write(line.strip().split(",")[0] + "\n")
    
    outFileOutside.close()
    outFileInside.close()