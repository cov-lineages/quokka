#Extracts the sequences from the UShER tree in a clade of interest and checks whether samples within a given list are within this clade
#Used to check whether sequences designated to a lineage are in the lineage's clade in the UShER tree
#Takes a sample-paths file, created from a protobuf file using matUtils: matUtils extract -i tree.pb -S sample-paths
#Prints sample names that are in the designations for a lineage but are not in the lineage in the UShER tree

import argparse
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help = "sample-paths containing sequences and their paths in the tree")
    parser.add_argument("-n", nargs = "+", help = "Names of sequences to examine")
    parser.add_argument("-m", help = "Mutation to define the root node of the clade")
    parser.add_argument("-d", help = "Sequence designations to be checked")
    parser.add_argument("-o", help = "File to which sequences in clade of interest will be written")
    args = parser.parse_args()

    #Sequence names
    sequences = args.n

    #Mutation of interest
    mutation = args.m

    outFile = open(args.o, "w")

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
        #Iterate through the sequences, check if they contain the node in their path and add to the list of sequences if so
        with open(args.s) as fileobject:
            for line in fileobject:
                if node + ":" in line.strip():
                    sN = re.sub("\\|.*", "", line.strip().split("\t")[0])
                    clade.append(sN)
                    outFile.write(sN)
    
    #Iterate through the designated sequences, check if they are in the clade and print if not
    with open(args.d) as fileobject:
        for line in fileobject:
            if line.strip() not in clade:
                print(line.strip())
    
    outFile.close()