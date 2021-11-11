#Extracts the defining mutations for a given lineage
#These mutations are those that are acquired leading to the root of the lineage and are present in
#at least a given percentage of sequences designated to the lineage, default 70%

import argparse
import re
from collections import Counter

#Extracts the mutations up to a node from a node string
def extractMutationsToNode(path, node):
    mutations = list()
    
    #Iterate through the mutations upstream of the node and add to mutations
    for m in path.split(node)[0].split():
        for eM in m.split(":")[1].split(","):
            mutations.append(eM)
    
    #Iterate through the mutations on the node and add to mutations
    for nM in path.split(node + ":")[1].split()[0].split(","):
        mutations.append(nM)
    
    return(mutations)

#Extracts the mutations that have occurred up to a node of interest
def extractMutations(paths, sequences, node, mutation):
    #If a node is given, extract mutations up to that node
    if node:
        with open(paths) as fileobject:
            for line in fileobject:
                if (node + ":") in line:
                    mutations = extractMutationsToNode(line.strip().split("\t")[1], node)
                    break
    #If the node is unknown, infer the root node of the lineage from the given mutation
    else:
        #Check the mutation occurs in the same place in 5 designated sequences
        s = 0
        nodes = list()
        with open(paths) as fileobject:
            for line in fileobject:
                if s == 5:
                    break
                else:
                    if re.sub("\\|.*", "", line.strip().split("\t")[0]) in sequences:
                        if mutation not in line:
                            raise RuntimeError("Mutation " + mutation + " not in the path to " + re.sub("\\|.*", "", line.strip().split("\t")[0]))
                        else:
                            nodes.append(line.strip().split(mutation)[0].split(" ")[-1].split(":")[0])
                        s += 1
        
        #Check if the mutation occurs on the same node in each sequence, if so extract the mutations to this node
        if len(set(nodes)) != 1:
            raise RuntimeError("The given mutation occurs at different nodes in different sequences")
        else:
            rootNode = nodes[0]
            with open(paths) as fileobject:
                for line in fileobject:
                    if (rootNode + ":") in line:
                        mutations = extractMutationsToNode(line.strip().split("\t")[1], rootNode)
                        break
    
    return(mutations)

#Filters positions that have mutated multiple times to keep only the most recent mutation
def filterMutations(mutations):
    #Examined positions
    p = list()
    #Mutations to keep
    mL = list()

    #Iterate through the mutations, check if their position is already in p, if so the position
    #mutates later in the path so the earlier mutation is excluded. If the position is not already
    #present, add to mL
    for m in mutations[::-1]:
        if m[1:-1] not in p:
            mL.append(m)
        p.append(m[1:-1])
    
    return(mL[::-1])

#Calculates the proportion of a given set of sequences that have each mutation
def extractNucleotides(paths, sequences, mutations):
    #Mutations as keys, sequences with the mutation as values
    sM = dict()
    for m in mutations:
        sM[int(m[1:-1])] = [m, m[-1], []]
    
    #Iterate through the sequences, check if they are in the designated sequences and extract their
    #nucleotide at each position if so
    with open(paths) as fileobject:
        for line in fileobject:
            if re.sub("\\|.*", "", line.strip().split("\t")[0]) in sequences:
                for p in sM:
                    #Extract the nucleotide of the last match to the position in the line
                    sM[p][2].append(re.findall("[A-Z]" + str(p) + "[A-Z]", line.strip())[-1][-1])
    
    return(sM)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",
                        "--sample_paths",
                        dest = "sample_paths",
                        required = True,
                        help = "sample-paths file from matUtils")
    parser.add_argument("-l",
                        "--lineages",
                        dest = "lineages",
                        required = True,
                        help = "lineages.csv containing designated sequences and their lineage")
    parser.add_argument("-d",
                        "--designation",
                        dest = "designation",
                        required = True,
                        help = "Name of the lineage to be examined")
    parser.add_argument("-p",
                        "--proportion",
                        dest = "proportion",
                        help = "Minimum proportion of designated sequences required for the mutation " + 
                        "to be included in the consensus, deafult 70",
                        default = "70")
    
    mutations = parser.add_mutually_exclusive_group(required = True)
    mutations.add_argument("-n",
                            "--node",
                            dest = "node",
                            help = "The node number on which the lineage starts. The mutations up to and " + 
                            "including this node will be used for the consensus. Either -n or -m is required",
                            default = None)
    mutations.add_argument("-m",
                            "--mutation",
                            dest = "mutation",
                            help = "The mutation that defines the start of the lineage. Use this if the node at " + 
                            "which the lineage starts is unknown. This mutation will be searched for in the path. If the " + 
                            "mutation occurs multiple times in the path, the last occurrence will be used. If there are multiple " + 
                            "mutations on the branch leading to the lineage, only provide one of those. Either -n or -m is required",
                            default = None)
    
    parser.add_argument("-o",
                        "--outfile",
                        dest = "outFile",
                        required = True,
                        help = "Prefix of output files to which the defining mutations will be written")
    
    args = parser.parse_args()

    #Extract the sequences designated to the lineage of interest
    lineage = args.designation
    sequences = list()
    with open(args.lineages) as fileobject:
        for line in fileobject:
            if line.strip().split(",")[1] == lineage:
                sequences.append(line.strip().split(",")[0])
    
    #Extract the mutations to the root of the lineages
    mutations = extractMutations(args.sample_paths, sequences, args.node, args.mutation)

    #Filter the mutations so positions that have mutated multiple times keep the most recent mutation
    fM = filterMutations(mutations)

    #Calculate the proportion of designated sequences with the mutation
    positions = extractNucleotides(args.sample_paths, sequences, fM)

    outFile = open(args.outFile + "_consensus.txt", "w")
    outFileProportions = open(args.outFile + "_mutation_proportions.csv", "w")
    outFileProportions.write("Mutation,Proportion_of_designated_sequences\n")

    #Iterate through the mutations, check if they are in a sufficiently high proportion of designated sequences and write if so
    for p in sorted(positions.keys()):
        sC = Counter(positions[p][2])
        mP = (float(sC[positions[p][1]])/float(len(positions[p][2]))) * float(100)
        if mP >= float(args.proportion):
            outFile.write(positions[p][0] + "\n")
        outFileProportions.write(positions[p][0] + "," + str(mP) + "\n")

    outFile.close()
    outFileProportions.close()