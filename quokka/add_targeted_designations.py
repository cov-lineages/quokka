#Adds a new set of designations to a given lineage
#Extracts sequences whose UShER path contains a given node but not one or more other nodes
#Filters these sequence names based on the grinch metadata file
#Randomly subsamples a given number of these and designates to a given lineage

import random
import argparse

#Identify the mutation(s) leading to the given node
#From designate_sequences_usher.py
def get_root_mutations(samplePaths, node):
    with open(samplePaths) as f:
        for l in f:
            if node in l:
                nm = list(l.strip().split(node)[1].split(" ")[0].split(","))
                break
    
    return(nm)

#Reverts a given list of mutations
#From designate_sequences_usher.py
def revert_mutations(mutations):
    rm = list()

    for m in mutations:
        rm.append(m[-1] + m[1:-1] + m[0])
    
    return(rm)

#Extract the sequence names in a given clade, exclude sequences with a reversion of one or more of the
#defining mutations
#Adapted from get_sequence_names in designate_sequences_usher.py by adding an additional filter to exclude sequences in any of a given
#set of downstream nodes
def extract_sequence_names(samplePaths, node, exclude, nodeMutations):
    #Revert defining mutations
    mutations = revert_mutations(nodeMutations)

    #Sequence names
    sn = list()
    
    with open(samplePaths) as f:
        for l in f:
            if node in l:
                #Exclude sequences with a reversion or in one of the clades to exclude
                if (not any(m in l for m in mutations)) and (not any(n in l for n in exclude)):
                    sn.append(l.strip().split("|")[0])
    
    return(sn)

#Filters a set of sequence names to keep those present in a given list
#Adapted from designate_sequences_usher.py to remove the final step of converting fn to a set and random.sample takes a list
def filter_sequence_names(sequenceNames, metadata):
    #Filtered sequence names
    fn = list()

    sequenceNames = set(sequenceNames)

    with open(metadata) as f:
        for l in f:
            if l.strip().split(",")[0] in sequenceNames:
                fn.append(l.strip().split(",")[0])
    
    return(fn)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help = "sample-paths file from UShER")
    parser.add_argument("-m", help = "all_samples file from grinch containing filtered sequences that can be designated")
    parser.add_argument("-d", help = "lineages.csv containing current designations")
    parser.add_argument("-n", help = "Name of root node of the lineage")
    parser.add_argument("-e", nargs = "+", help = "One or more nodes to exclude from the designations")
    parser.add_argument("-l", help = "Lineage to add new designations to")
    parser.add_argument("-r", help = "Number of sequences to designate. This number of sequences from the filtered list will be randomly sampled")
    parser.add_argument("-o", help = "Name of new lineages.csv containing updated designations")
    args = parser.parse_args()

    #Node to be extracted, append : if not already present to not match nodes with the node name within their name
    node = args.n
    if node[-1] != ":":
        node = node + ":"
    #Nodes to be excluded, append : if not already present
    exclude = list()
    for e in args.e:
        if e[-1] != ":":
            exclude.append(e + ":")
        else:
            exclude.append(e)
    
    #Lineage to which sequences will be designated
    lineage = args.l

    #Identify the mutation(s) leading to the root node
    print("Identifying mutations on root node")
    nodeMutations = get_root_mutations(args.s, node)
    print("Mutations leading to " + args.n + ": " + ",".join(nodeMutations))

    #Extract the sequence names in the clade, exclude sequences with any of the exclude nodes or with a reversion
    #of one or more of the defining mutations
    print("Extracting sequences")
    sequenceNames = extract_sequence_names(args.s, node, exclude, nodeMutations)
    print(str(len(sequenceNames)) + " sequences in clade " + args.n + " excluding sequences with a reversion and sequences in " + ",".join(args.e))

    #Filter the sequence names that are in the metadata
    print("Filtering sequence names")
    filteredNames = filter_sequence_names(sequenceNames, args.m)
    print(str(len(filteredNames)) + " sequences post-filter")

    #Randomly sample sequences to designate
    ss = random.sample(filteredNames, int(args.r))
    ss = set(ss)
    print("Randomly sampled " + args.r + " sequences")

    #Sequences already in lineages.csv, used to ensure sequences are not written twice
    included = list()
    #New designations
    nd = 0
    #Updated designations
    ud = 0
    #Previous designations of updated
    prevD = list()

    #Add new designations
    print("Designating sequences")
    outFile = open(args.o, "w")
    with open(args.d) as f:
        for l in f:
            sName = l.strip().split(",")[0]
            if sName in ss:
                prevD.append(l.strip().split(",")[1])
                outFile.write(sName + "," + lineage + "\n")
                included.append(sName)
                ud += 1
            #The sequence is outside of the lineage clade
            else:
                outFile.write(l)
    
    #Convert included to set
    included = set(included)
    
    #Iterate through the filtered sequence names and add any that haven't already been added
    for s in ss:
        if s not in included:
            outFile.write(s + "," + lineage + "\n")
            nd += 1
    
    outFile.close()

    print(str(nd + ud) + " sequences designated to " + lineage + " of which " + str(nd) + " are new designations and " + str(ud) + " were updated")

    print("Summary of updated designations:")

    for eL in set(prevD):
        print(str(prevD.count(eL)) + " designations updated from " + eL)