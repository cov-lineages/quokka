#Extracts the defining mutations for a given lineage
#These mutations are those that are acquired leading to the root of the lineage and are present in
#at least a given percentage of sequences designated to the lineage, default 70%

import argparse
import re
from collections import Counter
from Bio.Seq import Seq
import json
import array
from type_constellations import load_feature_coordinates


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

#Adds each of a series of mutations to the reference sequence
def addToReference(ref, mutations):
    #Don't update the original reference to enable comparison
    uRef = ref

    #Iterate through the mutations and add to the reference
    for m in mutations:
        p = int(m[1:-1])
        uRef = uRef[:(p - 1)] + m[-1] + uRef[p:]
    
    return(uRef)

#Translates each of a set of genes given a nucleotide sequence and gene coordinates
def translate(sequence, genes, nsp = False):
    geneDict = dict()

    for g in genes:
        geneDict[g] = Seq(sequence[(genes[g][0] - 1):genes[g][1]]).translate()
    
    return(geneDict)

#Identifies whether each mutation in a list is synonymous or nonsynonymous
#This is determined in the context of all mutations in the list so 2 or more mutations in
#the same codon are examined together
def getEffect(ref, uRef, genes, mutations, nsp = False):
    #Extract the protein sequences for the original and updated references
    refGenes = translate(ref, genes)
    uRefGenes = translate(uRef, genes)

    #Positions as keys, effect as values
    pDict = dict()

    #Iterate through the mutations, get their amino acid position and check if the amino acid has changed
    for m in mutations:
        #Will change to False if the mutation is in a gene
        intergenic = True

        p = int(m[1:-1])
        
        for g in genes:
            if "nsp" not in g:
                if (p >= genes[g][0]) and (p <= genes[g][1]):
                    intergenic = False
                    #Nucleotide position of the mutation in the gene, zero based
                    gP = p - genes[g][0]
                    #Amino acid position in the gene, zero based
                    aaP = int(gP/3)

                    rAA = refGenes[g][aaP]
                    uAA = uRefGenes[g][aaP]

                    if rAA == uAA:
                        pDict[p] = "synonymous"
                    else:
                        pDict[p] = g + ":" + rAA + str(aaP + 1) + uAA
                    break
    
        #Add intergenic mutations
        if intergenic:
            pDict[p] = "synonymous"
    
    return(pDict)

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

#Creates a dictionary to write as a constellation
def getConstellation(mutations, lineage):
    con_dict = dict()
    con_dict["label"] = lineage
    con_dict["description"] = lineage + " defining mutations - acquired along the path to the root of the lineage and remain conserved within the lineage"
    con_dict["sources"] = list()
    con_dict["type"] = "lineage"
    con_dict["tags"] = [lineage]

    sites = list()
    for m in mutations:
        if "synonymous" in m:
            sites.append("nuc:" + m.split("(")[0])
        else:
            sites.append(m.split("(")[1].split(")")[0])
    con_dict["sites"] = sites

    return(con_dict)

def convertMutation(m):
    if m.split(":")[0] == "orf1a":
        return("Orf1ab:" + m.split(":")[1])
    elif m.split(":")[0] == "orf1b":
        return("Orf1ab:" + m.split(":")[1][0] + str(int(m.split(":")[1][1:-1]) + 4401) + m[-1])
    elif m.split(":")[0] == "s":
        return("S:" + m.split(":")[1])
    elif m.split(":")[0] == "orf3a":
        return("Orf3a:" + m.split(":")[1])
    elif m.split(":")[0] == "e":
        return("E:" + m.split(":")[1])
    elif m.split(":")[0] == "m":
        return("M:" + m.split(":")[1])
    elif m.split(":")[0] == "orf6":
        return("Orf6:" + m.split(":")[1])
    elif m.split(":")[0] == "orf7a":
        return("Orf7a:" + m.split(":")[1])
    elif m.split(":")[0] == "orf8":
        return("Orf8:" + m.split(":")[1])
    elif m.split(":")[0] == "n":
        return("N:" + m.split(":")[1])
    elif m.split(":")[0] == "orf10":
        return("Orf10:" + m.split(":")[1])
    else:
        return(m)

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
    parser.add_argument("-r",
                        "--reference",
                        dest = "ref_json",
                        help = "Reference json file containing genome sequence and gene coordinates, default SARS-CoV-2.json " + 
                        "assumed to be in the working directory",
                        default = "SARS-CoV-2.json")
    
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
                        help = "Prefix of output files to which the defining mutations will be written. 3 files are saved:\n" + 
                        "prefix_consensus.txt - contains the consensus mutations and their effects\n" + 
                        "prefix_consensus_constellation.txt - contains the consensus mutations in constellation format\n" + 
                        "prefix_mutation_proportions.csv - contains all of the mutations acquired leading to the root of the clade (excluding " + 
                        "back mutations), including those that don't reach the consensus threshold. Proportions are the percentage of sequences " + 
                        "designated to the lineage with the mutation")
    
    args = parser.parse_args()

    #Extract the sequences designated to the lineage of interest
    lineage = args.designation
    sequences = list()
    with open(args.lineages) as fileobject:
        for line in fileobject:
            if line.strip().split(",")[1] == lineage:
                sequences.append(line.strip().split(",")[0])
    
    #Extract the mutations to the root of the lineage
    mutations = extractMutations(args.sample_paths, sequences, args.node, args.mutation)

    #Filter the mutations so positions that have mutated multiple times keep the most recent mutation
    fM = filterMutations(mutations)

    #Extract the sequence and gene coordinates of the reference
    ref = load_feature_coordinates(args.ref_json)

    #Add the mutations to the reference sequence
    uRef = addToReference(ref[0], fM)

    #Identify whether each mutation is synonymous or nonsynonymous relative to the reference
    mEffect = getEffect(ref[0], uRef, ref[1], mutations)

    #Calculate the proportion of designated sequences with each mutation
    positions = extractNucleotides(args.sample_paths, sequences, fM)

    outFile = open(args.outFile + "_consensus.txt", "w")
    outFileProportions = open(args.outFile + "_mutation_proportions.csv", "w")
    outFileProportions.write("Mutation,Proportion_of_designated_sequences\n")
    outFileConstellation = open(args.outFile + "_consensus_constellation.json", "w")

    #Mutations to be written to the constellation
    cM = list()

    #Iterate through the mutations, check if they are in a sufficiently high proportion of designated sequences and write if so
    for p in sorted(positions.keys()):
        sC = Counter(positions[p][2])
        mP = (float(sC[positions[p][1]])/float(len(positions[p][2]))) * float(100)
        if mP >= float(args.proportion):
            conMut = convertMutation(mEffect[p])
            outFile.write(positions[p][0] + "(" + conMut + ")\n")
            cM.append(positions[p][0] + "(" + mEffect[p] + ")")
        outFileProportions.write(positions[p][0] + "(" + conMut + ")," + str(mP) + "\n")
    
    #Write the constellation file
    con_dict = getConstellation(cM, args.designation)
    json.dump(con_dict, outFileConstellation, indent = 4)

    outFile.close()
    outFileProportions.close()
    outFileConstellation.close()