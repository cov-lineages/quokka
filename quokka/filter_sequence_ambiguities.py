#Counts the number of ambiguities in each sequence in an alignment
#Outputs a csv file containing the ambiguity content of each sequence and an alignment filtered to keep only
#sequences with less than a given proportion of ambiguities, default 5%
#Gaps are treated as real deletions so are not treated as ambiguities
#Expects an alignment that has been created with pangolin so is padded with Ns at the 5' and 3' end
#Ambiguities are only counted in the coding region between positions 266 and 29674
#Output csv file contains: sequence name, number ambiguities, proportion ambiguities

import argparse
from re import A
from Bio import AlignIO
from collections import Counter

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help = "fasta alignment file from pangolin")
    parser.add_argument("-p", help = "Proportion ambiguities to filter a sequence, default 5", default = "5")
    parser.add_argument("-o", help = "Output prefix, will output this with .csv containing ambiguity stats and with " + 
                                    ".fasta containing a filtered alignment")
    args = parser.parse_args()

    #Import the alignment
    alignment = AlignIO.read(args.a, "fasta")

    outFile = open(args.o + ".csv", "w")
    outFile.write("sequence,number_ambiguities,proportion_ambiguities\n")

    outFasta = open(args.o + ".fasta", "w")

    #Proportion of ambiguities
    ambiguities = float(args.p)

    #Length of the coding region
    gL = float(29409)

    #Iterate through the sequences, count their ambiguities and write
    for s in alignment:
        nC = Counter(s.seq[265:29674])
        #Number of ambiguities
        nA = float(29409 - (nC["A"] + nC["C"] + nC["G"] + nC["T"] + nC["-"]))
        #Proportion of ambiguities
        pA = (nA/gL) * float(100)

        #Write to the csv file
        outFile.write(s.id + "," + str(nA) + "," + str(pA) + "\n")

        #If the sequence contains <args.p% gaps, write to the alignment
        if pA <= ambiguities:
            outFasta.write(">" + s.id + "\n" + str(s.seq) + "\n")

    #print(alignment)

    outFile.close()
    outFasta.close()