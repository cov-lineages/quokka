#Checks each of a given set of positions for mixtures in a bam file

import pysam
from collections import Counter
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", help = "BAM file to be examined")
    parser.add_argument("-p", nargs = "+", help = "Set of positions to be checked")
    parser.add_argument("-r", help = "Name of reference the reads were mapped against, default MN908947.3", default = "MN908947.3")
    parser.add_argument("-o", help = "Name of output file to which position summary will be written")
    args = parser.parse_args()

    #Import BAM
    bam = pysam.Samfile(args.b, "rb")

    #Reference sequence name
    reference = args.r

    outFile = open(args.o, "w")
    outFile.write("Position,Total_reads,A_reads,C_reads,G_reads,T_reads,A_proportion,C_proportion,G_proportion,T_proportion\n")

    #Iterate through the positions, identify the nucleotide in each read and count each base
    for p in args.p:
        position = int(p)

        bases = ""

        for c in bam.pileup(reference, (position - 1), position):
            for r in c.pileups:
                if c.pos == (position - 1):
                    if r.query_position is not None:
                        bases += r.alignment.query_sequence[r.query_position]
        
        bc = Counter(bases)

        tr = float(sum(bc.values()))

        outFile.write(p + "," + str(tr) + ",")
        if "A" in bc:
            outFile.write(str(bc["A"]) + ",")
        else:
            outFile.write("0,")
        if "C" in bc:
            outFile.write(str(bc["C"]) + ",")
        else:
            outFile.write("0,")
        if "G" in bc:
            outFile.write(str(bc["G"]) + ",")
        else:
            outFile.write("0,")
        if "T" in bc:
            outFile.write(str(bc["T"]) + ",")
        else:
            outFile.write("0,")
        
        if "A" in bc:
            outFile.write(str(float(bc["A"])/tr) + ",")
        else:
            outFile.write("0,")
        if "C" in bc:
            outFile.write(str(float(bc["C"])/tr) + ",")
        else:
            outFile.write("0,")
        if "G" in bc:
            outFile.write(str(float(bc["G"])/tr) + ",")
        else:
            outFile.write("0,")
        if "T" in bc:
            outFile.write(str(float(bc["T"])/tr) + "\n")
        else:
            outFile.write("0\n")
    
    outFile.close()