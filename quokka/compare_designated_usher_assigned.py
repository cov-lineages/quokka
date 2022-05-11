#Compares the number of designated sequences with the number of UShER assigned sequences to each lineage

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", help = "UShER metadata file containing assigned lineages")
    parser.add_argument("-l", help = "lineages.csv containing current designations")
    parser.add_argument("-o", help = "Name of output file containing designation and assignment counts")
    args = parser.parse_args()

    #Verify that the headers in the UShER metadata file are in the correct locations
    with open(args.m) as f:
        h = f.readline()
        uC = 0
        for i, c in enumerate(h.strip().split("\t")):
            if c == "pango_lineage_usher":
                uC = i
        
        if uC != 10:
            raise RuntimeError("pango_lineage_usher column is not in the expected location in the metadata file")

    #Lineages as keys, lists as values containing number of designated and number of assigned sequences
    lineages = dict()

    #Count the number of sequences designated to each lineage
    with open(args.l) as f:
        next(f)
        for l in f:
            sL = l.strip().split(",")[1]
            if sL not in lineages:
                lineages[sL] = [0,0]
            lineages[sL][0] += 1
    
    #Lineages with assignments but no designations
    lD = set()
    
    #Count the number of sequences assigned to each lineage
    with open(args.m) as f:
        next(f)
        for l in f:
            sL = l.strip().split("\t")[uC]
            if sL not in lineages:
                lD.add(sL)
            else:
                lineages[sL][1] += 1
    
    #Write the lineage countrs
    outFile = open(args.o, "w")
    outFile.write("Lineage,Number_designations,Number_assignments,Designation_assignment_ratio\n")

    for eL in lineages:
        outFile.write(eL + "," + str(lineages[eL][0]) + "," + str(lineages[eL][1]) + ",")
        if lineages[eL][1] == 0:
            outFile.write("0\n")
        else:
            outFile.write(str(float(lineages[eL][0])/float(lineages[eL][1])) + "\n")
    
    outFile.close()

    print("Lineages with assignments but no designations:")
    print("\n".join(lD))