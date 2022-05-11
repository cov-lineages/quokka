#Compares pangoLEARN and UShER lineages in an UShER metadata file to identify
#lineages that need investigation

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", help = "Metadata file from UShER")
    parser.add_argument("-o", help = "Name of csv output file containing lineage comparisons")
    args = parser.parse_args()

    #Verify that the headers are in the correct locations
    with open(args.m) as f:
        h = f.readline()
        lC = 0
        uC = 0
        for i, c in enumerate(h.strip().split("\t")):
            if c == "pangolin_lineage":
                lC = i
            elif c == "pango_lineage_usher":
                uC = i
        
        if (lC != 8) or (uC != 10):
            raise RuntimeError("pangolin_lineage and/or pango_lineage_usher columns are not in the expected locations in the metadata file")
    
    #Lineages as keys, lists with shared assignments, pangolin-only assignments and UShER-only assignments as values
    lineages = dict()

    #Iterate through the sequences and add to lineages
    with open(args.m) as f:
        next(f)
        for l in f:
            pL = l.strip().split("\t")[lC]
            uL = l.strip().split("\t")[uC]

            if pL not in lineages:
                lineages[pL] = [0,0,0]
            if uL not in lineages:
                lineages[uL] = [0,0,0]

            if pL == uL:
                lineages[pL][0] += 1
            else:
                lineages[pL][1] += 1
                lineages[uL][2] += 1
    
    outFile = open(args.o, "w")
    outFile.write("Lineage,Assignments_shared_pangoLEARN_UShER,Assignments_specific_to_pangoLEARN,Assignments_specific_to_UShER,Proportion_pangoLEARN_specific,Proportion_UShER_specific\n")

    for l in lineages:
        outFile.write(l + "," + str(lineages[l][0]) + "," + str(lineages[l][1]) + "," + str(lineages[l][2]) + "\n")#+ "," + str(float(lineages[l][1])/(float(lineages[l][0]) + float(lineages[l][1]))) + "," + str(float(lineages[l][2])/(float(lineages[l][0]) + float(lineages[l][2]))) + "\n")
    
    outFile.close()