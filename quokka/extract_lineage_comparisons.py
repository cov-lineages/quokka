#Takes the latest lineages.csv file and metadata.tsv from GISAID
#Identifies the number of assigned and designated sequences for each lineage
#Also identifies the number of sequences from each lineage collected over the previous 4 months
#Outputs a csv file containing the number of designated, assigned and recent assigned sequences for each lineage
#To run: python3 extract_lineage_comparisons.py -l lineages.csv -m metadata,tsv -o designated_assigned_comparison.csv

import argparse
import csv
from datetime import date
from dateutil.relativedelta import relativedelta

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", help = "lineages.csv from pango-designation")
    parser.add_argument("-m", help = "metadata.tsv from GISAID")
    parser.add_argument("-o", help = "Name of output csv file")
    args = parser.parse_args()

    #Keys are lineages, values are number of designated sequences
    designated = dict()
    #Keys are lineages, values are number of assigned sequences
    assigned = dict()
    #Keys are lineages, values are sequences from the past 4 months
    recent = dict()

    #Cutoff date for recent sequences
    cutoff = date.today() - relativedelta(months = 4)

    i = 0

    #Iterate through the current designations and add to designated
    with open(args.l) as fileobject:
        lReader = csv.reader(fileobject)
        #Skip the header
        next(lReader, None)

        for row in lReader:
            if row[1] in designated:
                designated[row[1]] += 1
            else:
                designated[row[1]] = 1
    
    #Iterate through the assignments and add to assigned and recent
    with open(args.m) as fileobject:
        mReader = csv.reader(fileobject, delimiter = "\t")
        #Skip the header
        next(mReader, None)

        for row in mReader:
            if row[11] in assigned:
                assigned[row[11]] += 1
            else:
                assigned[row[11]] = 1
            
            #Check if the sequence is recent
            if row[3].count("-") == 2:
                if date(int(row[3].split("-")[0]), int(row[3].split("-")[1]), int(row[3].split("-")[2])) >= cutoff:
                    if row[11] in recent:
                        recent[row[11]] += 1
                    else:
                        recent[row[11]] = 1
    
    #Write the sequence counts
    outFile = open(args.o, "w")
    outFile.write("lineage,designated_sequences,assigned_sequences,assigned_designated_ratio,assigned_in_past_4_months\n")

    for lineage in designated:
        outFile.write(lineage + "," + str(designated[lineage]) + ",")

        if lineage in assigned:
            outFile.write(str(assigned[lineage]) + "," + str(float(assigned[lineage])/float(designated[lineage])) + ",")
        else:
            outFile.write("0,0,")
        
        if lineage in recent:
            outFile.write(str(recent[lineage]) + "\n")
        else:
            outFile.write("0\n")

    outFile.close()