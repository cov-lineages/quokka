#Extracts the path to one or more sequences from an UShER sample-paths file

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help = "Sample-paths file")
    parser.add_argument("-s", nargs = "+", help = "One or more sample names to be extracted")
    args = parser.parse_args()

    #Sample names to be extracted
    sn = args.s

    #Iterate through the sample-paths, check if the line contains any of the samples of interest and print if so
    with open(args.p) as f:
        for l in f:
            if any(s in l for s in sn):
                print(l.strip())