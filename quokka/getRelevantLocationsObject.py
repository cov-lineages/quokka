import Bio
import pandas as pd
import csv
import re
import pickle
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

## read in reference file
reference = SeqIO.read(sys.argv[1], "fasta")

## read in large file of designated sequences
sequences_index = Bio.SeqIO.index(sys.argv[2], "fasta")
## designation index

## load the csv of what lineage each sequence has been designated to
lineage_designations = pd.read_csv(sys.argv[3], ",")

lineages = set()
final_positions = set()

for index, row in lineage_designations.iterrows(): 
    lineages.add(row["lineage"])


## outputs a dataframe of mutations associated with a given lineage

## lineage of interest given as a string i.e. lineage_of interest = "B.1.1.7"
def lineage_associated_SNPs(designation_list, lineage_of_interest, designation_index, reference_sequence): 
    sequence_ID_list = []
    for index, row in designation_list.iterrows(): 
        if row["lineage"] == str(lineage_of_interest):
            sequence_ID_list.append(row["sequence_name"])
            
    ### output of list of sequence_IDs for the given lineage
    
    sequences_lineage = Bio.Align.MultipleSeqAlignment([])
    
    for sequence in designation_index: ## for a sequence in the indexed sequences
        if sequence in sequence_ID_list: ## is this sequence in the ID list? 
            sequences_lineage.append(designation_index[sequence]) ## if it is, add the sequence to the list
    
    
    info = SummaryInfo(sequences_lineage)
    consensus_sequence =  info.gap_consensus(
    threshold=0.95, 
    ambiguous='N')
        
    nuc_position = []
    reference_nuc = []
    lineage_nuc = []

    for i in range(len(consensus_sequence)):
        if reference_sequence[i] != "N" and consensus_sequence[i] != "N":
            if consensus_sequence[i] != reference_sequence[i]:
                nuc_position.append(i)
                reference_nuc.append(reference_sequence[i])
                lineage_nuc.append(consensus_sequence[i])
    
    defining_mutations = pd.DataFrame(
    {'nucleotide_position': nuc_position,
     'reference_nucleotide': reference_nuc,
     'lineage_nucleotide': lineage_nuc
    })
    
    return nuc_position



for l in lineages:
    print(l)
    for i in lineage_associated_SNPs(lineage_designations, l, sequences_index, reference):
        final_positions.add(i)
        print(i)
    print()

out_loc = os.path.join(sys.argv[4],"relevantPositions.pickle")

with open(out_loc, 'wb') as pickle_file:
    pickle.dump(final_positions, pickle_file)


def get_lineage_sequences(designation_list, lineage_of_interest, designation_index, reference_sequence): 
    sequence_ID_list = []
    for index, row in designation_list.iterrows(): 
        if row["lineage"] == str(lineage_of_interest):
            sequence_ID_list.append(row["taxon"])
            
    ### output of list of sequence_IDs for the given lineage
    
    sequences_lineage =  Bio.Align.MultipleSeqAlignment([])
    
    for sequence in designation_index: ## for a sequence in the indexed sequences
        if sequence in sequence_ID_list: ## is this sequence in the ID list? 
            sequences_lineage.append(designation_index[sequence]) ## if it is, add the sequence to the MultipleSeqAlignment object 
    return(lineage_sequences)


def mutation_count(lineage_sequences, defining_mutations):
    mutcount_list = []
    epi_ID = []

    for sequence in lineage_sequences: 
        epi_ID.append(sequence.name)

        mutation_count = 0

        for index, row in defining_mutations.iterrows():
            if row["lineage_nucleotide"] == sequence.seq[row['nucleotide_position']]:
                mutation_count  += 1

        mutcount_list.append(mutation_count)    
        
    mutation_counts = pd.DataFrame(
    {'Sequence_name': epi_ID,
     'mutation_count': mutcount_list})
    
    return(mutation_counts)
