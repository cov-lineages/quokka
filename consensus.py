from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

def get_95_cns(seq_file):
    alignment = AlignIO.read(seq_file, "fasta")
    info = SummaryInfo(alignment)
    cns =  info.gap_consensus(
    threshold=0.7, 
    ambiguous='N')
    return cns

print(get_95_cns("sequences.aln.fasta"))