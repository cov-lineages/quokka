from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

def get_95_cns(seq_file,cns_threshold=0.7):
    alignment = AlignIO.read(seq_file, "fasta")
    info = SummaryInfo(alignment)
    cns =  info.gap_consensus(
    threshold=cns_threshold, 
    ambiguous='N')
    return cns

print(get_95_cns("sequences.aln.fasta"))