from Bio import SeqIO

fasta_file = "/home/egonie/data/egonie/phd/reference_genomes_kike/GRCh38/gencode/intronic_sequences.fa"
sequences = []

# Intronic sequences
for record in SeqIO.parse(fasta_file, "fasta"):
    sequences.append(record)

from tqdm import tqdm
def get_kmers_from_seqrecord(seq_records, k, kmers_set):
    for record in tqdm(seq_records):
        sequence = str(record.seq)
        for i in range(len(sequence) - k + 1):
            kmers_set.add(sequence[i:i+k])

# Import exclusive/common lncRNAs/protein-coding genes. Read fasta file
for record in SeqIO.parse('sequences_exclusive_lncRNAs_hg_10k_PBMCs.fa', "fasta"):
    exclusive_lncRNAs.append(record)
for record in SeqIO.parse('sequences_common_lncRNAs_hg_10k_PBMCs.fa', "fasta"):
    common_lncRNAs.append(record)
for record in SeqIO.parse('sequences_exclusive_PCs_hg_10k_PBMCs.fa', "fasta"):
    exclusive_PCs.append(record)
for record in SeqIO.parse('sequences_common_PCs_hg_10k_PBMCs.fa', "fasta"):
    common_PCs.append(record)

############################################################################# k-mer overlap ###############################################################
#1.  k = 31
k = 31
kmers_set = set()
get_kmers_from_seqrecord(sequences, k, kmers_set)

exclusive_lncRNAs = []
common_lncRNAs = []
exclusive_PCs = []
common_PCs = []
exclusive_lncRNAs_kmers = set()
get_kmers_from_seqrecord(exclusive_lncRNAs, k, exclusive_lncRNAs_kmers)
common_lncRNAs_kmers = set()
get_kmers_from_seqrecord(common_lncRNAs, k, common_lncRNAs_kmers)
exclusive_PCs_kmers = set()
get_kmers_from_seqrecord(exclusive_PCs, k, exclusive_PCs_kmers)
common_PCs_kmers = set()
get_kmers_from_seqrecord(common_PCs, k, common_PCs_kmers)

#lncRNAs
len(kmers_set.intersection(common_lncRNAs_kmers)) # 66% kmers intersect
len(kmers_set.intersection(exclusive_lncRNAs_kmers))# 58% kmers intersect

# protein-coding genes
len(kmers_set.intersection(common_PCs_kmers)) # 44.63% kmers intersect
len(kmers_set.intersection(exclusive_PCs_kmers))# 45.63% kmers intersect

#2. k = 50
k = 50
kmers_set = set()
get_kmers_from_seqrecord(sequences, k, kmers_set)

exclusive_lncRNAs = []
common_lncRNAs = []
exclusive_PCs = []
common_PCs = []
exclusive_lncRNAs_kmers = set()
get_kmers_from_seqrecord(exclusive_lncRNAs, k, exclusive_lncRNAs_kmers)
common_lncRNAs_kmers = set()
get_kmers_from_seqrecord(common_lncRNAs, k, common_lncRNAs_kmers)
exclusive_PCs_kmers = set()
get_kmers_from_seqrecord(exclusive_PCs, k, exclusive_PCs_kmers)
common_PCs_kmers = set()
get_kmers_from_seqrecord(common_PCs, k, common_PCs_kmers)

#lncRNAs
len(kmers_set.intersection(common_lncRNAs_kmers)) # 63.07% kmers intersect
len(kmers_set.intersection(exclusive_lncRNAs_kmers))# 53.77% kmers intersect

# protein-coding genes
len(kmers_set.intersection(common_PCs_kmers)) # 41.75% kmers intersect
len(kmers_set.intersection(exclusive_PCs_kmers)) #41.28% kmers intersect

#3. k = 75
k = 75
kmers_set = set()
get_kmers_from_seqrecord(sequences, k, kmers_set)

exclusive_lncRNAs = []
common_lncRNAs = []
exclusive_PCs = []
common_PCs = []
exclusive_lncRNAs_kmers = set()
get_kmers_from_seqrecord(exclusive_lncRNAs, k, exclusive_lncRNAs_kmers)
common_lncRNAs_kmers = set()
get_kmers_from_seqrecord(common_lncRNAs, k, common_lncRNAs_kmers)
exclusive_PCs_kmers = set()
get_kmers_from_seqrecord(exclusive_PCs, k, exclusive_PCs_kmers)
common_PCs_kmers = set()
get_kmers_from_seqrecord(common_PCs, k, common_PCs_kmers)

#lncRNAs
len(kmers_set.intersection(common_lncRNAs_kmers)) # 60.05% kmers intersect
len(kmers_set.intersection(exclusive_lncRNAs_kmers))# 50.29% kmers intersect

# protein-coding genes
len(kmers_set.intersection(common_PCs_kmers)) # 38.48% kmers intersect
len(kmers_set.intersection(exclusive_PCs_kmers)) #38.19% kmers intersect

#4. k = 91
k = 91
kmers_set = set()
get_kmers_from_seqrecord(sequences, k, kmers_set)

exclusive_lncRNAs = []
common_lncRNAs = []
exclusive_PCs = []
common_PCs = []
exclusive_lncRNAs_kmers = set()
get_kmers_from_seqrecord(exclusive_lncRNAs, k, exclusive_lncRNAs_kmers)
common_lncRNAs_kmers = set()
get_kmers_from_seqrecord(common_lncRNAs, k, common_lncRNAs_kmers)
exclusive_PCs_kmers = set()
get_kmers_from_seqrecord(exclusive_PCs, k, exclusive_PCs_kmers)
common_PCs_kmers = set()
get_kmers_from_seqrecord(common_PCs, k, common_PCs_kmers)

#lncRNAs
len(kmers_set.intersection(common_lncRNAs_kmers)) # 58.50% kmers intersect
len(kmers_set.intersection(exclusive_lncRNAs_kmers))# 48.67% kmers intersect

# protein-coding genes
len(kmers_set.intersection(common_PCs_kmers)) # 36.72% kmers intersect
len(kmers_set.intersection(exclusive_PCs_kmers)) #36.61% kmers intersect

