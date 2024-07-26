# Download auxiliary files for running scRNA-seq preprocessing pipelines
cd /home/egonie/kike/phd/test_data/auxiliary_data/
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

cd ~/dato-activo/reference.genomes_kike/GRCh38/gencode/
cat gencode.v37.annotation.gtf | sed -n -e 's/^.*gene_id //p' | sed -n -e 's/gene_type.*//p' | sed 's/;//g' | sed 's/transcript_id //g'  | sed 's/"//g' | sed -e 's/ /\t/g' | awk '{print $2 "\t" $1}'  | grep '^ENSMUST' | uniq > tr2g.tsv

# for the intronic annotation (in R)
introns=read.table("~/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_introns_t2c.txt")
exons=read.table("~/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_introns_cdna_t2c.txt")
gtf=rtracklayer::import("~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf")
gene_name_exons_match <- gtf$gene_id[match(exons$V1, gtf$transcript_id)]
gene_name_introns_match <- gtf$gene_id[match(gsub("-.*","",introns$V1), gtf$transcript_id)]
new_tr2g <- cbind( c(exons$V1, introns$V1), c(gene_name_exons_match, gene_name_introns_match) )
write.table(new_tr2g,"~/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_introns_t2g_final.txt",quote = F, col.names = F, row.names = F, sep = "\t")

cd ~/dato-activo/reference.genomes_kike/GRCm39/gencode/
cat gencode.v37.annotation.gtf | sed -n -e 's/^.*gene_id //p' | sed -n -e 's/gene_type.*//p' | sed 's/;//g' | sed 's/transcript_id //g'  | sed 's/"//g' | sed -e 's/ /\t/g' | awk '{print $2 "\t" $1}'  | grep '^ENSM' | uniq > tr2g.tsv

# Analyze repeat content (RepeatMasker v4.1.5)
RepeatMasker_human="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38.primary_assembly_GENCODE.genome.fa.out_cleaned.gff"
./RepeatMasker /home/egonie/data/egonie/phd/reference_genomes_kike/GRCh38/gencode/GRCh38.primary_assembly_GENCODE.genome.fa -species human -gff -html -dir human
cat GRCh38.primary_assembly_GENCODE.genome.fa.out.gff | grep "chr" | grep -v ":(" | grep -v "gff-version" | grep -v "A-rich" | grep -v "G-rich" > $name

RepeatMasker_mouse="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/GRCm39.primary_assembly.genome.fa.out_cleaned.gff"
./RepeatMasker /home/egonie/data/egonie/phd/reference_genomes_kike/GRCm39/gencode/GRCm39.primary_assembly.genome.fa -species mouse -gff -html -dir mouse
cat GRCm39.primary_assembly.genome.fa.out.gff | grep "chr" | grep -v ":(" | grep -v "gff-version" | grep -v "A-rich" | grep -v "G-rich" > $RepeatMasker_mouse


# Seekr for grouping functionally-related lncRNAs by k-mer content (following  https://github.com/CalabreseLab/seekr)
# Installation:
pip install seekr

# human
cd /home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode
seekr_download_gencode lncRNA -s human -r 37
seekr_canonical_gencode v37_lncRNA.fa v37_lncRNA_canonical.fa
# Add ENST00000441208.1: It is the isoform AL121895.1-204 (2701 bp) from AL121895.1 
cat v37_lncRNA.fa | sed -n '/ENST00000441208/,/>/p' | sed '$d' >> v37_lncRNA_canonical.fa
seekr_norm_vectors v37_lncRNA_canonical.fa 
seekr_kmer_counts v37_lncRNA_canonical.fa -o 6mers.csv -mv mean.npy -sv std.npy  # the 6mers.csv is the input used for clustering AL121895.1 together with proven cis-repressors, cis-activators
seekr_pearson 6mers.csv 6mers.csv -o correlations_all_lncRNAs_kmers.csv 
seekr_graph correlations_all_lncRNAs_kmers.csv 0.13 -g correlations_communities_6mers.gml -c SEEKR_communities_6mers.csv 

# mouse
cd /home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode
seekr_download_gencode lncRNA -s mouse -r M27
seekr_canonical_gencode vM27_lncRNA.fa vM27_lncRNA_canonical.fa
seekr_norm_vectors vM27_lncRNA_canonical.fa
seekr_kmer_counts vM27_lncRNA_canonical.fa -o 6mers_mouse.csv -mv mean.npy -sv std.npy
seekr_pearson 6mers_mouse.csv 6mers_mouse.csv -o correlations_all_lncRNAs_mouse_kmers.csv 
seekr_graph correlations_all_lncRNAs_mouse_kmers.csv 0.13 -g correlations_all_lncRNAs_mouse_kmers.gml -c SEEKR_communities_6mers.csv  # CREO QUE AHORA HA SALIDO! Igualmente hago el graph-based clustering


######################################################################################################################################
##############################################  Parse NONCODE annotation  ############################################################
######################################################################################################################################
# Part of the parsing is indicated in auxiliary_data.r while other is in auxiliary_data.sh
cd /home/egonie/dato-activo/reference.genomes_kike/GRCh38/NONCODE
# Download from http://v5.noncode.org/download.php
# I need to add the geneid to the gtf (in R)
# remove all that have strand "."
cat NONCODEv5_human_hg38_lncRNA_complete.gtf | awk '$7 == "." { next } { print }' > NONCODEv5_human_hg38_lncRNA_complete_cleaned.gtf

# 2.1. Generate index
# 2.1.1. Cell Ranger
cellranger mkref --genome=GRCh38_CellRanger_NONCODE_complete_cleaned_ref --fasta=/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38.primary_assembly_GENCODE.genome.fa  --genes=/home/egonie/dato-activo/reference.genomes_kike/GRCh38/NONCODE/NONCODEv5_human_hg38_lncRNA_complete_cleaned.gtf

#2.1.2. Kallisto
kallisto index -i transcriptome_index_kallisto_NONCODE.idx NONCODEv5_human.fa

#2.2 Preprocess with Cell Ranger and Kallisto
cd /home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/NONCODE
#2.2.1. Cell Ranger
cd 01.CellRanger
cellranger count --id=Parent_NGSC3_DI_PBMC --transcriptome=/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/NONCODE/GRCh38_CellRanger_NONCODE_complete_cleaned_ref  --fastqs=/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/NONCODE/01.CellRanger --expect-cells=10000 --sample=Parent_NGSC3_DI_PBMC --localcores=8 --localmem=56320

#2.2.2. Kallisto
cd 01.Kallisto
kallisto bus -i /home/egonie/dato-activo/reference.genomes_kike/GRCh38/NONCODE/transcriptome_index_kallisto_NONCODE.idx -o output_bus_transcriptome -x 10xv3 -t 8 Parent_NGSC3_DI_PBMC_R1_alldata.fastq.gz Parent_NGSC3_DI_PBMC_R2_alldata.fastq.gz
# transcript to genes
cat ~/dato-activo/reference.genomes_kike/GRCh38/NONCODE/NONCODEv5_human_hg38_lncRNA_complete.gtf | sed -n -e 's/^.*gene_id //p' | sed -n -e 's/FPKM.*//p' | sed 's/;//g' | sed 's/transcript_id //g'  | sed 's/"//g' | sed -e 's/ /\t/g' | awk '{print $2 "\t" $1}' | uniq > tr2g.tsv
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g ../tr2g.tsv -e ../matrix.ec -t ../transcripts.txt --genecounts -


