# ScRNA-seq Extended benchmark of Human healthy intestine pool of stromal and epithelium cells from 12,17 and 19 PCW. GSM4808348. Reads are 150bp in R1 and R2.
# Downloaded fastq from GSM4808348
cd /home/egonie/kike/phd/test_data/GSE158702_intestine/GSM4808348

# 1. Cat lanes 
cat SRR12735757_1.fastq.gz SRR12735758_1.fastq.gz SRR12735759_1.fastq.gz SRR12735760_1.fastq.gz > FES_GSM4808348_S1_L001_R1_001.fastq.gz
cat SRR12735757_2.fastq.gz SRR12735758_2.fastq.gz SRR12735759_2.fastq.gz SRR12735760_2.fastq.gz > FES_GSM4808348_S1_L001_R2_001.fastq.gz

#2.1. Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/phd/test_data/GSE158702_intestine/01.CellRanger"
cellranger count --id=FES_GSM4808348 --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=FES_GSM4808348 --localcores=8 --localmem=56320

#2.2. Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome_150bp_onlyR2 -x 10xv3 -t 8 FES_GSM4808348_S1_L001_R1_001.fastq.gz FES_GSM4808348_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome_91bp_reads
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_transcriptome_91bp_reads
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -
