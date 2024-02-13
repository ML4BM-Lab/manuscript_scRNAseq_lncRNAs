# ScRNA-seq Extended benchmark of Human healthy intestine pool of epithelium cells from 8,9,13,20 PCW. GSM4808339. Reads are 150bp in R1 and R2.
# Downloaded fastq from GSM4808339
cd /home/egonie/kike/phd/test_data/GSE158702_intestine/

# 1. Cat lanes 
cat SRR12735687_1.fastq.gz SRR12735688_1.fastq.gz SRR12735689_1.fastq.gz SRR12735690_1.fastq.gz SRR12735691_1.fastq.gz SRR12735692_1.fastq.gz SRR12735693_1.fastq.gz SRR12735694_1.fastq.gz > ../FES_GSM4037320_S1_L001_R1_001.fastq.gz
cat SRR12735687_2.fastq.gz SRR12735688_2.fastq.gz SRR12735689_2.fastq.gz SRR12735690_2.fastq.gz SRR12735691_2.fastq.gz SRR12735692_2.fastq.gz SRR12735693_2.fastq.gz SRR12735694_2.fastq.gz > ../FES_GSM4037320_S1_L001_R2_001.fastq.gz

#2.1. Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/phd/test_data/GSE158702_intestine/01.CellRanger"
cellranger count --id=FES_GSM4037320 --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=FES_GSM4037320 --localcores=8 --localmem=56320

#2.2. Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome_150bp_onlyR2 -x 10xv3 -t 8 FES_GSM4037320_S1_L001_R1_001.fastq.gz FES_GSM4037320_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome_150bp_onlyR2
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -
