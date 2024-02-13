# ScRNA-seq Extended benchmark of Human primary pulmonar arterial cells from a healthy donor. GSM5020383. 
# Downloaded fastq from GSM4808348
cd /home/egonie/kike/phd/test_data/GSE164829_healthy_lung/GSM5020383

# 1. Cat lanes 
cat SRR13436369_1.fastq.gz SRR13436370_1.fastq.gz > FES_GSM5020383_S1_L001_R1_001.fastq.gz
cat SRR13436369_2.fastq.gz SRR13436370_2.fastq.gz > FES_GSM5020383_S1_L001_R2_001.fastq.gz

#2.1. Cell Ranger
cd /home/egonie/kike/phd/test_data/GSE164829_healthy_lung/01.CellRanger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/phd/test_data/GSE164829_healthy_lung/01.CellRanger"
cellranger count --id=FES_GSM5020383 --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=FES_GSM5020383 --localcores=8 --localmem=56320

#2.2. Kallisto
cd /home/egonie/kike/phd/test_data/GSE164829_healthy_lung/01.Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 FES_GSM5020383_S1_L001_R1_001.fastq.gz FES_GSM5020383_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_transcriptome
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -
