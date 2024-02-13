# ScRNA-seq Extended benchmark of Human pulmonary fibrosis sample. GSM4037320. v2 chemistry
# Downloaded fastq from GSM4037316
cd /home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis

# 1. Cat lanes 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/015/SRR10368215/SRR10368215_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/015/SRR10368215/SRR10368215_2.fastq.gz
mv VUILD63_SRR9985406_1.fastq.gz VUILD63_SRR9985406_S1_L001_R1_001.fastq.gz
mv VUILD63_SRR9985406_2.fastq.gz VUILD63_SRR9985406_S1_L001_R2_001.fastq.gz

#2.1. Cell Ranger
cd /home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis/03.CellRanger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis/03.CellRanger"
cellranger count --id=VUILD63 --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=7500 --sample=VUILD63_SRR9985406 --localcores=8 --localmem=56320

#2.2. Kallisto
cd /home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis/03.Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome_split_cDNA -x 0,0,16:0,16,26:0,26,0,1,0,0 -t 8 VUILD63_SRR9985406_S1_L001_R1_001.fastq.gz VUILD63_SRR9985406_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome_split_cDNA
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_transcriptome_split_cDNA
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -
