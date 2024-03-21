# TNBC scRNA-seq analysis from 5 samples. Raw fastq files need to be downloaded. They are available in GEO: GSE246142 (token provided in the manuscript)
cd /home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/

##############################################################################################################################################
# Sample1. AA017
cd NextSeq2000.RUN62.20220309

# Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN62.20220309/"
cellranger count --id=AA017_alldata --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=AA017_alldata --localcores=8 --localmem=56320

# Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 AA017_alldata_S1_L001_R1_001.fastq.gz AA017_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results_no_multimappers
cd bustools_results_no_multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -


##############################################################################################################################################
# Sample2. AA024
cd /home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN74.20220512

# Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN62.20220309/"
cellranger count --id=AA024_alldata --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=AA024_alldata --localcores=8 --localmem=56320

# Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 AA024_alldata_S1_L001_R1_001.fastq.gz AA024_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results_no_multimappers
cd bustools_results_no_multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

##############################################################################################################################################
# Sample3. AA025
cd /home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN81.aacacho

# Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN62.20220309/"
cellranger count --id=AA025_alldata --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=AA025_alldata --localcores=8 --localmem=56320

# Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 AA025_alldata_S1_L001_R1_001.fastq.gz AA025_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results_no_multimappers
cd bustools_results_no_multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -


##############################################################################################################################################
# Sample4. AA038
cd /home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN109.aacacho

# Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN62.20220309/"
cellranger count --id=AA038_alldata --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=AA038_alldata --localcores=8 --localmem=56320

# Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 AA038_alldata_S1_L001_R1_001.fastq.gz AA038_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results_no_multimappers
cd bustools_results_no_multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -


##############################################################################################################################################
# Sample5. AA051
cd /home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN109.aacacho

# Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN62.20220309/"
cellranger count --id=AA051_alldata --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=AA051_alldata --localcores=8 --localmem=56320

# Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 AA051_alldata_S1_L001_R1_001.fastq.gz AA051_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results_no_multimappers
cd bustools_results_no_multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -