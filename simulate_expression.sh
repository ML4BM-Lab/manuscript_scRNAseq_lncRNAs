cd /home/egonie/kike/phd/test_data/Minju_ha_nat_comm/simulation_scRNAseq/
# Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/phd/test_data/Minju_ha_nat_comm/simulation_scRNAseq/"

cellranger count --id=simulation1_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu1 --localcores=8 --localmem=56320
cellranger count --id=simulation2_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu2 --localcores=8 --localmem=56320
cellranger count --id=simulation3_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu3 --localcores=8 --localmem=56320
cellranger count --id=simulation4_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu4 --localcores=8 --localmem=56320
cellranger count --id=simulation5_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu5 --localcores=8 --localmem=56320
cellranger count --id=simulation6_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu6 --localcores=8 --localmem=56320
cellranger count --id=simulation7_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu7 --localcores=8 --localmem=56320
cellranger count --id=simulation8_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu8 --localcores=8 --localmem=56320
cellranger count --id=simulation9_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu9 --localcores=8 --localmem=56320
cellranger count --id=simulation10_CellRanger --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --sample=simulationsimu10 --localcores=8 --localmem=56320

# Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"

kallisto bus -i $KALLISTO_INDEX -o simulation1_kallisto -x 10xv3 -t 8 simulationsimu1_S1_L001_R1_001.fastq.gz simulationsimu1_S1_L001_R2_001.fastq.gz
cd simulation1_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation2_kallisto -x 10xv3 -t 8 simulationsimu2_S1_L001_R1_001.fastq.gz simulationsimu2_S1_L001_R2_001.fastq.gz
cd simulation2_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation3_kallisto -x 10xv3 -t 8 simulationsimu3_S1_L001_R1_001.fastq.gz simulationsimu3_S1_L001_R2_001.fastq.gz
cd simulation3_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation4_kallisto -x 10xv3 -t 8 simulationsimu4_S1_L001_R1_001.fastq.gz simulationsimu4_S1_L001_R2_001.fastq.gz
cd simulation4_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -


kallisto bus -i $KALLISTO_INDEX -o simulation5_kallisto -x 10xv3 -t 8 simulationsimu5_S1_L001_R1_001.fastq.gz simulationsimu5_S1_L001_R2_001.fastq.gz
cd simulation5_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation6_kallisto -x 10xv3 -t 8 simulationsimu6_S1_L001_R1_001.fastq.gz simulationsimu6_S1_L001_R2_001.fastq.gz
cd simulation6_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation7_kallisto -x 10xv3 -t 8 simulationsimu7_S1_L001_R1_001.fastq.gz simulationsimu7_S1_L001_R2_001.fastq.gz
cd simulation7_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation8_kallisto -x 10xv3 -t 8 simulationsimu8_S1_L001_R1_001.fastq.gz simulationsimu8_S1_L001_R2_001.fastq.gz
cd simulation8_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation9_kallisto -x 10xv3 -t 8 simulationsimu9_S1_L001_R1_001.fastq.gz simulationsimu9_S1_L001_R2_001.fastq.gz
cd simulation9_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

kallisto bus -i $KALLISTO_INDEX -o simulation10_kallisto -x 10xv3 -t 8 simulationsimu10_S1_L001_R1_001.fastq.gz simulationsimu10_S1_L001_R2_001.fastq.gz
cd simulation10_kallisto
mv transcripts.txt transcripts_back.txt
cat transcripts_back.txt | sed "s/|.*//" > transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

