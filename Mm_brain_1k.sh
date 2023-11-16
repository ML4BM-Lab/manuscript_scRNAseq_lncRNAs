# ScRNA-seq benchmark of 1k mouse brain cells
# Downloaded from https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0
cd /home/egonie/kike/phd/test_data/10X/Mouse_Brain_1K_cells
# 1. Cat lanes 
cat neuron_1k_v3_fastqs/neuron_1k_v3_S1_L001_I1_001.fastq.gz neuron_1k_v3_fastqs/neuron_1k_v3_S1_L002_I1_001.fastq.gz > neuron_1k_v3_alldata_S1_L001_I1_001.fastq.gz
cat neuron_1k_v3_fastqs/neuron_1k_v3_S1_L001_R1_001.fastq.gz neuron_1k_v3_fastqs/neuron_1k_v3_S1_L002_R1_001.fastq.gz > neuron_1k_v3_alldata_S1_L001_R1_001.fastq.gz
cat neuron_1k_v3_fastqs/neuron_1k_v3_S1_L001_R2_001.fastq.gz neuron_1k_v3_fastqs/neuron_1k_v3_S1_L002_R2_001.fastq.gz > neuron_1k_v3_alldata_S1_L001_R2_001.fastq.gz

#2.1. Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCm39/gencode/CELLRANGER_REF"
FASTQ_DIR="/home/egonie/kike/phd/test_data/10X/Mouse_Brain_1K_cells/01.CellRanger/"
cellranger count --id=neuron_1k_v3_alldata --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=1000 --sample=neuron_1k_v3_alldata --localcores=8 --localmem=56320

#2.2. STARsolo
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
STAR_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/STAR_index"
STAR --genomeDir $STAR_INDEX --readFilesCommand zcat --readFilesIn neuron_1k_v3_alldata_S1_L001_R2_001.fastq.gz neuron_1k_v3_alldata_S1_L001_R1_001.fastq.gz --soloUMIfiltering MultiGeneUMI_CR --outSAMtype BAM SortedByCoordinate --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --outFilterScoreMin 30 --soloCellFilter  EmptyDrops_CR  --soloType CB_UMI_Simple --soloCBwhitelist $TENxV3_WHITELIST --soloUMIlen 12

#2.3. Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 neuron_1k_v3_alldata_S1_L001_R1_001.fastq.gz neuron_1k_v3_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_transcriptome
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

#2.4. Salmon
SALMON_INDEX_SELECTIVE_ALIGNMENT="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcriptome_selective_alignment_index_salmon"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/tr2g.tsv"
salmon alevin -lISR -1 neuron_1k_v3_alldata_S1_L001_R1_001.fastq.gz -2 neuron_1k_v3_alldata_S1_L001_R2_001.fastq.gz --chromiumV3 -i $SALMON_INDEX_SELECTIVE_ALIGNMENT -p 8 -o alevin_selective_alignment_output --tgMap $TRANSCRIPTS_TO_GENES
