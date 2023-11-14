# ScRNA-seq benchmark of healthy human 10k PBMCs
# Downloaded from https://www.10xgenomics.com/resources/datasets/pbm-cs-from-a-healthy-donor-whole-transcriptome-analysis-3-1-standard-4-0-0
cd /home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/
# 1. Cat lanes 
cat Parent_NGSC3_DI_PBMC_S1_L001_I1_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L002_I1_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L003_I1_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L004_I1_001.fastq.gz > Parent_NGSC3_DI_PBMC_I1_alldata.fastq.gz
cat Parent_NGSC3_DI_PBMC_S1_L001_I2_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L002_I2_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L003_I2_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L004_I2_001.fastq.gz > Parent_NGSC3_DI_PBMC_I2_alldata.fastq.gz
cat Parent_NGSC3_DI_PBMC_S1_L001_R1_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L002_R1_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L003_R1_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L004_R1_001.fastq.gz > Parent_NGSC3_DI_PBMC_R1_alldata.fastq.gz
cat Parent_NGSC3_DI_PBMC_S1_L001_R2_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L002_R2_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L003_R2_001.fastq.gz Parent_NGSC3_DI_PBMC_S1_L004_R2_001.fastq.gz > Parent_NGSC3_DI_PBMC_R2_alldata.fastq.gz

#2.1. Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/"
cellranger count --id=Parent_NGSC3_DI_PBMC --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=Parent_NGSC3_DI_PBMC --localcores=8 --localmem=56320

#2.2. STARsolo
STAR_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/star_indices"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
STAR --genomeDir $STAR_INDEX --readFilesCommand zcat --readFilesIn Parent_NGSC3_DI_PBMC_R2_alldata.fastq.gz Parent_NGSC3_DI_PBMC_R1_alldata.fastq.gz --soloUMIfiltering MultiGeneUMI_CR --outSAMtype BAM SortedByCoordinate --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --outFilterScoreMin 30 --soloCellFilter  EmptyDrops_CR  --soloType CB_UMI_Simple --soloCBwhitelist $TENxV3_WHITELIST --soloUMIlen 12 --limitOutSJcollapsed 2000000

#2.3. Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 Parent_NGSC3_DI_PBMC_R1_alldata.fastq.gz Parent_NGSC3_DI_PBMC_R2_alldata.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

#2.4. Salmon
SALMON_INDEX_SELECTIVE_ALIGNMENT="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_selective_alignment_index_salmon"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
salmon alevin -lISR -1 Parent_NGSC3_DI_PBMC_R1_alldata.fastq.gz -2 Parent_NGSC3_DI_PBMC_R2_alldata.fastq.gz --chromiumV3 -i $SALMON_INDEX_SELECTIVE_ALIGNMENT -p 8 -o alevin_selective_alignment_output --tgMap $TRANSCRIPTS_TO_GENES



