# ScRNA-seq Extended benchmark of healthy mouse 10k PBMCs
# Downloaded from https://support.10xgenomics.com/single-cell-gene-expression/datasets/6.0.0/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_Multiplex
cd /home/egonie/kike/phd/test_data/10X/Mouse_PBMCs_Multiplexed_fastqs/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_1_gex
wget https://cg.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_Multiplex_fastqs.tar

#1.1. Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCm39/gencode/CELLRANGER_REF"
FASTQ_DIR="/home/egonie/kike/phd/test_data/10X/Mouse_PBMCs_Multiplexed_fastqs/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_1_gex/03.CellRanger/"
cellranger count --id=Mouse_PBMC_10K_gex --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=10000 --sample=Mouse_PBMC_10K_gex --localcores=8 --localmem=56320

#1.2. Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_1_gex_S2_L002_R1_001.fastq.gz SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_1_gex_S2_L002_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_transcriptome
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

