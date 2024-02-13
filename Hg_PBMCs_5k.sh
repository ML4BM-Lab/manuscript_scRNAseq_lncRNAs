# ScRNA-seq Extended benchmark of healthy human 5k PBMCs
# Downloaded from https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-3-1-standard-3-0-2
cd /home/egonie/kike/phd/test_data/10X/PBMCS_5K_healthy_donor
wget https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-3-1-standard-3-0-2/5k_pbmc_v3_nextgem_fastqs.tar
tar -xvf 5k_pbmc_v3_nextgem_fastqs.tar

# 1. Cat lanes 
cd 5k_pbmc_v3_nextgem_fastqs
cat *I1* > ../5k_pbmc_v3_nextgem_alldata_S1_L001_I1_001.fastq.gz 
cat *R1*fastq.gz > ../5k_pbmc_v3_nextgem_alldata_S1_L001_R1_001.fastq.gz
cat *R2*fastq.gz > ../5k_pbmc_v3_nextgem_alldata_S1_L001_R2_001.fastq.gz

#2.1. Cell Ranger
REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/phd/test_data/10X/PBMCS_5K_healthy_donor"
cellranger count --id=5k_pbmc_v3_nextgem_alldata --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=5000 --sample=5k_pbmc_v3_nextgem_alldata --localcores=8 --localmem=56320

#2.2. Kallisto
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcriptome_index_kallisto.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/3M-february-2018.txt"
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/tr2g.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome -x 10xv3 -t 8 5k_pbmc_v3_nextgem_alldata_S1_L001_R1_001.fastq.gz 5k_pbmc_v3_nextgem_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
cd output_bus_total_transcriptome
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -




