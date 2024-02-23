# Single-cell multiome 3K PBMCs from a healthy female donor
# Downloaded from: https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0
cd /home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/pbmc_granulocyte_sorted_3k

# 1. scATAC-seq directly downloaded from: https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi

# 2. Download scRNA-seq fastq files
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_fastqs.tar
tar -xvf pbmc_granulocyte_sorted_3k_fastqs.tar
cd pbmc_granulocyte_sorted_3k/gex
cat *I1* > ../pbmc_granulocyte_sorted_3k_alldata_S1_L001_I1_001.fastq.gz
cat *I2* > ../pbmc_granulocyte_sorted_3k_alldata_S1_L001_I2_001.fastq.gz 
cat *R1*fastq.gz > ../pbmc_granulocyte_sorted_3k_alldata_S1_L001_R1_001.fastq.gz
cat *R2*fastq.gz > ../pbmc_granulocyte_sorted_3k_alldata_S1_L001_R2_001.fastq.gz

#2.1. Cell Ranger
# following to analyze only gene expression https://kb.10xgenomics.com/hc/en-us/articles/360059656912-Can-I-analyze-only-the-Gene-Expression-data-from-my-single-cell-multiome-experiment-
# version 5(include introns and chemistry = --chemistry=ARC-v1)
cd /home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.CellRanger
ln -s ../pbmc_granulocyte_sorted_3k/*S1*fastq.gz .

REF_TRANSCRIPTOME_CELLRANGER="/datos/huartelab/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38_CellRanger_GENCODE_ref"
FASTQ_DIR="/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.CellRanger"
cellranger count --id=pbmc_granulocyte_sorted_3k_alldata_cellRanger_modified_GEX_test2 --transcriptome=$REF_TRANSCRIPTOME_CELLRANGER  --fastqs=$FASTQ_DIR --expect-cells=3000 --sample=pbmc_granulocyte_sorted_3k_alldata --include-introns --chemistry=ARC-v1 --localcores=8 --localmem=56320

#2.2 Kallisto
# First it is needed to build the reference counting also introns (generate the index according to: https://www.biostars.org/p/468180/)
cd /home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode
pip install kb-python
kb ref -i kallisto_nuclei_introns_index.idx -g kallisto_nuclei_introns_t2g.txt -f1 kallisto_nuclei_introns_cdna.fa -f2 kallisto_nuclei_introns_intron.fa -c1 kallisto_nuclei_introns_cdna_t2c.txt -c2 kallisto_nuclei_introns_t2c.txt --workflow lamanno -n 8 GRCh38.primary_assembly_GENCODE.genome.fa gencode.v37.annotation.gtf
cat kallisto_nuclei_introns_cdna.fa kallisto_nuclei_introns_intron.fa > kallisto_nuclei_cDNA_introns_ALL.fa
kallisto index -i kallisto_nuclei_cDNA_introns.idx kallisto_nuclei_cDNA_introns_ALL.fa

cd /home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.Kallisto
ln -s ../pbmc_granulocyte_sorted_3k/*S1*fastq.gz .

KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_cDNA_introns.idx"
TENxV3_WHITELIST="/home/egonie/kike/phd/test_data/auxiliary_data/737K-arc-v1.txt" # Different whitelist (https://www.singlecellcourse.org/processing-raw-scrna-seq-sequencing-data-from-reads-to-a-count-matrix.html)
TRANSCRIPTS_TO_GENES="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_introns_t2g_final.txt" #Different t2g files (now with introns)
kallisto bus -i $KALLISTO_INDEX -o output_bus_transcriptome_introns -x 10xv3 -t 16 pbmc_granulocyte_sorted_3k_alldata_S1_L001_R1_001.fastq.gz pbmc_granulocyte_sorted_3k_alldata_S1_L001_R2_001.fastq.gz
cd output_bus_transcriptome_introns
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
mkdir bustools_results_no_multimappers
cd bustools_results_no_multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -




