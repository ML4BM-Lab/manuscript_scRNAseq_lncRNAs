# Index generation for each scRNA-seq preprocessing pipeline
hg38_gencode_genome_fasta="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38.primary_assembly_GENCODE.genome.fa"
hg38_gencode_gtf="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_gencode_transcriptome_fasta="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.transcripts.fa"

cd /home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/
# Cell Ranger
cellranger mkref --genome=GRCh38_CellRanger_GENCODE_ref --fasta=$hg38_gencode_genome_fasta  --genes=$hg38_gencode_gtf

#STARsolo
STAR  --runMode genomeGenerate  --genomeDir $genomeDir  --genomeFastaFiles $hg38_gencode_genome_fasta --sjdbGTFfile $hg38_gencode_gtf   --sjdbOverhang 100

# Kallisto
kallisto index -i transcriptome_index_kallisto.idx $hg38_gencode_transcriptome_fasta

# Salmon (selective alignment, following https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
cat $hg38_gencode_genome_fasta | grep ">" |  cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt   # prepare decoys
cat $hg38_gencode_transcriptome_fasta $hg38_gencode_genome_fasta > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -p 12 -i transcriptome_selective_alignment_index_salmon --gencode

