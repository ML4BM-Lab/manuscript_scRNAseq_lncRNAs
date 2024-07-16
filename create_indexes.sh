# Index generation for each scRNA-seq preprocessing pipeline
# Hg38
hg38_gencode_genome_fasta="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38.primary_assembly_GENCODE.genome.fa"
hg38_gencode_gtf="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_gencode_transcriptome_fasta="/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.transcripts.fa"

cd /home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/
# Cell Ranger
cellranger mkref --genome=GRCh38_CellRanger_GENCODE_ref --fasta=$hg38_gencode_genome_fasta  --genes=$hg38_gencode_gtf

#STARsolo
STAR  --runMode genomeGenerate  --genomeDir star_indices  --genomeFastaFiles $hg38_gencode_genome_fasta --sjdbGTFfile $hg38_gencode_gtf   --sjdbOverhang 100

# Kallisto
kallisto index -i transcriptome_index_kallisto.idx $hg38_gencode_transcriptome_fasta

# Salmon (selective alignment, following https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
cat $hg38_gencode_genome_fasta | grep ">" |  cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt   # prepare decoys
cat $hg38_gencode_transcriptome_fasta $hg38_gencode_genome_fasta > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -p 12 -i transcriptome_selective_alignment_index_salmon --gencode


# Index generation for single-cell multiome (include intronic reads, see: https://www.biostars.org/p/468180/)
pip install kb-python
kb ref -i kallisto_nuclei_introns_index.idx -g kallisto_nuclei_introns_t2g.txt -f1 kallisto_nuclei_introns_cdna.fa -f2 kallisto_nuclei_introns_intron.fa -c1 kallisto_nuclei_introns_cdna_t2c.txt -c2 kallisto_nuclei_introns_t2c.txt --workflow lamanno -n 8 $hg38_gencode_genome_fasta $hg38_gencode_gtf
cat kallisto_nuclei_introns_cdna.fa kallisto_nuclei_introns_intron.fa > kallisto_nuclei_cDNA_introns_ALL.fa
kallisto index -i kallisto_nuclei_cDNA_introns.idx kallisto_nuclei_cDNA_introns_ALL.fa

# mm39
mm39_gencode_genome_fasta="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/GRCm39.primary_assembly.genome.fa"
mm39_gencode_gtf="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/gencode.vM27.annotation.gtf"
mm39_gencode_transcriptome_fasta="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/gencode.vM27.transcripts.fa"

cd /home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode
# Cell Ranger
cellranger mkref --genome=CELLRANGER_REF --fasta=$mm39_gencode_genome_fasta  --genes=$mm39_gencode_gtf

#STARsolo
STAR  --runMode genomeGenerate  --genomeDir STAR_index  --genomeFastaFiles $mm39_gencode_genome_fasta --sjdbGTFfile $mm39_gencode_gtf   --sjdbOverhang 100

# Kallisto
kallisto index -i transcriptome_index_kallisto.idx $mm39_gencode_transcriptome_fasta

# Salmon (selective alignment, following https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
cat $mm39_gencode_genome_fasta | grep ">" |  cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt   # prepare decoys
cat $mm39_gencode_transcriptome_fasta $mm39_gencode_genome_fasta > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -p 12 -i transcriptome_selective_alignment_index_salmon --gencode




