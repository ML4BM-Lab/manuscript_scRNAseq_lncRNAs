# Download auxiliary files for running scRNA-seq preprocessing pipelines
cd /home/egonie/kike/phd/test_data/auxiliary_data/
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

cd ~/dato-activo/reference.genomes_kike/GRCh38/gencode/
cat gencode.v37.annotation.gtf | sed -n -e 's/^.*gene_id //p' | sed -n -e 's/gene_type.*//p' | sed 's/;//g' | sed 's/transcript_id //g'  | sed 's/"//g' | sed -e 's/ /\t/g' | awk '{print $2 "\t" $1}'  | grep '^ENST' | uniq > tr2g.tsv

# for the intronic annotation (in R)
introns=read.table("~/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_introns_t2c.txt")
exons=read.table("~/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_introns_cdna_t2c.txt")
gtf=rtracklayer::import("~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf")
gene_name_exons_match <- gtf$gene_id[match(exons$V1, gtf$transcript_id)]
gene_name_introns_match <- gtf$gene_id[match(gsub("-.*","",introns$V1), gtf$transcript_id)]
new_tr2g <- cbind( c(exons$V1, introns$V1), c(gene_name_exons_match, gene_name_introns_match) )
write.table(new_tr2g,"~/dato-activo/reference.genomes_kike/GRCh38/gencode/kallisto_nuclei_introns_t2g_final.txt",quote = F, col.names = F, row.names = F, sep = "\t")

cd ~/dato-activo/reference.genomes_kike/GRCm39/gencode/
cat gencode.v37.annotation.gtf | sed -n -e 's/^.*gene_id //p' | sed -n -e 's/gene_type.*//p' | sed 's/;//g' | sed 's/transcript_id //g'  | sed 's/"//g' | sed -e 's/ /\t/g' | awk '{print $2 "\t" $1}'  | grep '^ENSM' | uniq > tr2g.tsv

