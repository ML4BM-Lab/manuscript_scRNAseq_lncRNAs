# Download auxiliary files for running scRNA-seq preprocessing pipelines
cd /home/egonie/kike/phd/test_data/auxiliary_data/
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

cd ~/dato-activo/reference.genomes_kike/GRCh38/gencode/
cat gencode.v37.annotation.gtf | sed -n -e 's/^.*gene_id //p' | sed -n -e 's/gene_type.*//p' | sed 's/;//g' | sed 's/transcript_id //g'  | sed 's/"//g' | sed -e 's/ /\t/g' | awk '{print $2 "\t" $1}'  | grep '^ENST' | uniq > tr2g.tsv
