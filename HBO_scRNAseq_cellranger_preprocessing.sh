
##Cellranger sequence processing##

#download human reference genome
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz

#load cellranger
module load cellranger

#run count
cellranger count --id=run_uti89_count_hbo1_65194 \
   --fastqs=/mount/patras/clare/sc-seq/hbo/hbo1/raw-data \
   --sample=BLT-65194 \
   --transcriptome=/mount/patras/clare/sc-seq/GCF_000013265.1_ASM1326v1_rna_from_genomic.fna \
    --create-bam=true

cellranger count --id=run_uti89_count_hbo1_65195 \
   --fastqs=/mount/patras/clare/sc-seq/hbo/hbo1/raw-data \
   --sample=BLT-65195 \
   --transcriptome=/mount/patras/clare/sc-seq/GCF_000013265.1_ASM1326v1_rna_from_genomic.fna \
    --create-bam=true
