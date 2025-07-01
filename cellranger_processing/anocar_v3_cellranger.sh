#!/bin/bash

#SBATCH --account=aealmada_561
#SBATCH --partition=main
#SBATCH --job-name=cell_ranger
#SBATCH --mail-user=difeizhu@usc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=cell_ranger_slurm.out
#SBATCH --error=cell_ranger_slurm.err


# construct cell ranger reference
REF_DIR = "/project/aealmada_561/judyz/single_cell/cell_ranger/AnoCar_3_1"

if [ -d "$REF_DIR" ]; then
  echo "Reference genome already exists. Skipping mkref step."
else
  echo "Reference genome not found. Running cellranger mkref."
  # construct cell ranger reference
  cellranger mkref \
    --genome=AnoCar_3_1 \
    --fasta=/project/aealmada_561/judyz/single_cell/cell_ranger/reference/GCF_035594765.1.fa \
    --genes=/project/aealmada_561/judyz/single_cell/cell_ranger/reference/GCF_035594765.1_rAnoCar3.1.pri.ncbiRefSeq.gtf
fi


# run quads and tail cell ranger count pipeline
## quads
cellranger count --id=Quads \
           --transcriptome=/project/aealmada_561/judyz/single_cell/anocar_v3/cell_ranger/AnoCar_3_1 \
           --fastqs=/project/aealmada_561/judyz/single_cell/cell_ranger/fastqs \
           --sample=Quads \
           --create-bam=true \
           --localcores=8 \
           --localmem=64

## tail
cellranger count --id=Tail \
           --transcriptome=/project/aealmada_561/judyz/single_cell/anocar_v3/cell_ranger/AnoCar_3_1 \
           --fastqs=/project/aealmada_561/judyz/single_cell/cell_ranger/fastqs \
           --sample=Tail \
           --create-bam=true \
           --localcores=8 \
           --localmem=64