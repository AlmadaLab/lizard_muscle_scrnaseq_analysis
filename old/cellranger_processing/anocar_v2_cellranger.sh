#!/bin/bash

#SBATCH --account=aealmada_561
#SBATCH --partition=main
#SBATCH --job-name=cell_ranger_v2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=cell_ranger_v2_slurm.out
#SBATCH --error=cell_ranger_v_slurm.err

cellranger mkref \
  --genome=AnoCar_2_0 \
  --fasta=/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/reference/Anolis_carolinensis.AnoCar2.0v2.dna.complete.fa \
  --genes=/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/reference/Anolis_carolinensis.AnoCar2.0v2.113.gtf


cellranger count --id=Quads \
           --transcriptome=/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/AnoCar_2_0 \
           --fastqs=/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/fastqs \
           --sample=Quads \
           --create-bam=true \
           --localcores=8 \
           --localmem=64

cellranger count --id=Tail \
           --transcriptome=/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/AnoCar_2_0 \
           --fastqs=/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/fastqs \
           --sample=Tail \
           --create-bam=true \
           --localcores=8 \
           --localmem=64

