#!/bin/bash

#SBATCH --account=def-rshap
#SBATCH --cpus-per-task=2
#SBATCH --mem=90G               # memory per node
#SBATCH --time=0-24:30
#SBATCH --output=./logs/slurm-%j.out

module load samtools

samtools mpileup ./bam_align/frs_1_align.reordered.bam --positions samtools_het_pos_list.txt -o ./pileups/frs1.all_het.pileup --no-output-ins --no-output-del --no-output-ends -f GCF_000182965.3_ASM18296v3_genomic.fna
samtools mpileup ./bam_align/frs_242_align.reordered.bam --positions samtools_het_pos_list.txt -o ./pileups/frs242.all_het.pileup --no-output-ins --no-output-del --no-output-ends -f GCF_000182965.3_ASM18296v3_genomic.fna
samtools mpileup ./bam_align/frs_246_align.reordered.bam --positions samtools_het_pos_list.txt -o ./pileups/frs246.all_het.pileup --no-output-ins --no-output-del --no-output-ends -f GCF_000182965.3_ASM18296v3_genomic.fna

