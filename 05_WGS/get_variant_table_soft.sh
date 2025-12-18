#!/bin/bash
#SBATCH --account=def-rshap
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G               # memory per node
#SBATCH --time=0-00:30
#SBATCH --output=logs/slurm-%j.out

export JAVA_TOOL_OPTIONS="-Xms256m -Xmx6g"

module load picard
module load gatk
module load bwa
module load samtools
module load bcftools


gatk --java-options "-Xmx4g" VariantsToTable \
     -V ./vcf/frs_1_soft_filtered.nonbinary.vcf \
     -F CHROM -F POS -F FILTER -F QUAL -F TYPE -F QD -F AF -GF AD --show-filtered \
     -O frs_1_soft.table
     
     
gatk --java-options "-Xmx4g" VariantsToTable \
     -V ./vcf/frs_242_soft_filtered.nonbinary.vcf \
     -F CHROM -F POS -F FILTER -F QUAL -F TYPE -F QD -F AF -GF AD --show-filtered \
     -O frs_242_soft.table
     
     
gatk --java-options "-Xmx4g" VariantsToTable \
     -V ./vcf/frs_246_soft_filtered.nonbinary.vcf \
     -F CHROM -F POS -F FILTER -F QUAL -F TYPE -F QD -F AF -GF AD --show-filtered \
     -O frs_246_soft.table
     
     