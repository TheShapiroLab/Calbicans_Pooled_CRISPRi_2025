#!/bin/bash
#SBATCH --account=def-rshap
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G               # memory per node
#SBATCH --time=0-03:00
#SBATCH --output=logs/slurm-%j.out

export JAVA_TOOL_OPTIONS="-Xms256m -Xmx6g"

module load picard
module load gatk
module load bwa
module load samtools
module load bcftools

####################### FRS1

# convert fastq
java -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=./raw_data/fRS1_S175_R1_001.fastq FASTQ2=./raw_data/fRS1_S175_R2_001.fastq OUTPUT=./u_bam/frs_1_u.bam READ_GROUP_NAME=LH00162 SAMPLE_NAME=FRS1 LIBRARY_NAME=FRS1 PLATFORM_UNIT=LH00162 PLATFORM=illumina

#mark illumina adapters
gatk --java-options "-Xmx4g" MarkIlluminaAdapters -I ./u_bam/frs_1_u.bam -O ./u_bam_clip/frs_1_u_clipped.bam -M ./logs/frs1_u_markilluminaadapters_metrics.txt

# switch back to fastq with trimmed adapters
gatk --java-options "-Xmx4g" SamToFastq -I ./u_bam_clip/frs_1_u_clipped.bam --FASTQ ./fastq_clip/frs_1_u_clip.fastq  -CLIP_ATTR XT  -CLIP_ACT 2  -INTER true -NON_PF true

# align on ref with bwa mem
bwa mem -M -t 4 -p GCF_000182965.3_ASM18296v3_genomic.fna ./fastq_clip/frs_1_u_clip.fastq -o sam_align/frs_1_.sam

# map aligned fastq reads back to bam
# reference sequence must also have a sequence dict created using picard
gatk --java-options "-Xmx4g" MergeBamAlignment -ALIGNED sam_align/frs_1_.sam -UNMAPPED u_bam_clip/frs_1_u_clipped.bam -R GCF_000182965.3_ASM18296v3_genomic.fna -O ./bam_align/frs_1_align.bam

# mark duplicate reads
gatk --java-options "-Xmx4g" MarkDuplicates -I ./bam_align/frs_1_align.bam -O ./bam_align/frs_1_align.marked_duplicates.bam -M frs_1_.marked_duplicates.metrics

# reorder bam
# -SD is a hidden (in the docs) variable for the reference genome
gatk --java-options "-Xmx4g" ReorderSam -I ./bam_align/frs_1_align.marked_duplicates.bam -O ./bam_align/frs_1_align.reordered.bam -SD GCF_000182965.3_ASM18296v3_genomic.fna

# build bam index
gatk --java-options "-Xmx4g" BuildBamIndex -I ./bam_align/frs_1_align.reordered.bam -O ./bam_align/frs_1_align.reordered.bai

# call haplotypes
gatk --java-options "-Xmx4g" HaplotypeCaller -R GCF_000182965.3_ASM18296v3_genomic.fna -I ./bam_align/frs_1_align.reordered.bam -O ./vcf/frs_1.g.vcf.gz -ERC GVCF -ploidy 2
# takes about 20 minutes/WGS candida albicans sample

# after all samples are processed, merge vcfs
##gatk --java-options "-Xmx4G" CombineGVCFs -R ./reference_genome/data/GCA_002759435.2_Cand_auris_B8441_V2_genomic.fna -O ./vcf/variants.merged.gvcf --variant ./vcf/SRR3883436.g.vcf.gz

# Perform joint genotyping
gatk --java-options "-Xmx4G" GenotypeGVCFs -R GCF_000182965.3_ASM18296v3_genomic.fna -O ./vcf/frs_1.out.gvcf -V ./vcf/frs_1.g.vcf.gz

#Select variants, both SNPs and raw_indels. Mixed sites are discarted.
gatk --java-options "-Xmx4G" SelectVariants -V ./vcf/frs_1.out.gvcf  -R GCF_000182965.3_ASM18296v3_genomic.fna -select-type SNP -O ./vcf/frs_1.snps.vcf.gz

gatk --java-options "-Xmx4G" SelectVariants -V ./vcf/frs_1.out.gvcf  -R GCF_000182965.3_ASM18296v3_genomic.fna -select-type INDEL -O ./vcf/frs_1.indels.vcf.gz

#Filter variants, different paramters for Indels and SNPs. Recommendations by GATK.
gatk --java-options "-Xmx4G" VariantFiltration -V ./vcf/frs_1.snps.vcf.gz -R GCF_000182965.3_ASM18296v3_genomic.fna \
-filter "QD < 12.5" --filter-name "QD12.5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-O ./vcf/frs_1.soft.snps_filtered.vcf.gz

gatk --java-options "-Xmx4G" VariantFiltration -V ./vcf/frs_1.indels.vcf.gz -R GCF_000182965.3_ASM18296v3_genomic.fna \
-filter "QD < 20.0" --filter-name "QD20" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-O ./vcf/frs_1.soft.indels_filtered.vcf.gz

#Combine output vcfs into a single output. The output of VariantFiltration is sorted.
gatk --java-options "-Xmx4G" MergeVcfs -R GCF_000182965.3_ASM18296v3_genomic.fna \
-I ./vcf/frs_1.soft.snps_filtered.vcf.gz \
-I ./vcf/frs_1.soft.indels_filtered.vcf.gz \
-O ./vcf/frs_1_soft_filtered.vcf.gz


# make vcf human readable
bcftools view ./vcf/frs_1_soft_filtered.vcf.gz > ./vcf/frs_1_soft_filtered.nonbinary.vcf
bcftools view ./vcf/frs_1_soft_filtered.nonbinary.vcf -f "PASS" -o ./vcf/frs_1_soft_filtered.nonbinary.dropped.vcf


####################### FRS242
# convert fastq
java -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=./raw_data/242_S176_R1_001.fastq FASTQ2=./raw_data/242_S176_R2_001.fastq OUTPUT=./u_bam/frs_242_u.bam READ_GROUP_NAME=LH00162 SAMPLE_NAME=FRS242 LIBRARY_NAME=FRS242 PLATFORM_UNIT=LH00162 PLATFORM=illumina

#mark illumina adapters
gatk --java-options "-Xmx4g" MarkIlluminaAdapters -I ./u_bam/frs_242_u.bam -O ./u_bam_clip/frs_242_u_clipped.bam -M ./logs/frs242_u_markilluminaadapters_metrics.txt

# switch back to fastq with trimmed adapters
gatk --java-options "-Xmx4g" SamToFastq -I ./u_bam_clip/frs_242_u_clipped.bam --FASTQ ./fastq_clip/frs_242_u_clip.fastq  -CLIP_ATTR XT  -CLIP_ACT 2  -INTER true -NON_PF true

# align on ref with bwa mem
bwa mem -M -t 4 -p GCF_000182965.3_ASM18296v3_genomic.fna ./fastq_clip/frs_242_u_clip.fastq -o sam_align/frs_242_.sam

# map aligned fastq reads back to bam
# reference sequence must also have a sequence dict created using picard
gatk --java-options "-Xmx4g" MergeBamAlignment -ALIGNED sam_align/frs_242_.sam -UNMAPPED u_bam_clip/frs_242_u_clipped.bam -R GCF_000182965.3_ASM18296v3_genomic.fna -O ./bam_align/frs_242_align.bam

# mark duplicate reads
gatk --java-options "-Xmx4g" MarkDuplicates -I ./bam_align/frs_242_align.bam -O ./bam_align/frs_242_align.marked_duplicates.bam -M frs_242_.marked_duplicates.metrics

# reorder bam
# -SD is a hidden (in the docs) variable for the reference genome
gatk --java-options "-Xmx4g" ReorderSam -I ./bam_align/frs_242_align.marked_duplicates.bam -O ./bam_align/frs_242_align.reordered.bam -SD GCF_000182965.3_ASM18296v3_genomic.fna

# build bam index
gatk --java-options "-Xmx4g" BuildBamIndex -I ./bam_align/frs_242_align.reordered.bam -O ./bam_align/frs_242_align.reordered.bai

# call haplotypes
gatk --java-options "-Xmx4g" HaplotypeCaller -R GCF_000182965.3_ASM18296v3_genomic.fna -I ./bam_align/frs_242_align.reordered.bam -O ./vcf/frs_242.g.vcf.gz -ERC GVCF -ploidy 2
# takes about 20 minutes/WGS candida albicans sample

# after all samples are processed, merge vcfs
##gatk --java-options "-Xmx4G" CombineGVCFs -R ./reference_genome/data/GCA_002759435.2_Cand_auris_B8441_V2_genomic.fna -O ./vcf/variants.merged.gvcf --variant ./vcf/SRR3883436.g.vcf.gz

# Perform joint genotyping
gatk --java-options "-Xmx4G" GenotypeGVCFs -R GCF_000182965.3_ASM18296v3_genomic.fna -O ./vcf/frs_242.out.gvcf -V ./vcf/frs_242.g.vcf.gz

#Select variants, both SNPs and raw_indels. Mixed sites are discarted.
gatk --java-options "-Xmx4G" SelectVariants -V ./vcf/frs_242.out.gvcf  -R GCF_000182965.3_ASM18296v3_genomic.fna -select-type SNP -O ./vcf/frs_242.snps.vcf.gz

gatk --java-options "-Xmx4G" SelectVariants -V ./vcf/frs_242.out.gvcf  -R GCF_000182965.3_ASM18296v3_genomic.fna -select-type INDEL -O ./vcf/frs_242.indels.vcf.gz

#Filter variants, different paramters for Indels and SNPs. Recommendations by GATK.
gatk --java-options "-Xmx4G" VariantFiltration -V ./vcf/frs_242.snps.vcf.gz -R GCF_000182965.3_ASM18296v3_genomic.fna \
-filter "QD < 12.5" --filter-name "QD12.5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-O ./vcf/frs_242.soft.snps_filtered.vcf.gz

gatk --java-options "-Xmx4G" VariantFiltration -V ./vcf/frs_242.indels.vcf.gz -R GCF_000182965.3_ASM18296v3_genomic.fna \
-filter "QD < 12.5" --filter-name "QD12.5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-O ./vcf/frs_242.soft.indels_filtered.vcf.gz

#Combine output vcfs into a single output. The output of VariantFiltration is sorted.
gatk --java-options "-Xmx4G" MergeVcfs -R GCF_000182965.3_ASM18296v3_genomic.fna \
-I ./vcf/frs_242.soft.snps_filtered.vcf.gz \
-I ./vcf/frs_242.soft.indels_filtered.vcf.gz \
-O ./vcf/frs_242_soft_filtered.vcf.gz


# make vcf human readable
bcftools view ./vcf/frs_242_soft_filtered.vcf.gz > ./vcf/frs_242_soft_filtered.nonbinary.vcf
bcftools view ./vcf/frs_242_soft_filtered.nonbinary.vcf -f "PASS" -o ./vcf/frs_242_soft_filtered.nonbinary.dropped.vcf




####################### FRS246

# convert fastq
java -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=./raw_data/246_S177_R1_001.fastq FASTQ2=./raw_data/246_S177_R2_001.fastq OUTPUT=./u_bam/frs_246_u.bam READ_GROUP_NAME=LH00162 SAMPLE_NAME=FRS246 LIBRARY_NAME=FRS246 PLATFORM_UNIT=LH00162 PLATFORM=illumina

#mark illumina adapters
gatk --java-options "-Xmx4g" MarkIlluminaAdapters -I ./u_bam/frs_246_u.bam -O ./u_bam_clip/frs_246_u_clipped.bam -M ./logs/frs246_u_markilluminaadapters_metrics.txt

# switch back to fastq with trimmed adapters
gatk --java-options "-Xmx4g" SamToFastq -I ./u_bam_clip/frs_246_u_clipped.bam --FASTQ ./fastq_clip/frs_246_u_clip.fastq  -CLIP_ATTR XT  -CLIP_ACT 2  -INTER true -NON_PF true

# align on ref with bwa mem
bwa mem -M -t 4 -p GCF_000182965.3_ASM18296v3_genomic.fna ./fastq_clip/frs_246_u_clip.fastq -o sam_align/frs_246_.sam

# map aligned fastq reads back to bam
# reference sequence must also have a sequence dict created using picard
gatk --java-options "-Xmx4g" MergeBamAlignment -ALIGNED sam_align/frs_246_.sam -UNMAPPED u_bam_clip/frs_246_u_clipped.bam -R GCF_000182965.3_ASM18296v3_genomic.fna -O ./bam_align/frs_246_align.bam

# mark duplicate reads
gatk --java-options "-Xmx4g" MarkDuplicates -I ./bam_align/frs_246_align.bam -O ./bam_align/frs_246_align.marked_duplicates.bam -M frs_246_.marked_duplicates.metrics

# reorder bam
# -SD is a hidden (in the docs) variable for the reference genome
gatk --java-options "-Xmx4g" ReorderSam -I ./bam_align/frs_246_align.marked_duplicates.bam -O ./bam_align/frs_246_align.reordered.bam -SD GCF_000182965.3_ASM18296v3_genomic.fna

# build bam index
gatk --java-options "-Xmx4g" BuildBamIndex -I ./bam_align/frs_246_align.reordered.bam -O ./bam_align/frs_246_align.reordered.bai

# call haplotypes
gatk --java-options "-Xmx4g" HaplotypeCaller -R GCF_000182965.3_ASM18296v3_genomic.fna -I ./bam_align/frs_246_align.reordered.bam -O ./vcf/frs_246.g.vcf.gz -ERC GVCF -ploidy 2
# takes about 20 minutes/WGS candida albicans sample

# after all samples are processed, merge vcfs
##gatk --java-options "-Xmx4G" CombineGVCFs -R ./reference_genome/data/GCA_002759435.2_Cand_auris_B8441_V2_genomic.fna -O ./vcf/variants.merged.gvcf --variant ./vcf/SRR3883436.g.vcf.gz

# Perform joint genotyping
gatk --java-options "-Xmx4G" GenotypeGVCFs -R GCF_000182965.3_ASM18296v3_genomic.fna -O ./vcf/frs_246.out.gvcf -V ./vcf/frs_246.g.vcf.gz

#Select variants, both SNPs and raw_indels. Mixed sites are discarted.
gatk --java-options "-Xmx4G" SelectVariants -V ./vcf/frs_246.out.gvcf  -R GCF_000182965.3_ASM18296v3_genomic.fna -select-type SNP -O ./vcf/frs_246.snps.vcf.gz

gatk --java-options "-Xmx4G" SelectVariants -V ./vcf/frs_246.out.gvcf  -R GCF_000182965.3_ASM18296v3_genomic.fna -select-type INDEL -O ./vcf/frs_246.indels.vcf.gz

#Filter variants, different paramters for Indels and SNPs. Recommendations by GATK.
gatk --java-options "-Xmx4G" VariantFiltration -V ./vcf/frs_246.snps.vcf.gz -R GCF_000182965.3_ASM18296v3_genomic.fna \
-filter "QD < 12.5" --filter-name "QD12.5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-O ./vcf/frs_246.soft.snps_filtered.vcf.gz

gatk --java-options "-Xmx4G" VariantFiltration -V ./vcf/frs_246.indels.vcf.gz -R GCF_000182965.3_ASM18296v3_genomic.fna \
-filter "QD < 12.5" --filter-name "QD12.5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-O ./vcf/frs_246.soft.indels_filtered.vcf.gz

#Combine output vcfs into a single output. The output of VariantFiltration is sorted.
gatk --java-options "-Xmx4G" MergeVcfs -R GCF_000182965.3_ASM18296v3_genomic.fna \
-I ./vcf/frs_246.soft.snps_filtered.vcf.gz \
-I ./vcf/frs_246.soft.indels_filtered.vcf.gz \
-O ./vcf/frs_246_soft_filtered.vcf.gz


# make vcf human readable
bcftools view ./vcf/frs_246_soft_filtered.vcf.gz > ./vcf/frs_246_soft_filtered.nonbinary.vcf
bcftools view ./vcf/frs_246_soft_filtered.nonbinary.vcf -f "PASS" -o ./vcf/frs_246_soft_filtered.nonbinary.dropped.vcf


