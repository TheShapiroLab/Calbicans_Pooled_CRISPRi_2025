#!/bin/bash

module load bedtools

bedtools genomecov -ibam ./bam_align/frs_1_align.reordered.bam > ./coverage_hists/frs_1_align.cov.txt
bedtools genomecov -ibam ./bam_align/frs_242_align.reordered.bam > ./coverage_hists/frs_242_align.cov.txt
bedtools genomecov -ibam ./bam_align/frs_246_align.reordered.bam > ./coverage_hists/frs_246_align.cov.txt
