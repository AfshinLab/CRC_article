#!/bin/bash


#merge
less chr1_out.vcf | head -1000 | grep "^#" | sed '40s/20\t20/normal\ttumor/g' > fixed_header.vcf &&
grep -vh '^#' chr{1..22}_out.vcf >> fixed_header.vcf
grep -vh '^#' chrX_out.vcf >> fixed_header.vcf
grep -vh '^#' chrY_out.vcf >> fixed_header.vcf

# filter
less fixed_header.vcf | grep "^#" > P19.onlyPASS_header_fixed.vcf &&
less fixed_header.vcf | grep PASS >> P19.onlyPASS_header_fixed.vcf

