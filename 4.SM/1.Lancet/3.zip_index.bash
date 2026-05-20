#!/bin/bash

source /Path/to/env/bin/activate 

GzipIndex () {
bgzip -@ 10 $1
tabix -p vcf $1".gz"
}

for vcf in fixed_header.vcf P19.onlyPASS_header_fixed.vcf; do
GzipIndex $vcf
done
