#!/bin/bash
# GENOTYPE TABLE (Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier
# to work with than a VCF)
# Required positional arguments:
# - vcf file to process

# stop if a command fails (non-zero exit status)
set -e

VCFFILE=$1

root=`basename $VCFFILE .filt.vcf`
gatk --java-options "-Djava.io.tmpdir=." VariantsToTable -V $VCFFILE -F POS -F TYPE -F REF -F ALT -GF AD -O $root.tab
