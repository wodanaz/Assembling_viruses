#!/usr/bin/env bash
# Investigating escape mutations in SARS-CoV-2
# - Indexes MT246667.fasta
# - Processes all *.fastq.gz files in the current directory.
#
#SBATCH --job-name=ev-pipeline
#SBATCH --output=logs/ev-pipeline-%j.out

# stop if a command fails (non-zero exit status)
set -e

echo "Pipeline - Starting"


echo "Step 1 - Index Genome Reference"
sbatch --wait scripts/index-reference-genome.sh
echo "Step 1 - Done"
echo ""


echo "Step 2 - Remove Nextera Adapters"
# create the list of fastq.gz filenames to process
ls *.fastq.gz > reads.list
./scripts/sbatch-array.sh scripts/remove-nextera-adapters.sh reads.list
echo "Step 2 - Done"
echo ""


echo "Step 3 - Map using BWA with the cleaned libraries"
# create the list of trimmed reads to process
ls *_trimmed.fq.gz > reads2.list
./scripts/sbatch-array.sh scripts/map-bwa-cleaned-libs.sh reads2.list
echo "Step 3 - Done"
echo ""


echo "Step 4 - Create BAM from SAM and make an index"
# create the list of sam files to process
ls *.sam > sams.list
./scripts/sbatch-array.sh scripts/create-bam-from-sam.sh sams.list
echo "Step 4 - Done"
echo ""


echo "Step 5 - Create a dictionary file for using picard tools"
sbatch --wait scripts/create-picard-dictionary.sh
echo "Step 5 - Done"
echo ""


echo "GATK Step 1 - 'Deduplicate' or mark PCR duplicates"
ls *.bam > bams.list
./scripts/sbatch-array.sh scripts/mark-duplicates.sh bams.list
echo "GATK Step 1 - Done"
echo ""


echo "GATK Step 2 - BAM TO VCF USING BCFTOOLS TO CREATE A PRELIMINAR VCF FILE"
ls *dedup.bam > dedup.list
./scripts/sbatch-array.sh scripts/convert-bam-to-vcf.sh dedup.list
echo "GATK Step 2 - Done"
echo ""


echo "GATK Step 3a - In this step, we filter the raw SNPs using bcftools as a template of known variable sites"
ls *raw.vcf > vcfs.list
./scripts/sbatch-array.sh scripts/filter-raw-snps.sh vcfs.list
echo "GATK Step 3a - Done"
echo ""


echo "GATK Step 3b - Add sme info for the read groups"
ls *dedup.bam > dedup.list
./scripts/sbatch-array.sh scripts/update-read-groups.sh  dedup.list
echo "GATK Step 3b - Done"
echo ""


echo "GATK Step 4 - Base recalibration. First pass of the Base Quality Score Recalibration (BQSR)"
ls *bam2 > bams2.list
./scripts/sbatch-array.sh scripts/run-base-recalibration.sh bams2.list
echo "GATK Step 4 - Done"
echo ""


echo "GATK Step 5 - APPLY BQSR (Apply a linear base quality recalibration model trained with the BaseRecalibrator tool)"
./scripts/sbatch-array.sh scripts/apply-bqsr.sh bams2.list
echo "GATK Step 5 - Done"
echo ""


echo "GATK Step 6a - Collect statistics: Produces a summary of alignment metrics from a SAM or BAM file"
ls *bqsr.bam > bqsrs.list
./scripts/sbatch-array.sh scripts/collect-statistics.sh bqsrs.list
echo "GATK Step 6a - Done"
echo ""


echo "GATK Step 6b How much of the reference genome is covered by more than 1 read?"
./scripts/sbatch-array.sh scripts/run-depth-coverage.sh bqsrs.list
echo "GATK Step 6b - Done"
echo ""


echo "GATK Step 6c - make a table with coverage information"
ls *depth.bed > depths.list
sbatch --wait  ./scripts/create-coverage-table.sh depths.list
echo "GATK Step 6c - Done"
echo ""


echo "GATK Step 7 - HAPLOTYPE CALLER: Call germline SNPs and indels via local re-assembly of haplotypes"
ls *bqsr.bam > bqsrs.list
./scripts/sbatch-array.sh scripts/run-haplotype-caller.sh bqsrs.list
echo "GATK Step 7 - Done"
echo ""


echo "GATK Step 8 - FILTER VCFs (Filter variant calls based on INFO and/or FORMAT annotations)"
ls *gatk.vcf > vcfs.list
./scripts/sbatch-array.sh scripts/filter-vcfs.sh vcfs.list
echo "GATK Step 8 - Done"
echo ""


echo "GATK Step 9 - GENOTYPE TABLE (Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier to work with than a VCF)"
ls *filt.vcf > filts.list
./scripts/sbatch-array.sh scripts/create-genotype-table.sh filts.list
echo "GATK Step 9 - Done"
echo ""


echo "GATK Step 10a - Make final consensus fasta sequence using the SNPs in the vcf file"
ls *depth.bed > depths.list
./scripts/sbatch-array.sh scripts/run-bedtools-merge.sh depths.list
echo "GATK Step 10a - Done"
echo ""


echo "GATK Step 10b - Fix merged bed positions"
sbatch --wait  ./scripts/fix-merged-bed-positions.sh
echo "GATK Step 10b - Done"
echo ""


echo "GATK Step 10c - Run bcftools consensus"
ls *gatk.filt.vcf > gfilts.list
./scripts/sbatch-array.sh scripts/run-bcftools-consensus.sh gfilts.list
echo "GATK Step 10c - Done"
echo ""


echo "GATK Step 11a - bcftools query for '%POS %ALT' and '%POS [%AD]'"
ls *.gatk.filt.vcf > vcfs2.list
sbatch --wait  ./scripts/run-bcftools-query-alt-ad.sh
echo "GATK Step 11a - Done"
echo ""


echo "GATK Step 11b - Compile all tab tables for depth"
sbatch --wait  ./scripts/run-depth-compiler.sh
echo "GATK Step 11b - Done"
echo ""


echo "GATK Step 11c - Compile all tab tables for genotype"
sbatch --wait  ./scripts/run-genotype-compiler.sh
echo "GATK Step 11c - Done"
echo ""


echo "Pipeline - Done"
