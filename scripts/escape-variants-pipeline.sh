#!/usr/bin/env bash
# Investigating escape mutations in SARS-CoV-2
# - Processes all *.fastq.gz files in the current directory.
#
#SBATCH --job-name=ev-pipeline

# try to record the hash of the latest commit of this script
GIT_COMMIT=$(git log --pretty="%H" -n 1 $EVBASEDIR)

# stop if a command fails (non-zero exit status)
set -e

# show an error message when a step fails
trap "echo 'ERROR: Script failed'" ERR

echo "Pipeline - Starting - Version $GIT_COMMIT"

echo "Creating a temp directory"
# create temp directory starting within the $WORKDIR
export EVTMPDIR=$(mktemp -d --tmpdir=$WORKDIR)
# change into the temp directory
cd $EVTMPDIR
echo "Running in temp directory $EVTMPDIR"

echo "Step 1 - Remove Nextera Adapters"
# create the list of input fastq.gz filenames to process
ls ${INPUTDIR}/*.fastq.gz > reads.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/remove-nextera-adapters.sh reads.list
echo "Step 1 - Done"
echo ""


echo "Step 2 - Map using BWA with the cleaned libraries"
# create the list of trimmed reads to process
ls *_trimmed.fq.gz > reads2.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/map-bwa-cleaned-libs.sh reads2.list
echo "Step 2 - Done"
echo ""


echo "Step 3 - Create BAM from SAM and make an index"
# create the list of sam files to process
ls *.sam > sams.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/create-bam-from-sam.sh sams.list
echo "Step 3 - Done"
echo ""


echo "GATK Step 1 - 'Deduplicate' or mark PCR duplicates"
ls *.bam > bams.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/mark-duplicates.sh bams.list
echo "GATK Step 1 - Done"
echo ""


echo "GATK Step 2 - BAM TO VCF USING BCFTOOLS TO CREATE A PRELIMINAR VCF FILE"
ls *dedup.bam > dedup.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/convert-bam-to-vcf.sh dedup.list
echo "GATK Step 2 - Done"
echo ""


echo "GATK Step 3a - In this step, we filter the raw SNPs using bcftools as a template of known variable sites"
ls *raw.vcf > vcfs.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/filter-raw-snps.sh vcfs.list
echo "GATK Step 3a - Done"
echo ""


echo "GATK Step 3b - Add sme info for the read groups"
ls *dedup.bam > dedup.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/update-read-groups.sh  dedup.list
echo "GATK Step 3b - Done"
echo ""


echo "GATK Step 4 - Base recalibration. First pass of the Base Quality Score Recalibration (BQSR)"
ls *bam2 > bams2.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/run-base-recalibration.sh bams2.list
echo "GATK Step 4 - Done"
echo ""


echo "GATK Step 5 - APPLY BQSR (Apply a linear base quality recalibration model trained with the BaseRecalibrator tool)"
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/apply-bqsr.sh bams2.list
echo "GATK Step 5 - Done"
echo ""


echo "GATK Step 6a - Collect statistics: Produces a summary of alignment metrics from a SAM or BAM file"
ls *bqsr.bam > bqsrs.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/collect-statistics.sh bqsrs.list
echo "GATK Step 6a - Done"
echo ""


echo "GATK Step 6b How much of the reference genome is covered by more than 1 read?"
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/run-depth-coverage.sh bqsrs.list
echo "GATK Step 6b - Done"
echo ""


echo "GATK Step 6c - make a table with coverage information"
ls *depth.bed > depths.list
sbatch --wait "--output=${LOGDIR}/create-coverage-table-%j.out" $EVSCRIPTS/create-coverage-table.sh depths.list
echo "GATK Step 6c - Done"
echo ""


echo "GATK Step 7 - HAPLOTYPE CALLER: Call germline SNPs and indels via local re-assembly of haplotypes"
ls *bqsr.bam > bqsrs.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/run-haplotype-caller.sh bqsrs.list
echo "GATK Step 7 - Done"
echo ""


echo "GATK Step 8 - FILTER VCFs (Filter variant calls based on INFO and/or FORMAT annotations)"
ls *gatk.vcf > vcfs.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/filter-vcfs.sh vcfs.list
echo "GATK Step 8 - Done"
echo ""


echo "GATK Step 9 - GENOTYPE TABLE (Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier to work with than a VCF)"
ls *filt.vcf > filts.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/create-genotype-table.sh filts.list
echo "GATK Step 9 - Done"
echo ""


echo "GATK Step 10a - Make final consensus fasta sequence using the SNPs in the vcf file"
ls *depth.bed > depths.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/run-bedtools-merge.sh depths.list
echo "GATK Step 10a - Done"
echo ""


#echo "GATK Step 10b - Fix merged bed positions"
#sbatch --wait "--output=${LOGDIR}/fix-merged-bed-positions-%j.out" $EVSCRIPTS/fix-merged-bed-positions.sh
#echo "GATK Step 10b - Done"
#echo ""


echo "GATK Step 10b - Run bcftools consensus"
ls *gatk.filt.vcf > gfilts.list
$EVSCRIPTS/sbatch-array.sh $EVSCRIPTS/run-bcftools-consensus.sh gfilts.list
echo "GATK Step 10b - Done"
echo ""


if [ "$SURVEILLANCE_MODE" == "Y" ]
then
    echo "Skipping GATK Step 11 since running in surveilance mode."
else
    echo "GATK Step 11a - bcftools query for '%POS %ALT' and '%POS [%AD]'"
    ls *.gatk.filt.vcf > vcfs2.list
    sbatch --wait "--output=${LOGDIR}/run-bcftools-query-alt-ad-%j.out" $EVSCRIPTS/run-bcftools-query-alt-ad.sh
    echo "GATK Step 11a - Done"
    echo ""


    echo "GATK Step 11b - Compile all tab tables for depth"
    sbatch --wait "--output=${LOGDIR}/run-depth-compiler-%j.out" $EVSCRIPTS/run-depth-compiler.sh
    echo "GATK Step 11b - Done"
    echo ""


    echo "GATK Step 11c - Compile all tab tables for genotype"
    sbatch --wait "--output=${LOGDIR}/run-genotype-compiler-%j.out" $EVSCRIPTS/run-genotype-compiler.sh
    echo "GATK Step 11c - Done"
    echo ""
fi


echo "Copying output files to $OUTDIR"
mkdir -p $OUTDIR
# saving *masked.fasta, the *gatk.tab file, the *gatk.filt.vcf.gz, and table.sort.tab output files
cp *gatk.tab $OUTDIR/.
cp *gatk.filt.vcf.gz $OUTDIR/.
cp coverage.gatk.tab $OUTDIR/.
cp coverage.raw.tab $OUTDIR/.
cp *cleaned.fasta $OUTDIR/.


if [ "$DELETE_EVTMPDIR" == "Y" ]
    then
    echo "Deleting temporary $EVTMPDIR directory"
    # go into the parent directory
    cd ..
    # then we can delete the EVTMPDIR
    rm -rf $EVTMPDIR
fi


echo "Pipeline - Done"
