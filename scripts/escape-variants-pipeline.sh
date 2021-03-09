#!/usr/bin/env bash
# Investigating escape mutations in SARS-CoV-2
# - Processes all *.fastq.gz files in the current directory.
#
#SBATCH --job-name=ev-pipeline

# stop if a command fails (non-zero exit status)
set -e

# show an error message when a step fails
trap "echo 'ERROR: Script failed'" ERR

# record the hash of the latest commit of this script
export EV_GIT_COMMIT=$(git log --pretty="%H" -n 1)

echo "Pipeline - Starting - Version $EV_GIT_COMMIT"
echo ""

echo "Current directory is $(pwd)"
echo ""

echo "Creating a temp directory"
# create temp directory starting within the $WORKDIR
export EVDIR=$(mktemp -d --tmpdir=$WORKDIR)
echo "Running using directory $EVDIR"
echo ""


echo "Exported Environment Variables"
echo "    export EVDIR=$EVDIR"
echo "    export GENOME=$GENOME"
echo "    export INPUTDIR=$INPUTDIR"
echo ""


echo "Step 1 - Remove Nextera Adapters"
# create the list of input fastq.gz filenames to process
ls ${INPUTDIR}/*.fastq.gz > $EVDIR/reads.list
./scripts/sbatch-array.sh \
    ./scripts/remove-nextera-adapters.sh $EVDIR/reads.list
echo "Step 1 - Done"
echo ""


echo "Step 2 - Map using BWA with the cleaned libraries"
# create the list of trimmed reads to process
ls $EVDIR/*_trimmed.fq.gz > $EVDIR/reads2.list
./scripts/sbatch-array.sh \
    ./scripts/map-bwa-cleaned-libs.sh $EVDIR/reads2.list
echo "Step 2 - Done"
echo ""


echo "Step 3 - Create BAM from SAM and make an index"
# create the list of sam files to process
ls $EVDIR/*.sam > $EVDIR/sams.list
./scripts/sbatch-array.sh \
    ./scripts/create-bam-from-sam.sh $EVDIR/sams.list
echo "Step 3 - Done"
echo ""


echo "GATK Step 1 - 'Deduplicate' or mark PCR duplicates"
ls $EVDIR/*.bam > $EVDIR/bams.list
./scripts/sbatch-array.sh \
    ./scripts/mark-duplicates.sh $EVDIR/bams.list
echo "GATK Step 1 - Done"
echo ""


echo "GATK Step 2 - BAM TO VCF USING BCFTOOLS TO CREATE A PRELIMINAR VCF FILE"
ls $EVDIR/*dedup.bam > $EVDIR/dedup.list
./scripts/sbatch-array.sh \
    ./scripts/convert-bam-to-vcf.sh $EVDIR/dedup.list
echo "GATK Step 2 - Done"
echo ""


echo "GATK Step 3a - In this step, we filter the raw SNPs using bcftools as a template of known variable sites"
ls $EVDIR/*raw.vcf > $EVDIR/vcfs.list
./scripts/sbatch-array.sh \
    ./scripts/filter-raw-snps.sh $EVDIR/vcfs.list
echo "GATK Step 3a - Done"
echo ""


echo "GATK Step 3b - Add sme info for the read groups"
ls $EVDIR/*dedup.bam > $EVDIR/dedup.list
./scripts/sbatch-array.sh \
    ./scripts/update-read-groups.sh $EVDIR/dedup.list
echo "GATK Step 3b - Done"
echo ""


echo "GATK Step 4 - Base recalibration. First pass of the Base Quality Score Recalibration (BQSR)"
ls $EVDIR/*bam2 > $EVDIR/bams2.list
./scripts/sbatch-array.sh \
    ./scripts/run-base-recalibration.sh $EVDIR/bams2.list
echo "GATK Step 4 - Done"
echo ""


echo "GATK Step 5 - APPLY BQSR (Apply a linear base quality recalibration model trained with the BaseRecalibrator tool)"
./scripts/sbatch-array.sh \
    ./scripts/apply-bqsr.sh $EVDIR/bams2.list
echo "GATK Step 5 - Done"
echo ""


echo "GATK Step 6a - Collect statistics: Produces a summary of alignment metrics from a SAM or BAM file"
ls $EVDIR/*bqsr.bam > $EVDIR/bqsrs.list
./scripts/sbatch-array.sh \
    ./scripts/collect-statistics.sh $EVDIR/bqsrs.list
echo "GATK Step 6a - Done"
echo ""


echo "GATK Step 6b How much of the reference genome is covered by more than 1 read?"
./scripts/sbatch-array.sh \
    ./scripts/run-depth-coverage.sh $EVDIR/bqsrs.list
echo "GATK Step 6b - Done"
echo ""


echo "GATK Step 6c - make a table with coverage information"
echo "Running:"
echo "    ./scripts/create-coverage-table.sh"
sbatch --wait "--output=${LOGDIR}/create-coverage-table-%j.out" \
    ./scripts/create-coverage-table.sh
echo "GATK Step 6c - Done"
echo ""


echo "GATK Step 7 - HAPLOTYPE CALLER: Call germline SNPs and indels via local re-assembly of haplotypes"
ls $EVDIR/*bqsr.bam > $EVDIR/bqsrs.list
./scripts/sbatch-array.sh \
    ./scripts/run-haplotype-caller.sh $EVDIR/bqsrs.list
echo "GATK Step 7 - Done"
echo ""


echo "GATK Step 8 - FILTER VCFs (Filter variant calls based on INFO and/or FORMAT annotations)"
ls $EVDIR/*gatk.vcf > $EVDIR/vcfs.list
./scripts/sbatch-array.sh \
    ./scripts/filter-vcfs.sh $EVDIR/vcfs.list
echo "GATK Step 8 - Done"
echo ""


echo "GATK Step 9 - GENOTYPE TABLE (Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier to work with than a VCF)"
ls $EVDIR/*filt.vcf > $EVDIR/filts.list
./scripts/sbatch-array.sh \
    ./scripts/create-genotype-table.sh $EVDIR/filts.list
echo "GATK Step 9 - Done"
echo ""


echo "GATK Step 10a - Make final consensus fasta sequence using the SNPs in the vcf file"
ls $EVDIR/*depth.bed > $EVDIR/depths.list
./scripts/sbatch-array.sh \
    ./scripts/run-bedtools-merge.sh $EVDIR/depths.list
echo "GATK Step 10a - Done"
echo ""


echo "GATK Step 10b - Run bcftools consensus"
ls $EVDIR/*gatk.filt.vcf > $EVDIR/gfilts.list
./scripts/sbatch-array.sh \
    ./scripts/run-bcftools-consensus.sh $EVDIR/gfilts.list
echo "GATK Step 10b - Done"
echo ""


if [ "$SURVEILLANCE_MODE" == "Y" ]
then
    echo "Skipping GATK Step 11 since running in surveillance mode."
else
    echo "GATK Step 11a - bcftools query for '%POS %ALT' and '%POS [%AD]'"
    ls $EVDIR/*.gatk.filt.vcf > $EVDIR/vcfs2.list
    echo "Running:"
    echo "    ./scripts/run-bcftools-query-alt-ad.sh"
    sbatch --wait "--output=${LOGDIR}/run-bcftools-query-alt-ad-%j.out" \
        ./scripts/run-bcftools-query-alt-ad.sh
    echo "GATK Step 11a - Done"
    echo ""


    echo "GATK Step 11b - Compile all tab tables for depth"
    echo "Running:"
    echo "    ./scripts/run-depth-compiler.sh"
    sbatch --wait "--output=${LOGDIR}/run-depth-compiler-%j.out" \
        ./scripts/run-depth-compiler.sh
    echo "GATK Step 11b - Done"
    echo ""


    echo "GATK Step 11c - Compile all tab tables for genotype"
    echo "Running:"
    echo "    ./scripts/run-genotype-compiler.sh"
    sbatch --wait "--output=${LOGDIR}/run-genotype-compiler-%j.out" \
        ./scripts/run-genotype-compiler.sh
    echo "GATK Step 11c - Done"
    echo ""
fi


echo "Copying output files to $OUTDIR"
mkdir -p $OUTDIR
# saving standard output files
cp $EVDIR/*gatk.tab $OUTDIR/.
cp $EVDIR/*gatk.filt.vcf.gz $OUTDIR/.
cp $EVDIR/coverage.gatk.tab $OUTDIR/.
cp $EVDIR/coverage.raw.tab $OUTDIR/.
cp $EVDIR/*cleaned.fasta $OUTDIR/.
# when not in surveillance save GATK Step 11 files
if [ "$SURVEILLANCE_MODE" != "Y" ]
then
   cp $EVDIR/*.filt.tab $OUTDIR/.
   cp $EVDIR/*.depth.tab $OUTDIR/.
   cp $EVDIR/alldepths.final.tab $OUTDIR/.
   cp $EVDIR/allgenotypes.final.tab $OUTDIR/.
fi


if [ "$DELETE_EVTMPDIR" == "Y" ]
    then
    echo "Deleting temporary $EVDIR directory"
    # go into the parent directory
    cd ..
    # then we can delete the EVTMPDIR
    rm -rf $EVDIR
fi

echo ""
echo "Pipeline - Done - Results in $OUTDIR"
