# Assembling Viruses

### Genome Mapping Experiment 4 with BWA, trimmed TrimGalore and GATK

1. Index Genome Reference

```bash

module load bwa
bwa index -a is sars_cov_2.fasta

```

2. Remove Nextera Adapters


```bash

module load cutadapt
module load TrimGalore/0.6.5-fasrc01



ls *.fastq.gz > reads.list
for i in `cat reads.list`; do
root=`basename $i .fastq.gz`;
echo '#!/usr/bin/env bash' > $root.trimgalore.sh;
echo "trim_galore --fastqc --nextera $i  " >> $root.trimgalore.sh
done


for file in *trimgalore.sh ; do sbatch $file ; done
```


3. Map using BWA with the cleaned libraries


```bash
module load bwa

ls *_trimmed.fq.gz > reads2.list

for i in `cat reads2.list`; do
root=`basename $i _trimmed.fq.gz`;
root2=`basename $root _R1_001`;
echo '#!/usr/bin/env bash' > $root.bwa.sh;
echo "#SBATCH -N 1" >> $root.bwa.sh;
echo "bwa mem MT246667.fasta -R '@RG\tID:ID_${root2}\tPU:PU_${root2}\tSM:${root2}\tLB:${root}' $i > $root2.sam"  >> $root.bwa.sh;
done

for file in *bwa.sh ; do sbatch $file ; done

```


4. Create BAM from SAM and make an index

```bash
module load samtools

ls *.sam > sams.list
for i in `cat sams.list`; do
root=`basename $i .sam`;
echo '#!/usr/bin/env bash' > $root.sam2bam.sh;
echo "#SBATCH -N 1" >> $root.sam2bam.sh;
echo "samtools view -Sb $i | samtools sort - > $root.bam" >> $root.sam2bam.sh;
echo "samtools index  $root.bam" >> $root.sam2bam.sh;

done


for file in *sam2bam.sh ; do sbatch $file ; done
```

5. Create a dictionary file for using picard tools

```bash
module load picard-tools
java -jar /nfs/software/helmod/apps/Core/picard-tools/2.4.1-gcb01/picard.jar  CreateSequenceDictionary R=sars_cov_2.fasta O=sars_cov_2.dict 


module load samtools

samtools faidx sars_cov_2.fasta
```

# GATK starts here

1. 'Deduplicate' or mark PCR duplicates


```bash
module load picard-tools
rm bams.list
ls *.bam > bams.list
wc -l bams.list

for i in `cat bams.list`; do
root=`basename $i .bam`;
echo '#!/usr/bin/env bash' > $root.bam2dedup.sh;
echo "#SBATCH -N 1" >> $root.bam2dedup.sh;
echo "#SBATCH --mem 10G" >> $root.bam2dedup.sh;
echo "java -Xmx7g -jar /nfs/software/helmod/apps/Core/picard-tools/2.4.1-gcb01/picard.jar MarkDuplicates I=$i O=$root.dedup.bam M=$root.metric.txt" >> $root.bam2dedup.sh;
done

for file in *bam2dedup.sh ; do sbatch $file ; done

```


2. BAM TO VCF USING BCFTOOLS TO CREATE A PRELIMINAR VCF FILE


```bash

module load samtools
module load bcftools



ls *dedup.bam > dedup.list

for i in `cat dedup.list`; do
root=`basename $i .dedup.bam`;
echo '#!/usr/bin/env bash' > $root.bam2vcf.sh;
echo "#SBATCH -N 1" >> $root.bam2vcf.sh;
echo "bcftools mpileup -Ou -f sars_cov_2.fasta $i --annotate FORMAT/DPR > $root.bcf " >> $root.bam2vcf.sh;
echo "bcftools call -vm --ploidy 1 $root.bcf > $root.raw.vcf " >> $root.bam2vcf.sh;
done



for file in *bam2vcf.sh ; do sbatch $file ; done
```

3. In this step, we filter the raw SNPs using bcftools as a template of known variable sites

```bash
ls *raw.vcf > vcfs.list
wc -l vcfs.list


for i in `cat vcfs.list`; do
root=`basename $i .raw.vcf`;
echo '#!/usr/bin/env bash' > $root.vcf2filt.sh;
echo "#SBATCH -N 1" >> $root.vcf2filt.sh;
echo "bcftools view -i '%QUAL>=20 && DP>5' -Oz $i > $root.filt.vcf.gz " >> $root.vcf2filt.sh;
echo "bcftools index $root.filt.vcf.gz " >> $root.vcf2filt.sh;
done



for file in *vcf2filt.sh ; do sbatch $file ; done

```

and now, lets ad sme info for the read groups

```bash
module load picard-tools
module load GATK



ls *dedup.bam > dedup.list

for i in `cat dedup.list`; do
root=`basename $i .dedup.bam`;
echo '#!/usr/bin/env bash' > $root.bam2bam2.sh;
echo "#SBATCH -N 1" >> $root.bam2bam2.sh;
echo "java -Xmx7g -jar /nfs/software/helmod/apps/Core/picard-tools/2.4.1-gcb01/picard.jar AddOrReplaceReadGroups I=$i O=$root.bam2 RGSM=$root RGPU=unit1 RGLB=lib_${root} RGPL=ILLUMINA" >> $root.bam2bam2.sh;
done 

for file in *bam2bam2.sh ; do sbatch $file ; done


```

2. Base recalibration. First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various
user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).


```bash


module load GATK
module load tabix


ls *bam2 > bams2.list
wc -l bams2.list



for i in `cat bams2.list`; do
root=`basename $i .bam2`;
echo '#!/usr/bin/env bash' > $root.bam2br.sh;
echo "#SBATCH -N 1" >> $root.bam2br.sh;
echo "tabix -p vcf $root.filt.vcf.gz -f" >> $root.bam2br.sh;
echo "gatk --java-options -Xmx8G BaseRecalibrator -I $i -R sars_cov_2.fasta --known-sites $root.filt.vcf.gz -O $root.table " >> $root.bam2br.sh;
done



for file in *bam2br.sh ; do sbatch $file ; done

```


3. APPLY BQSR (Apply a linear base quality recalibration model trained with the BaseRecalibrator tool)

```bash
for i in `cat bams2.list`; do
root=`basename $i .bam2`;
echo '#!/usr/bin/env bash' > $root.bam2bqsr.sh;
echo "#SBATCH -N 1" >> $root.bam2bqsr.sh;
echo "gatk --java-options -Xmx8G  ApplyBQSR -I $i -R sars_cov_2.fasta --bqsr-recal-file $root.table  -O $root.bqsr.bam " >> $root.bam2bqsr.sh;
done


for file in *bam2bqsr.sh ; do sbatch $file ; done


```


4. Collect statistics: Produces a summary of alignment metrics from a SAM or BAM file


```bash

module load picard-tools
module load GATK



ls *bqsr.bam > bqsrs.list

for i in `cat bqsrs.list`; do
root=`basename $i .bqsr.bam`;
echo '#!/usr/bin/env bash' > $root.bqsr2stat.sh;
echo "#SBATCH -N 1" >> $root.bqsr2stat.sh;
echo "java -Xmx7g -jar /nfs/software/helmod/apps/Core/picard-tools/2.4.1-gcb01/picard.jar CollectAlignmentSummaryMetrics R=sars_cov_2.fasta I=$i O=$root.stat.txt" >> $root.bqsr2stat.sh;
done 

for file in *bqsr2stat.sh ; do sbatch $file ; done


```

How much of the reference genome is covered by more than 1 read?

```bash

module load samtools

for i in `cat bqsrs.list`; do
root=`basename $i .bqsr.bam`;
echo '#!/usr/bin/env bash' > $root.bqsr2depthbed.sh;
echo "#SBATCH -N 1" >> $root.bqsr2depthbed.sh;
echo "samtools depth $i -a > $root.depth.bed  " >> $root.bqsr2depthbed.sh;
done

for file in *bqsr2depthbed.sh ; do sbatch $file ; done
```

Now, let's make a table:

```bash
ls *depth.bed > depths.list
wc -l depths.list


for i in `cat depths.list`; do
root=`basename $i .depth.bed`;
percent=`awk '{ if ( $3 == 0 )  count++ } END { print 100 - ( count*100 / 29903 ) }' $i `
echo $root $percent >> table.tab
done 

sort -k1,1 -V table.tab > table.sort.tab

```

5. HAPLOTYPE CALLER: Call germline SNPs and indels via local re-assembly of haplotypes


```bash
ls *bqsr.bam > bqsrs.list

for i in `cat bqsrs.list`; do
root=`basename $i .bqsr.bam`;
echo '#!/usr/bin/env bash' > $root.bqsr2vcf.sh;
echo "#SBATCH -N 1" >> $root.bqsr2vcf.sh;
echo "gatk --java-options -Xmx8G HaplotypeCaller -I $i -R sars_cov_2.fasta -ploidy 1 -O $root.gatk.vcf" >> $root.bqsr2vcf.sh;
done


for file in *bqsr2vcf.sh ; do sbatch $file ; done

```


6. FILTER VCFs (Filter variant calls based on INFO and/or FORMAT annotations)

```
ls *gatk.vcf > vcfs.list

for i in `cat vcfs.list`; do
root=`basename $i .gatk.vcf`;
echo '#!/usr/bin/env bash' > $root.vcf2filt.sh;
echo "#SBATCH -N 1" >> $root.vcf2filt.sh;
echo "gatk VariantFiltration -R sars_cov_2.fasta -V $i -O $root.gatk.filt.vcf --filter-expression 'QD < 2.0' --filter-name QD2 --filter-expression 'FS > 60.0'  --filter-name FS60 " >> $root.vcf2filt.sh;
done


for file in *vcf2filt.sh ; do sbatch $file ; done

```


7. GENOTYPE TABLE (Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier to work with than
a VCF)


```bash

ls *filt.vcf > filts.list

for i in `cat filts.list`; do
root=`basename $i .filt.vcf`;
echo '#!/usr/bin/env bash' > $root.vcf2tab.sh;
echo "#SBATCH -N 1" >> $root.vcf2tab.sh;
echo "gatk VariantsToTable -V $i -F POS -F TYPE -F REF -F ALT -GF AD -O $root.tab " >> $root.vcf2tab.sh;
done


for file in *vcf2tab.sh ; do sbatch $file ; done

```


8. Make final consensus fasta sequence using the SNPs in the vcf file.



```bash
ls *gatk.filt.vcf > gfilts.list

module load bedtools2
ls *depth.bed > depths.list
wc -l depths.list



for i in `cat depths.list`; do
root=`basename $i .depth.bed`;
echo '#!/usr/bin/env bash' > $root.depth2maskbed.sh;
echo "#SBATCH -N 1" >> $root.depth2maskbed.sh;
echo "awk '{ if ( \$3 == 0 )  print \"NC_045512\" \"\t\" \$2 \"\t\" \$2 }' ${i} > $root.mask.bed "  >> $root.depth2maskbed.sh;
echo "bedtools merge -i $root.mask.bed | awk '{ print \$1 \"\t\" \$2  \"\t\" \$3 - 1  }' > $root.merged.bed " >> $root.depth2maskbed.sh;
#echo "bedtools maskfasta  -fi $root.fasta -bed $root.merged.bed  -fo $root.masked.fasta " >> $root.depth2maskbed.sh;
done 

for file in *depth2maskbed.sh ; do sbatch $file ; done


```

##### Note: be careful with the merged bed files. Bedtools has this annoying behavior that adds a 1 to the end part of the BED file if end is not 1. Just run: 



```bash

head *merged.bed

```
and check there aren't any bed files that have a 0 in end position in the first line:

for example:

NC_045512 1 0
NC_045512 200 290

If so, just edit 0 for 1 using nano



```bash
ls *gatk.filt.vcf > gfilts.list


for i in `cat gfilts.list`; do
root=`basename $i .gatk.filt.vcf`;
echo '#!/usr/bin/env bash' > $root.vcf2fasta.sh;
echo "#SBATCH -N 1" >> $root.vcf2fasta.sh;
echo "bcftools view -Oz $i > $i.gz " >> $root.vcf2fasta.sh;
echo "bcftools index $i.gz " >> $root.vcf2fasta.sh;
echo "bcftools norm -f sars_cov_2.fasta $i.gz -Ob -o $root.norm.bcf " >> $root.vcf2fasta.sh;
echo "bcftools consensus -m $root.merged.bed -f sars_cov_2.fasta  -p ${root}_ -s ${root} -H A $i.gz > $root.masked.fasta " >> $root.vcf2fasta.sh;
done

for file in *vcf2fasta.sh ; do sbatch $file ; done


```


9. Compile all tab tables into one for depth and genotype


make sure you doanloaded genotype_compiler.pl and depth_compiler.pl

```bash
perl depth_compiler.pl *.depth.tab > alldepths.tab
cp alldepths.tab alldepths.backup.tab
sed -r 's/ /\t/g'  alldepths.tab | sed -r 's/.depth.tab//g' | sort -k1,1n > alldepths.final.tab 


perl genotype_compiler.pl *.filt.tab > allgenotypes.tab
cp allgenotypes.tab allgenotypes.backup.tab
sed -r 's/ /\t/g'  allgenotypes.tab | sed -r 's/.filt.tab//g' | sort -k1,1n > allgenotypes.final.tab 

```









