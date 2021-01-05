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

1. 'Deduplicate or mark PCR duplicates'


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

