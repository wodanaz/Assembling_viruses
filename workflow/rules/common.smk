import os
import glob

# create dictionary of sample name to filename
SAMPLE_TO_FILENAME = dict()
# find FASTQ files based on config file inputdir
for sample_path in glob.glob(config["inputdir"] + "/*.fastq.gz"):
    # sample name is the beginning of the filename before the first "_"
    sample_name = os.path.basename(sample_path).split("_")[0]
    SAMPLE_TO_FILENAME[sample_name] = sample_path
SAMPLES = SAMPLE_TO_FILENAME.keys()
if not SAMPLES:
    raise ValueError( "ERROR: No samples files found in {}.".format(config["inputdir"]))

# determine genome index filenames
GENOME_INDEX_FILES = multiext(config["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
GENOME_INDEX_FILES.append(config["genome"].replace(".fasta", ".dict"))

# Defines the files to be created by the workflow
def get_final_output():
    final_output = []
    final_output += expand("results/{sample}.gatk.tab", sample=SAMPLES)
    final_output += expand("results/{sample}.gatk.filt.vcf", sample=SAMPLES)
    final_output.append("results/coverage.gatk.tab")
    final_output.append("results/coverage.raw.tab")
    final_output += expand("results/{sample}.cleaned.fasta", sample=SAMPLES)
    final_output += expand("results/{sample}.spike.tab", sample=SAMPLES)
    final_output.append("results/spike_depths.final.tab")
    final_output += expand("results/{sample}.depth.tab", sample=SAMPLES)
    final_output.append("results/spike_genotypes.final.tab")
    final_output.append("results/" + config["project"] + ".csv")
    return final_output
