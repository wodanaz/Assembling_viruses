# find sample names based on config file inputdir and readsuffix
SAMPLES, = glob_wildcards(config["inputdir"] + "{fastqfname}" + config["readsuffix"])
if not SAMPLES:
    raise ValueError(
        "ERROR: No samples files found in {} with suffix {}.".format(config["inputdir"], config["readsuffix"]))

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

    if "datetab" in config:
        final_output.append("results/" + config["project"] + ".csv")
        final_output.append("results/" + config["project"] + "_lineages_of_concern.csv")
        final_output.append("results/supermetadata.tab")
        final_output.append("results/" + config["project"] + ".fasta")
        final_output.append("results/results.xlsx")

    if config["mode"] != 'experimental':
        final_output += expand("results/{sample}.spike.tab", sample=SAMPLES)
        final_output.append("results/spike_depths.final.tab")
        final_output += expand("results/{sample}.depth.tab", sample=SAMPLES)
        final_output.append("results/spike_genotypes.final.tab")

    return final_output
