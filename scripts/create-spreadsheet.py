#!/usr/bin/env python3
# Create a spreadsheet combining supermetadata and genotype data
#
# Required Environment Variables:
#  EVDIR - working directory containing supermetadata.tab and spike_genotypes.final.tab
#
import pandas as pd
import os

# Setup input files and output path based on $EVDIR
evdir = os.environ.get("EVDIR") + "/"
supermetadata_path = evdir + "supermetadata.tab"
spike_genotypes_path = evdir + "spike_genotypes.final.tab"
output_path = evdir + "results.xlsx"

# read supermetadata using Sample as an index column
supermetadata = pd.read_csv(supermetadata_path, sep="$", index_col=0,
                            names=["Sample", "Lineage", "Coverage", "Date"])

# read supermetadata using Sample as an index column
spike_genotypes = pd.read_csv(spike_genotypes_path, sep="\t", index_col=0)

# drop na columns and transpose
transposed_spike_genotypes = spike_genotypes.dropna(axis=1).transpose()

# join and save as spreadsheet
joined_df = supermetadata.join(transposed_spike_genotypes)
joined_df.to_excel(output_path, sheet_name='Spike Genotypes')
