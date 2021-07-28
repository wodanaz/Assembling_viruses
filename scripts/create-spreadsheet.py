#!/usr/bin/env python3
# Create a spreadsheet combining supermetadata and genotype data
#
# Required Environment Variables:
#  EVDIR - working directory containing supermetadata.tab and spike_genotypes.final.tab
#  PROJECTNAME - name of the project - used to determine the pangolin output csv
#
import pandas as pd
import os
from dateutil.parser import parse

# Setup input files and output path based on $EVDIR
evdir = os.environ.get("EVDIR") + "/"
supermetadata_path = evdir + "supermetadata.tab"
spike_genotypes_path = evdir + "spike_genotypes.final.tab"
output_path = evdir + "results.xlsx"
spike_depths_path = evdir + "spike_depths.final.tab"
pangolin_csv = evdir + os.environ.get("PROJECTNAME") + ".csv"

# Specifies when "week" 1 begins "every new week starts on a Friday"
week_start_date = parse("Jan 4 2021")

#############
# Functions #
#############

def create_metadata_month(onset_date):
    return parse(onset_date).strftime("%b_%Y")

def create_metadata_week(onset_date):
    # determine the number of days elapsed since the start date
    delta = parse(onset_date) - week_start_date
    # divide by 7 to get weeks since start date
    week_num = int(delta.days / 7)
    # add 1 so the first week is "week_1"
    week_num +=1
    # only return a week str if we have a valid input value
    if week_num >= 1:
        return "week_{}".format(week_num)
    return None

#####################################
# Create Spike Genotypes data frame #
#####################################
# read supermetadata using Sample as an index column
supermetadata = pd.read_csv(supermetadata_path, sep="$", index_col=0,
                            names=["Sample", "Lineage", "Coverage", "Date"])

# read supermetadata using Sample as an index column
spike_genotypes = pd.read_csv(spike_genotypes_path, sep="\t", index_col=0)

# drop na columns and transpose
transposed_spike_genotypes = spike_genotypes.dropna(axis=1).transpose()

# join
spike_genotypes_df = supermetadata.join(transposed_spike_genotypes)

#################################
# Create Read Depths data frame #
#################################
# read spike_depths.final.tab into a data frame
raw_read_depths_df = pd.read_csv(spike_depths_path, sep="\t", index_col=0)
# drop na columns and transpose
read_depths_df = raw_read_depths_df.dropna(axis=1).transpose()

#######################################
# Create Pangolin Lineages data frame #
#######################################
# read <projectname>.csv (pangolin output file) into a data frame
pangolin_lineages_df = pd.read_csv(pangolin_csv, index_col="taxon")

##############################
# Create Metadata data frame #
##############################
metadata_df = supermetadata
# Add Week and Month columns to supermetdata dataframe
metadata_df['Week'] = supermetadata['Date'].apply(create_metadata_week)
metadata_df['Month'] = supermetadata['Date'].apply(create_metadata_month)

# save as spreadsheet
with pd.ExcelWriter(output_path, mode='w') as writer:
    spike_genotypes_df.to_excel(writer, sheet_name='Spike Genotypes')
    read_depths_df.to_excel(writer, sheet_name='Read Depths')
    pangolin_lineages_df.to_excel(writer, sheet_name='Pangolin Lineages')
    metadata_df.to_excel(writer, sheet_name='Metadata')
