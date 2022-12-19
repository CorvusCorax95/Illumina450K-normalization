import math

import pandas as pd

import numpy as np


### FUNCTIONS ###

# making a help-dataframe with just probe-IDs and
# probe types (I and II)
def probetype_to_dataframe(input):
    # make dataframe from csv-file
    df = pd.read_csv(input, sep=',')

    # set probes as index
    df_i = df.set_index('probe')
    return df_i


# adding probetypes (I and II) to my dataframe
# to distinguish between Infinium I and Infinium II probes
def add_probetypes(df_dest, df_target):
    # search for fitting probe-ids
    df_merged = pd.merge(df_dest, df_target, on="probe")

    return df_merged


def dataframe_log(x):
    if (isinstance(x, int) and x != 0):
        math.log(x, math.e)
    return x


### MAIN ###
# rebuild illumina probe tables to dataframe with just probes and types
def prepare_data():

    # returns dataframe
    df_450 = probetype_to_dataframe('illumina-table-450.csv')

    # ## ORIGINAL FILES ##
    # # convert methylation-tables to dataframes
    # df_values_meth = pd.read_csv('methylation-rb/methylation_meth.wide.tsv', sep='\t')
    # df_values_unmeth = pd.read_csv('methylation-rb/methylation_unmeth.wide.tsv', sep='\t')

    ## SHORT FILES ##
    # convert methylation-tables to dataframes
    df_values_meth = pd.read_csv('methylation-rb/short_meth.tsv', sep='\t')
    df_values_unmeth = pd.read_csv('methylation-rb/short_unmeth.tsv', sep='\t')

    # set index to probe-id
    df_values_meth = df_values_meth.set_index('probe')
    df_values_unmeth = df_values_unmeth.set_index('probe')

    # match probes to type in methylation file
    # returns new dataframe including probe types
    df_meth_w_types = add_probetypes(df_450, df_values_meth)
    df_unmeth_w_types = add_probetypes(df_450, df_values_unmeth)

    df_meth_w_types.to_csv('short_methylated_w_types.csv', sep='\t')
    df_unmeth_w_types.to_csv('short_unmethylated_w_types.csv', sep='\t')


# reads pandas dataframe from csv
def make_dataframe(file):
    df = pd.read_csv(file, sep='\t')
    return df


# applying log to all values in the dataframe
# stores dataframe as csv in new file
def log_data(df, filename):
    # log values of dataframe
    df.set_index("probe", inplace=True)

    sample_list = df.columns.values.tolist()[2:]

    df_log = np.log2(df[sample_list])
    df_log.to_csv(filename)

    # find and remove inf values after log
    df_log.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_log = df_log.dropna()
    return df_log


# calculates means of all samples
def output_measures(df, title):
    sample_list = df.columns.values.tolist()[1:]
    print("Dataframe: " + title)
    means = []
    for x in sample_list:
        means.append(df[x].mean())
        # print(x, ": ", str(df[x].mean())
    print("Mean: ", np.average(means))


# beta-value: methylated / methylated + unmethylated + 100
def beta_value(df_meth, df_unmeth, offset):
    df = df_meth
    sample_list = df_meth.columns.values.tolist()[1:]
    probes_meth = df_meth.index.values.tolist()
    # loops through all columns
    for x in sample_list:
        for y in probes_meth:
            df[x][y] = df_meth[x][y] / (df_meth[x][y] + df_unmeth[x][y] + offset)
    return df


# M-value: log(methylated / unmethylated)
def m_value(df_meth, df_unmeth):
    df = df_meth
    sample_list = df_meth.columns.values.tolist()[1:]
    probes_meth = df_meth.index.values.tolist()
    # loops through all columns
    for x in sample_list:
        # loops through all rows
        for y in probes_meth:
            inner = df_meth[x][y] / df_unmeth[x][y]
            df[x][y] = np.log2(inner)
    return df
