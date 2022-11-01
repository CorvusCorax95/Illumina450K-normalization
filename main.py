import math

import numpy as np
import pandas as pd
import qnorm as qn
import matplotlib.pyplot as plt
import scipy

### FUNCTIONS ###
def probetype_to_dataframe(input):
    # make dataframe from csv-file
    df = pd.read_csv(input, sep=',')

    # set probes as index
    df_i = df.set_index('probe')
    return df_i

def add_probetypes(df_dest, df_target):

    # search for fitting probe-ids
    df_merged = pd.merge(df_dest, df_target, on="probe")

    return df_merged

def dataframe_log(x):
    if(isinstance(x, int) and x != 0):
        math.log(x, math.e)
    return x


### MAIN ###

# rebuild illumina probe tables to dataframe with just probes and types
# returns dataframe
df_450 = probetype_to_dataframe('illumina-table-450.csv')


# convert methylation-tables to dataframes
df_values_meth = pd.read_csv('methylation-rb/short.tsv', sep='\t')

#df_values_meth = pd.read_csv('methylation-rb/methylation_meth.wide.tsv', sep='\t')
df_values_unmeth = pd.read_csv('methylation-rb/methylation_unmeth.wide.tsv', sep='\t')

# set index to probe-id
df_values_meth = df_values_meth.set_index('probe')
df_values_unmeth = df_values_unmeth.set_index('probe')

# match probes to type in methylation file
# returns new dataframe including probe types
df_meth_w_types = add_probetypes(df_450, df_values_meth)
df_unmeth_w_types = add_probetypes(df_450, df_values_unmeth)

### QUANTILE NORMALIZATION ###

import qnorm

df_values_meth = df_values_meth.head(100)

fix, ax = plt.subplots(3)

ax[0] = df_values_meth.plot.density(linewidth=1)
ax[0] = plt.title("Raw density plot")

df_log = df_values_meth
df_log = df_log.applymap(dataframe_log)

ax[1] = df_log.plot.density(linewidth=1)
ax[1] = plt.title("Log density plot")

df_qn = qnorm.quantile_normalize(df_values_meth, axis=1)

ax[2] = df_qn.plot.density(linewidth=1)
ax[2] = plt.title("Quantile normalized plot")

plt.show()


print("All done")