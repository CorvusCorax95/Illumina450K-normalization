import math
import pandas as pd
import numpy as np

"""PREPARE DATA
This file transforms all the nasty inputs we get to pretty little dataframes.
"""


def _probetype_to_dataframe(input):
	"""helping dataframe for probe-types matched to probes"""
	# make dataframe from csv-file
	df = pd.read_csv(input, sep=',')

	# set probes as index

	df_i = df.set_index('probe')
	return df_i


def _add_probetypes(df_dest, df_target):
	"""adding probetypes (I and II) to my dataframe	to distinguish between
	Infinium I and Infinium II probes"""

	df_merged = pd.merge(df_dest, df_target, on="probe")

	return df_merged


def dataframe_log(x):
	if (isinstance(x, int) and x != 0):
		math.log(x, math.e)
	return x


def get_values_as_dataframe_w_types():
	"""rebuild illumina probe tables to dataframe with just probes and types"""

	df_450 = _probetype_to_dataframe('resources/illumina-table-450.csv')

	"""ORIGINAL FILES"""
	# # convert methylation-tables to dataframes
	# df_values_meth = pd.read_csv('methylation-rb/methylation_meth.wide.tsv', sep='\t')
	# df_values_unmeth = pd.read_csv('methylation-rb/methylation_unmeth.wide.tsv', sep='\t')

	"""SHORT FILES"""
	df_values_meth = pd.read_csv('data/short_meth.tsv', sep='\t').set_index(
		'probe')
	df_values_unmeth = pd.read_csv('data/short_unmeth.tsv', sep='\t').set_index(
		'probe')

	# match probes to type in methylation file
	# returns new dataframe including probe types
	df_meth_w_types = _add_probetypes(df_450, df_values_meth)
	df_unmeth_w_types = _add_probetypes(df_450, df_values_unmeth)

	df_meth_w_types.to_csv('short_methylated_w_types.csv', sep='\t')
	df_unmeth_w_types.to_csv('short_unmethylated_w_types.csv', sep='\t')

	return df_meth_w_types, df_unmeth_w_types


# applying log to all values in the dataframe
# stores dataframe as csv in new file
def log_data():
	"""Uses log on every value in the dataframe."""
	df_meth, df_unmeth = get_values_as_dataframe_w_types()
	# log values of dataframe
	index = df_meth.index

	sample_list = df_meth.columns.values.tolist()[2:]
	df_meth.replace(0, 0.01, inplace=True)
	df_log_meth = np.log2(df_meth[sample_list])

	# find and remove inf values after log
	df_log_meth.replace([np.inf, -np.inf], np.nan, inplace=True)
	df_log_meth = df_log_meth.dropna()

	sample_list = df_unmeth.columns.values.tolist()[2:]
	df_unmeth.replace(0, 0.01, inplace=True)
	df_log_unmeth = np.log2(df_unmeth[sample_list])

	# find and remove inf values after log
	df_log_unmeth.replace([np.inf, -np.inf], np.nan, inplace=True)
	df_log_unmeth = df_log_unmeth.dropna()

	return df_log_meth, df_log_unmeth


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
			df[x][y] = df_meth[x][y] / (
					df_meth[x][y] + df_unmeth[x][y] + offset)
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
