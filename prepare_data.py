import pandas as pd
import numpy as np
import h5py as h5

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


def add_probetypes(df_target):
	"""adding probetypes (I and II) to my dataframe	to distinguish between
	Infinium I and Infinium II probes"""
	# df_450 = _probetype_to_dataframe('resources/illumina-table-450.csv')
	df_450 = _probetype_to_dataframe('illumina-table-450.csv')
	df_target.index.name = 'probe'
	df_merged = pd.merge(df_450, df_target, on="probe")

	return df_merged


def get_dataframe(with_types):
	"""Rebuild illumina probe tables to dataframe with just probes and types"""

	"""ORIGINAL FILES"""
	# convert methylation-tables to dataframes
	# df_values_meth = pd.read_csv('data/methylation_meth.wide.tsv',
	#                             sep='\t').set_index('probe')
	# df_values_unmeth = pd.read_csv('data/methylation_unmeth.wide.tsv',
	# sep='\t').set_index('probe')

	# # """SHORT FILES"""
	# df_values_meth = pd.read_csv('data/short_meth.tsv', sep='\t').set_index(
	# 	'probe')
	# df_values_unmeth = pd.read_csv('data/short_unmeth.tsv',
	#                                sep='\t').set_index('probe')

	# """SHORT FILES""" needed for thesis plots
	df_meth = pd.read_csv('short_meth.tsv', sep='\t').set_index('probe')
	df_unmeth = pd.read_csv('short_unmeth.tsv', sep='\t').set_index('probe')

	"""	Matches probe types (Infinium I or Infinium II) from the df_450 file (
	file from Illumina with	mapping of probename to probetype) """
	if with_types:
		df_meth = add_probetypes(df_meth)
		df_unmeth = add_probetypes(df_unmeth)

	"""Saves in case of problems"""
	# df_meth_w_types.to_csv('methylated_w_types.csv', sep='\t')
	# df_unmeth_w_types.to_csv('unmethylated_w_types.csv', sep='\t')
	return df_meth, df_unmeth


# applying log to all values in the dataframe
# stores dataframe as csv in new file
def log_data():
	"""Logs on every value in the dataframe."""
	df_meth, df_unmeth = get_dataframe(True)

	# log values of dataframe
	sample_list = df_meth.columns.values.tolist()[1:]

	"""Methylated File"""
	df_meth.replace(0, 0.01, inplace=True)
	df_log_meth = np.log2(df_meth[sample_list])
	# find and remove inf values after log
	df_log_meth.replace([np.inf, -np.inf], np.nan, inplace=True)
	df_log_meth = df_log_meth.dropna()

	"""Unmethylated File"""
	df_unmeth.replace(0, 0.01, inplace=True)
	df_log_unmeth = np.log2(df_unmeth[sample_list])
	# find and remove inf values after log
	df_log_unmeth.replace([np.inf, -np.inf], np.nan, inplace=True)
	df_log_unmeth = df_log_unmeth.dropna()
	return df_log_meth, df_log_unmeth


def output_measures(df):
	""""Useful output measures."""
	sample_list = df.columns.values.tolist()
	means = []
	for x in sample_list:
		means.append(df[x].mean())
	return np.average(means), np.var(means), np.std(means)


def beta_value(df_meth, df_unmeth, offset):
	"""Makes beta-values from logged values.
	beta-value: methylated / methylated + unmethylated + 100"""

	sample_list = df_meth.columns.values.tolist()[1:]
	probe_list = df_meth.index.values.tolist()
	df = pd.DataFrame(index=probe_list, columns=sample_list)

	for sample in sample_list:
		for probe in probe_list:
			df[sample][probe] = (df_meth[sample][probe] + (offset / 2)) / (
					df_meth[sample][probe] + df_unmeth[sample][
				probe] + offset)
	return df


def m_value(df_meth, df_unmeth):
	"""Calculates M-Value. Different approach than beta-value."""
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


def split_types(df):
	"""BMIQ Normalization needs separate dataframes per type so this is the
	function that splits them."""
	df_t1 = df.loc[df["type"] == 'I']
	df_t2 = df.loc[df["type"] == 'II']
	del df_t1[df_t1.columns[0]]
	del df_t2[df_t2.columns[0]]
	return df_t1, df_t2


def get_parameters(path, n):
	'''After using the betamix-tool we need to get the esimated parameters
	from the h5 file provided.'''
	num = str(n)
	f = h5.File(path, "r")
	grp = f['estimation']
	subgrp = grp[num]
	data = subgrp["ab"]
	return pd.DataFrame(data, index=["U", "H", "M"], columns=["a", "b"])


def get_classes(path, probe_list, n):
	'''gets the class to probe and returns dataframe'''
	'''0: unmethylated, 1: hemimethylated, 2: methylated'''
	num = str(n)
	f = h5.File(path, "r")
	grp = f['evaluation']
	subgrp = grp[num]
	subsubgrp = subgrp['mix']
	data = subsubgrp['classes']
	return pd.DataFrame(data, index=probe_list, columns=["sample"])

def df_to_h5(df, filename):
	"""As preparation to use the betamix-tool by Schröder, Rahmann the
	dataframe is converted to a fitting h5 file."""
	# delete existing file first by hand!
	sample_list = df.columns.values.tolist()[1:]
	df.reset_index(drop=True)
	filename = filename + ".h5"
	hf = h5.File(filename, "w")
	group = hf.create_group("data")
	i = 0
	for sample in sample_list:
		array = np.array(df[sample], dtype=np.float64)
		group.create_dataset(str(i), data=np.sort(array))
		i = i + 1
	hf.close()
	print("Done converting.")