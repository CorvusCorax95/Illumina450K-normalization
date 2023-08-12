import pandas as pd
import qnorm as qn
import streamlit
import numpy as np

import bmiq as b

import prepare_data as prep


@streamlit.cache_data
def mean_normalization(df):
	sample_list = df.columns.values.tolist()
	if "type" in sample_list:
		sample_list.remove("type")
	norm_df = (df[sample_list] - df[sample_list].mean()) / df[sample_list].std()
	return norm_df


@streamlit.cache_data
def min_max_normalization(df):
	sample_list = df.columns.values.tolist()
	probe_list = df.index.values.tolist()
	mm_df = None
	for sample in sample_list:
		lst = []
		max = df[sample].max()
		min = df[sample].min()
		for probe in probe_list:
			mm = (df[sample][probe] - min)/(max - min)
			lst.append(mm)
		if mm_df is None:
			mm_df = pd.DataFrame(lst, index=probe_list, columns=[sample])
		else:
			mm_df[sample] = lst
	return mm_df


# QN muss auf intensities durchgeführt werden
@streamlit.cache_data
def quantile_normalization(df, reference):
	df['Median'] = df.median(axis=1)
	df_qn = qn.quantile_normalize(df, target=df[reference])
	return df_qn


@streamlit.cache_data
def qn_for_meth(df_meth, df_unmeth, reference, sample_list):
	df_meth = df_meth[sample_list]
	df_unmeth = df_unmeth[sample_list]

	df_meth['Median'] = df_meth.median(axis=1)
	df_unmeth['Median'] = df_unmeth.median(axis=1)

	df_qn_meth = qn.quantile_normalize(df_meth, target=df_meth[reference])
	df_qn_unmeth = qn.quantile_normalize(df_unmeth, target=df_unmeth[reference])
	df_qn = prep.beta_value(df_qn_meth, df_qn_unmeth, 100)
	return df_qn


def assign_probes_to_state(df):
	df = df.applymap(lambda x: set_states(x) if type(x) == float else x)
	return df


# cutoffs from Schröder/Rahmann 2017
def set_states(x):
	# unmethylated
	if float(x) <= 0.25:
		return 'M'
	# fully-methylated
	elif float(x) >= 0.75:
		return 'U'
	# hemi-methylated
	else:
		return 'H'

def _beta_for_bmiq(df_w, df_pi, sample_list, probe_list):

	df = pd.DataFrame(np.zeros((len(probe_list), len(sample_list))),
	                  index=probe_list, columns=sample_list)
	for sample in sample_list:
		for probe in probe_list:
			df[sample][probe] = df_pi["pi"]["U"]*df_w["U"][probe] + df_pi[
				"pi"]["H"]*df_w["H"][probe] + df_pi["pi"]["M"]*df_w["M"][probe]
	return df

@streamlit.cache_data
def bmiq(df_meth, df_sample_to_numbers):
	'''STEP 1'''
	'''Fitting of a three-state (unmethylated-U, hemimethylated-H, fully
	methylated-M) beta mixture model to the type1 and type2 probes
	separately. For sake of convenience we refer to intermediate allelic
	methylation as hemimethylation even though hemimethylation
	is most often used in the context of strand-specific methylation.
	-> Realized with betamix (Schroeder, Rahmann)'''
	print(df_meth)
	df_meth_t1, df_meth_t2 = prep.split_types(df_meth)
	'''list with all names of probes with type 2'''
	sample_list = df_meth.columns.values.tolist()[1:]
	probe_list_t1 = df_meth_t1.index.values.tolist()
	probe_list_t2 = df_meth_t2.index.values.tolist()

	'''dataframe with all classes to all samples and type 2 probes'''
	df_classes_t2 = None
	df_bmiq = None
	for sample in sample_list:
		index = df_sample_to_numbers.index[df_sample_to_numbers.samples ==
		                                   sample]
		n = index[0]
		'''parameters per sample (per type)'''
		df_t1_parameters = prep.get_parameters(
			"./betamix-results/type1_probes-est.h5", n)
		df_t2_parameters = prep.get_parameters(
			"./betamix-results/type2_probes-est.h5", n)
		'''classes per sample (type 2)'''
		df_t1_pi = prep.get_pi(
			"./betamix-results/type1_probes-est.h5", n)
		df_t2_pi = prep.get_pi(
			"./betamix-results/type2_probes-est.h5", n)
		'''classes per sample (type 2)'''
		df_classes = prep.get_classes("./betamix-results/type2_probes-eval.h5"
			 ,probe_list_t2, n)
		df_t1_w = prep.get_w("./betamix-results/type1_probes-est.h5", n,
		                     sample_list, probe_list_t1)
		df_t2_w = prep.get_w("./betamix-results/type2_probes-est.h5", n,
		                     sample_list, probe_list_t2)

		'''list with classes for type 2 probes (per sample)'''
		a = df_classes["sample"].to_numpy()
		'''extends df_classes_t2 for each sample'''
		if df_classes_t2 is None:
			df_classes_t2 = pd.DataFrame(a, index=probe_list_t2,
			                             columns=[sample])
		else:
			df_classes_t2[sample] = a

		df_beta_t1 = _beta_for_bmiq(df_t1_w, df_t1_pi, sample_list,
		                           probe_list_t1)
		df_beta_t2 = _beta_for_bmiq(df_t2_w, df_t2_pi, sample_list,
		                           probe_list_t2)

		'''lists with probes per sample'''
		unmethylated_probes = df_classes_t2.loc[df_classes_t2[sample] ==
		                                        0].index.values.tolist()
		hemimethylated_probes = df_classes_t2.loc[df_classes_t2[sample] ==
		                                          1].index.values.tolist()
		methylated_probes = df_classes_t2.loc[df_classes_t2[sample] ==
		                                      2].index.values.tolist()
		'''dataframes with corresponding probes as classified in betamix'''
		df_t2_unmethylated = df_beta_t2[sample].filter(unmethylated_probes,
		                                               axis=0)
		df_t2_methylated = df_beta_t2[sample].filter(methylated_probes,
		                                             axis=0)
		df_t2_hemimethylated = df_beta_t2[sample].filter(hemimethylated_probes,
		                                                 axis=0)
		'''Prep Step 2 & 3'''
		'''Calculate mean_type2_unmethylated, then let U2L(U2R) = set of U2
		probes with beta-values smaller(larger) than mean_type2_unmethylated'''
		''' Calculate mean_type2_methylated, then let M2L(M2R) = set of M2
		probes with beta-values smaller(larger) than mean_type2_methylated'''
		# ratios are ok
		mean_type2_unmethylated = df_t2_parameters['a']['U'] / (
				df_t2_parameters[
					'a']['U'] +
				df_t2_parameters[
					'b']['U'])
		mean_type2_methylated = df_t2_parameters['a']['M'] / (df_t2_parameters[
			                                                      'a']['M'] +
		                                                      df_t2_parameters[
			                                                      'b']['M'])
		'''list with bmiq-normalized, unmethylated values of type 2 probes
		order equal to unmethylated_probes'''
		'''
		SPECIAL CASE:
		If there are no values of one of the three classes, e.g. no 
		methylated values (RB_E_027) then we cannot perform BMIQ'''
		if len(unmethylated_probes) == 0 or len(methylated_probes) == 0 or len(
				hemimethylated_probes) == 0:
			print("Sample: ", sample, "\nNormalization with BMIQ not possible.")
			pass
		else:
			df_unmethylated_values, eta_u_list = b.bmiq_unmethylated_case(
				df_t2_unmethylated,
				unmethylated_probes,
				mean_type2_unmethylated,
				df_t2_parameters,
				df_t1_parameters)
			df_methylated_values, eta_m_list = b.bmiq_methylated(
				df_t2_methylated,
				methylated_probes,
				mean_type2_methylated,
				df_t2_parameters,
				df_t1_parameters)
			df_hemimethylated_values = b.bmiq_hemimethylated(
				df_t2_hemimethylated,
				df_t2_unmethylated,
				df_t2_methylated,
				eta_u_list,
				eta_m_list,
				hemimethylated_probes)
			'''to save everything in one dataframe I have to stick everything 
			back together, such that all probes have their respective normalized 
			values.
			so i have to rebuild the columns:
			- unmethylate_values, methylated_valued and hemimethylated_values 
			sticked together
			- sort by probe name
			- append to df_bmiq'''
			frames = [df_unmethylated_values, df_methylated_values,
			          df_hemimethylated_values]
			normalized_values = pd.concat(frames)
			normalized_values.columns = [sample]
			if df_bmiq is None:
				df_bmiq = normalized_values[[sample]].copy()
			else:
				df_bmiq[sample] = normalized_values
	frames = [df_bmiq, df_beta_t1]
	df = pd.concat(frames)
	df_res = df.sort_index()
	return df_res
