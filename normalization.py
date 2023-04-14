import numpy as np
import pandas as pd
import qnorm as qn
import scipy.stats as stats

import prepare_data as prep


def mean_normalization(df):
	sample_list = df.columns.values.tolist()[2:]
	norm_df = (df[sample_list] - df[sample_list].mean()) / df[sample_list].std()

	return norm_df


def min_max_normalization(df):
	sample_list = df.columns.values.tolist()[2:]
	mm_df = (df[sample_list] - df[sample_list].min()) / (
			df[sample_list].max() - df[sample_list].min())

	return mm_df


def quantile_normaliziation(df, reference):
	df_qn = qn.quantile_normalize(df, target=df[reference])

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


def likelihood(params, data):
	return stats.norm.logpdf(data, loc=params[0], scale=params[1]).sum()


def neglikelihood(params, data):
	return -1 * likelihood(params, data)


def bmiq_unmethylated(df_beta_t2, probe_list_t2, sample_list,
                      mean_type2_unmethylated, df_t2_parameters,
                      df_t1_parameters):
	'''STEP 2'''
	'''for type2 probes with U-state: transform their probabilities of 
	belonging to the U-state to quantiles using the inverse of the cumulative
	beta-distribution with beta parameters (aU2, bU2)'''
	'''Transform: Inverse Transform sampling?'''
	u2l_list = []
	u2r_list = []
	df_bmiq_unmethylated = None
	for sample in sample_list:
		q_u_list = []
		for probe in probe_list_t2:
			value = df_beta_t2.loc[probe][sample]
			if value <= mean_type2_unmethylated:
				u2l_list.append(probe)
			else:
				u2r_list.append(probe)
			if probe in u2l_list:
				rb = np.array(df_beta_t2[sample][probe])
				'''p: probability of probe belonging to the U state'''
				p_u = stats.beta.cdf(rb, df_t2_parameters['a']['U'],
				                     df_t2_parameters['b']['U'])
				q_u = stats.beta.ppf(p_u, df_t1_parameters['a']['U'],
				                     df_t1_parameters['b']['U'])
				q_u_list.append(q_u)
			else:
				rb = np.array(df_beta_t2[sample][probe])
				'''p: probability of probe belonging to the U state'''
				p_u = 1 - stats.beta.cdf(rb, df_t2_parameters['a']['U'],
				                         df_t2_parameters['b']['U'])
				q_u = 1 - stats.beta.ppf(p_u, df_t1_parameters['a']['U'],
				                         df_t1_parameters['b']['U'])
				q_u_list.append(q_u)
		if df_bmiq_unmethylated is None:
			df_bmiq_unmethylated = pd.DataFrame(q_u_list, index=probe_list_t2,
			                                    columns=[
				                                    sample])
		else:
			df_bmiq_unmethylated[sample] = q_u_list
	return df_bmiq_unmethylated


def bmiq_methylated(df_beta_t2, probe_list_t2, sample_list,
                    mean_type2_methylated,
                    df_t2_parameters, df_t1_parameters):
	'''STEP 3'''
	'''for type2 probes with M-state: transform their probabilities of 
	belonging to the M-state to quantiles using the inverse of the cumulative
	beta-distribution with beta parameters (aM2, bM2)'''
	m2l_list = []
	m2r_list = []
	df_bmiq_methylated = None
	for sample in sample_list:
		q_m_list = []
		for probe in probe_list_t2:
			value = df_beta_t2.loc[probe][sample]
			if value <= mean_type2_methylated:
				m2l_list.append(probe)
			else:
				m2r_list.append(probe)
			if probe in m2l_list:
				rb = np.array(df_beta_t2[sample][probe])
				'''p: probability of probe belonging to the U state'''
				p_m = stats.beta.cdf(rb, df_t2_parameters['a']['M'],
				                     df_t2_parameters['b']['M'])
				if p_m < 1.0e-100:
					p_m = 1.0e-100
				# print("Sample: ", sample, "Probe: ", probe, "pM: ", p_m)
				q_m = stats.beta.ppf(p_m, df_t1_parameters['a']['M'],
				                     df_t1_parameters['b']['M'], loc=0, scale=1)
				# print("p_M: ", p_m, "q_M: ", q_m)
				q_m_list.append(q_m)
			else:
				rb = np.array(df_beta_t2[sample][probe])
				'''p: probability of probe belonging to the U state'''
				p_m = 1 - stats.beta.cdf(rb, df_t2_parameters['a']['M'],
				                         df_t2_parameters['b']['M'])
				if p_m < 1.0e-100:
					p_m = 1.0e-100
				q_m = 1 - stats.beta.ppf(p_m, df_t1_parameters['a']['M'],
				                         df_t1_parameters['b']['M'], loc=0,
				                         scale=1)
				q_m_list.append(q_m)

		if df_bmiq_methylated is None:
			df_bmiq_methylated = pd.DataFrame(q_m_list, index=probe_list_t2,
			                                  columns=[
				                                  sample])
		else:
			df_bmiq_methylated[sample] = q_m_list
	return df_bmiq_methylated


def bmiq_hemimethylated(df_t2_hemimethylated,df_t2_unmethylated,
                        df_t2_methylated,
                        probe_list_t2_h,
                        sample_list,
                        mean_type2_methylated):
	'''STEP 4'''
	'''for type2 probes with H-state: perform a dilation (scale) 
	transformation to "fit" the data into the "gap" with endpoints defined by
	max(eta2U) and min(eta2M)'''
	df_bmiq_methylated = None
	for sample in sample_list:
		q_m_list = []
		maxH = df_t2_hemimethylated[sample].max()
		minH = df_t2_hemimethylated[sample].min()
		minM = df_t2_methylated[sample].min()
		maxU = df_t2_unmethylated[sample].max()
		deltaH = maxH-minH
		deltaUH = minH - maxU
		deltaHM = minM - maxH


	#return df_bmiq_hemimethylated


def bmiq():
	'''STEP 1'''
	'''Fitting of a three-state (unmethylated-U, hemimethylated-H, fully
	methylated-M) beta mixture model to the type1 and type2 probes
	separately. For sake of convenience we refer to intermediate allelic
	methylation as hemimethylation even though hemimethylation
	is most often used in the context of strand-specific methylation.
	-> Realized with betamix (Schroeder, Rahmann)'''
	df_meth, df_unmeth = prep.get_values_as_dataframe_w_types()

	'''Prepping for betamix'''
	# df_beta = prep.beta_value(df_meth, df_unmeth, 100)
	# df_beta.to_csv('df_beta.csv', sep='\t')
	df_beta = pd.read_csv('df_beta.csv', sep='\t', index_col=0)
	df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
	sample_list = df_beta_t2.columns.values.tolist()
	# prep.df_to_h5(df_beta_t1, "type1_probes")
	# prep.df_to_h5(df_beta_t2, "type2_probes")
	probe_list_t2 = df_meth.loc[df_meth["type"] == "II"].index.values.tolist()

	df_classes_t2 = None # all classes to all samples
	for sample in sample_list:
		n = 0

		df_t1_parameters = prep.get_parameters(
			"C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
			"\\betamix-results\\type1_probes-est.h5", n)
		df_t2_parameters = prep.get_parameters(
			"C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
			"\\betamix-results\\type2_probes-est.h5", n)
		df_classes = prep.get_classes(
			"C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
			"\\betamix-results\\type2_probes-eval.h5", probe_list_t2, n)
		a = []
		a = df_classes["sample"].to_numpy()
		if df_classes_t2 is None:
			df_classes_t2 = pd.DataFrame(a, index=probe_list_t2,
			                                  columns=[sample])
		else:
			df_classes_t2[sample] = a
		n = n + 1

	#TODO: Redo everything with fitting classes.
	df_t2_unmethylated = df_beta_t2.loc[df_beta_t2["class"] == 0]
	df_t2_methylated = df_beta_t2.loc[df_beta_t2["class"] == 2]
	df_t2_hemimethylated = df_beta_t2.loc[df_beta_t2["class"] == 1]

	probe_list_t2_u = df_t2_unmethylated.index.values.tolist()
	probe_list_t2_m = df_t2_methylated.index.values.tolist()
	probe_list_t2_h = df_t2_hemimethylated.index.values.tolist()
	'''Prep Step 2 & 3'''
	'''Calculate mean_type2_unmethylated, then let U2L(U2R) = set of U2 probes 
	with beta-values smaller(larger) than mean_type2_unmethylated '''
	''' Calculate mean_type2_methylated, then let M2L(M2R) = set of M2 probes with 
	beta-values smaller(larger) than mean_type2_methylated'''

	mean_type2_unmethylated = df_t2_parameters['a']['U'] / (df_t2_parameters[
		                                                        'a']['U'] +
	                                                        df_t2_parameters[
		                                                        'b']['U'])
	mean_type2_methylated = df_t2_parameters['a']['M'] / (df_t2_parameters[
		                                                      'a']['M'] +
	                                                      df_t2_parameters[
		                                                      'b']['M'])
	mean_type2_hemimethylated = df_t2_parameters['a']['H'] / (df_t2_parameters[
		                                                          'a']['H'] +
	                                                          df_t2_parameters[
		                                                          'b']['H'])

	df_bmiq_unmethylated = bmiq_unmethylated(df_t2_unmethylated,
	                                         probe_list_t2_u,
	                                         sample_list,
	                                         mean_type2_unmethylated,
	                                         df_t2_parameters, df_t1_parameters)
	df_bmiq_methylated = bmiq_methylated(df_t2_methylated, probe_list_t2_m,
	                                     sample_list,
	                                     mean_type2_methylated,
	                                     df_t2_parameters, df_t1_parameters)
	df_bmiq_hemimethylated = bmiq_hemimethylated(df_t2_hemimethylated,
	                                             df_t2_unmethylated,
	                                             df_t2_methylated,
	                                             probe_list_t2_h,
	                                             sample_list,
	                                             mean_type2_methylated)

	return df_bmiq_unmethylated, df_bmiq_methylated#, df_bmiq_hemimethylated
