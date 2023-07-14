import numpy as np
import pandas as pd
import qnorm as qn
import scipy.stats as stats
import streamlit

import prepare_data as prep


@streamlit.cache_data
def mean_normalization(df):
	sample_list = df.columns.values.tolist()
	norm_df = (df[sample_list] - df[sample_list].mean()) / df[sample_list].std()

	return norm_df


@streamlit.cache_data
def min_max_normalization(df):
	sample_list = df.columns.values.tolist()
	mm_df = (df[sample_list] - df[sample_list].min()) / (
			df[sample_list].max() - df[sample_list].min())
	return mm_df


# QN muss auf intensities durchgeführt werden
@streamlit.cache_data
def quantile_normalization(reference):
	df_meth, df_unmeth = prep.get_dataframe(False)

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


@streamlit.cache_data
def bmiq_unmethylated(df_t2_unmethylated, unmethylated_probes,
                      mean_type2_unmethylated, df_t2_parameters,
                      df_t1_parameters):
	'''STEP 2'''
	'''for type2 probes with U-state: transform their probabilities of 
	belonging to the U-state to quantiles using the inverse of the cumulative
	beta-distribution with beta parameters (aU2, bU2)'''
	'''Transform: Inverse Transform sampling?'''
	u2l_list = []
	u2r_list = []
	q_u_list = []
	for probe in unmethylated_probes:
		value = df_t2_unmethylated.loc[probe]
		if value <= mean_type2_unmethylated:
			u2l_list.append(probe)
		else:
			u2r_list.append(probe)
		if probe in u2l_list:
			rb = np.array(df_t2_unmethylated[probe])
			'''p: probability of probe belonging to the U state'''
			p_u = stats.beta.cdf(rb, df_t2_parameters['a']['U'],
			                     df_t2_parameters['b']['U'])
			q_u = stats.beta.ppf(p_u, df_t1_parameters['a']['U'],
			                     df_t1_parameters['b']['U'])
			q_u_list.append(q_u)
		else:
			rb = np.array(df_t2_unmethylated[probe])
			'''p: probability of probe belonging to the U state'''
			p_u = 1 - stats.beta.cdf(rb, df_t2_parameters['a']['U'],
			                         df_t2_parameters['b']['U'])
			q_u = 1 - stats.beta.ppf(p_u, df_t1_parameters['a']['U'],
			                         df_t1_parameters['b']['U'])
			q_u_list.append(q_u)
	df = pd.DataFrame(q_u_list, index=unmethylated_probes)
	return df, q_u_list


@streamlit.cache_data
def bmiq_methylated(df_t2_methylated, methylated_probes, mean_type2_methylated,
                    df_t2_parameters, df_t1_parameters):
	'''STEP 3'''
	'''for type2 probes with M-state: transform their probabilities of 
	belonging to the M-state to quantiles using the inverse of the cumulative
	beta-distribution with beta parameters (aM2, bM2)'''
	m2l_list = []
	m2r_list = []
	q_m_list = []
	for probe in methylated_probes:
		value = df_t2_methylated.loc[probe]
		if value <= mean_type2_methylated:
			m2l_list.append(probe)
		else:
			m2r_list.append(probe)
		if probe in m2l_list:
			rb = np.array(df_t2_methylated[probe])
			'''p: probability of probe belonging to the U state'''
			p_m = stats.beta.cdf(rb, df_t2_parameters['a']['M'],
			                     df_t2_parameters['b']['M'])
			if p_m < 1.0e-100:
				p_m = 1.0e-100
			q_m = stats.beta.ppf(p_m, df_t1_parameters['a']['M'],
			                     df_t1_parameters['b']['M'], loc=0, scale=1)
			q_m_list.append(q_m)
		else:
			rb = np.array(df_t2_methylated[probe])
			'''p: probability of probe belonging to the U state'''
			p_m = 1 - stats.beta.cdf(rb, df_t2_parameters['a']['M'],
			                         df_t2_parameters['b']['M'])
			if p_m < 1.0e-100:
				p_m = 1.0e-100
			q_m = 1 - stats.beta.ppf(p_m, df_t1_parameters['a']['M'],
			                         df_t1_parameters['b']['M'], loc=0,
			                         scale=1)
			q_m_list.append(q_m)

	df = pd.DataFrame(q_m_list, index=methylated_probes)
	return df, q_m_list


@streamlit.cache_data
def bmiq_hemimethylated(df_t2_hemimethylated, df_t2_unmethylated,
                        df_t2_methylated, eta_u_list,
                        eta_m_list, hemimethylated_probes):
	'''STEP 4'''
	'''for type2 probes with H-state: perform a dilation (scale) 
	transformation to "fit" the data into the "gap" with endpoints defined by
	max(eta2U) and min(eta2M)'''
	eta_2_H_list = []
	maxH = df_t2_hemimethylated.max()
	minH = df_t2_hemimethylated.min()
	minM = df_t2_methylated.min()
	maxU = df_t2_unmethylated.max()
	delta_beta_H = maxH - minH
	deltaUH = minH - maxU
	deltaHM = minM - maxH
	nminH = max(eta_u_list) - deltaUH
	nmaxH = min(eta_m_list) - deltaHM

	delta_eta_H = nmaxH - nminH
	dilation_factor = delta_eta_H / delta_beta_H
	for probe in hemimethylated_probes:
		eta_2_H = nminH + dilation_factor * (df_t2_hemimethylated[probe] - minH)
		eta_2_H_list.append(eta_2_H)
	df = pd.DataFrame(eta_2_H_list, index=hemimethylated_probes)
	return df


@streamlit.cache_data
def bmiq(df_beta):
	'''STEP 1'''
	'''Fitting of a three-state (unmethylated-U, hemimethylated-H, fully
	methylated-M) beta mixture model to the type1 and type2 probes
	separately. For sake of convenience we refer to intermediate allelic
	methylation as hemimethylation even though hemimethylation
	is most often used in the context of strand-specific methylation.
	-> Realized with betamix (Schroeder, Rahmann)'''
	df_beta_t1, df_beta_t2 = prep.split_types(df_beta)

	'''list with all names of probes with type 2'''
	sample_list = df_beta.columns.values.tolist()[1:]
	probe_list_t2 = df_beta.loc[df_beta["type"] == "II"].index.values.tolist()

	'''dataframe with all classes to all samples and type 2 probes'''
	df_classes_t2 = None
	df_bmiq = None
	for n, sample in enumerate(sample_list):
		'''parameters per sample (per type)'''
		df_t1_parameters = prep.get_parameters(
			"C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
			"\\betamix-results\\type1_probes-est.h5", n)
		df_t2_parameters = prep.get_parameters(
			"C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
			"\\betamix-results\\type2_probes-est.h5", n)
		'''classes per sample (type 2)'''
		df_classes = prep.get_classes(
			"C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
			"\\betamix-results\\type2_probes-eval.h5", probe_list_t2, n)
		'''list with classes for type 2 probes (per sample)'''
		a = df_classes["sample"].to_numpy()

		'''extends df_classes_t2 for each sample'''
		if df_classes_t2 is None:
			df_classes_t2 = pd.DataFrame(a, index=probe_list_t2,
			                             columns=[sample])
		else:
			df_classes_t2[sample] = a

		'''lists with probes per sample'''
		unmethylated_probes = df_classes_t2.loc[df_classes_t2[sample] ==
		                                        0].index.values.tolist()
		methylated_probes = df_classes_t2.loc[df_classes_t2[sample] ==
		                                      2].index.values.tolist()
		hemimethylated_probes = df_classes_t2.loc[df_classes_t2[sample] ==
		                                          1].index.values.tolist()
		'''dataframes with corresponding probes'''
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

		mean_type2_unmethylated = df_t2_parameters['a']['U'] / (
				df_t2_parameters[
					'a']['U'] +
				df_t2_parameters[
					'b']['U'])
		mean_type2_methylated = df_t2_parameters['a']['M'] / (df_t2_parameters[
			                                                      'a']['M'] +
		                                                      df_t2_parameters[

			                                                      'b']['M'])
		# not needed
		# mean_type2_hemimethylated = df_t2_parameters['a']['H'] / (
		# 		df_t2_parameters[
		# 			'a']['H'] +
		# 		df_t2_parameters[
		# 			'b']['H'])

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
			df_unmethylated_values, eta_u_list = (bmiq_unmethylated(
				df_t2_unmethylated,
				unmethylated_probes,
				mean_type2_unmethylated,
				df_t2_parameters,
				df_t1_parameters))
			df_methylated_values, eta_m_list = bmiq_methylated(df_t2_methylated,
			                                                   methylated_probes,
			                                                   mean_type2_methylated,
			                                                   df_t2_parameters,
			                                                   df_t1_parameters)

			df_hemimethylated_values = bmiq_hemimethylated(df_t2_hemimethylated,
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
			- append du df_bmiq'''
			frames = [df_unmethylated_values, df_methylated_values,
			          df_hemimethylated_values]
			normalized_values = pd.concat(frames)
			if df_bmiq is None:
				df_bmiq = pd.DataFrame(normalized_values.to_numpy(),
				                       index=probe_list_t2, columns=[sample])
			else:
				df_bmiq[sample] = normalized_values
	return df_bmiq
