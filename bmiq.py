import numpy as np
import pandas as pd
import scipy.stats as stats
import streamlit

@streamlit.cache_data
def bmiq_unmethylated_case(df_t2_unmethylated, unmethylated_probes,
                           mean_type2_unmethylated, df_t2_parameters,
                           df_t1_parameters):
	'''STEP 2'''
	'''for type2 probes with U-state: transform their probabilities of 
	belonging to the U-state to quantiles using the inverse of the cumulative
	beta-distribution with beta parameters (aU2, bU2)'''
	'''Transform: Inverse Transform sampling?'''

	# df_t2_unmethylated: DataFrame -  Dataframe with just the unmethylated
	# probes
	# unmethylated_probes: list - list of unmethylated probes (identical to
	# df_t2_unmethylated index)
	# mean_type2_unmethylated: int -  means of the estimated beta distribution
	# df_tx_parameters: parameters a and b for all three states from betamix
	u2l_list = [] #set of UII probes with beta-value smaller than mean_II_U
	u2r_list = [] #set of UII probes with beta-value larger than mean_II_U
	q_u_list = []
	x = 0
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
			if q_u > 1:
				x = x+1
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
			if q_m > 1:
				x = x+1
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


	x = 0
	list = []
	dilation_factor = delta_eta_H / delta_beta_H
	for probe in hemimethylated_probes:
		eta_2_H =  nminH + dilation_factor * (df_t2_hemimethylated[probe] -
		                                      minH)
		eta_2_H_list.append(eta_2_H)
		if eta_2_H > 1:
			x = x+1
			list.append(probe)
	print("Amount of hemimethylated probes with a problem:", x)
	df = pd.DataFrame(eta_2_H_list,
	                  index=df_t2_hemimethylated.index.values.tolist())
	return df
