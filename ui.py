import streamlit as st
import pandas as pd

import plotting as plot
import prepare_data as prep

""" UI.py
This class is dealing with the Streamlit Structure. Single Figures etc. are 
in Plotting.py.
"""

def _streamlit_config():
	st.set_page_config(page_title="450K normalization",
	                   page_icon=":chart_with_upwards_trend:",
	                   layout="wide")
	st.set_option('deprecation.showPyplotGlobalUse', False)


def make_header():
	"""Builds Header"""
	_streamlit_config()
	# HEADER #
	with st.container():  # for wrapping contents
		st.header("Normalization of Illumina HumanMethylation 450K Beadchips")
		st.subheader(
			"Methylation normalization because technical variability sucks.")

	# TABLE #
	with st.container():
		left_column, right_column = st.columns(2)
		with left_column:
			st.header("Who I am")
			st.write("""
            Lisa
            - Bioinformatics Student (Bachelor of Science)
            - struggling my way through university
            """
			         )
		with right_column:
			st.header("My University")
			st.write("""
            Saarland University
            - Center of Bioinformatics
            - Chair of Algorithmic Bioinformatics
            - university in a land just used for size comparisons
            """
			         )


def make_plots():
	"""Cares about all the plots and Dataframe-Views"""

	with st.container():
		left_column, right_column = st.columns(2)
		with left_column:
			st.header("Methylated Plots")
			# TODO: DONE
			plot.methylated_plots()

		with right_column:
			st.header("Unmethylated Plots")
			# TODO: DONE
			plot.unmethylated_plots()

	with st.container():
		plot.beta_value_plots()


# # TODO: UNDER CONSTRUCTION
# # BMIQ
#
# plot.beta_value_plots(df_beta)
#
# # trying to make data work in betamix
# df_beta.to_hdf('data.h5', key='df', mode='w')
#
# # FITTING BETA DISTRIBUTION TO OUR DATA
# m = df_beta.iloc[:, 1:66]  # Beta values in Array without classes
# # transform into numpy array (for fitting)
# naive_array = m.to_numpy()
#
# # fit: Return estimates of shape (if applicable), location, and scale parameters from data. The default
# # estimation method is Maximum Likelihood Estimation (MLE), but Method of Moments (MM) is also available.
# beta_params = stats.beta.fit(naive_array[0])
# # returns probability density function with given parameters (here, parameters from our array)
# st.write("Beta Distribution PDF")
# df_beta_pdf = pd.DataFrame(stats.beta.pdf(naive_array, beta_params[0], beta_params[1],
#                                           beta_params[2], beta_params[3]))
# fig_beta = plot.density_plot(df_beta_pdf, "beta pdf", -5, 10)
# st.pyplot(fig_beta)
#
# # Teschendorff
# # Fitting of a three-state (unmethylated-U, hemimethylated-H, fully
# # methylated-M) beta mixture model to the type1 and type2 probes
# # separately.
#
# # TODO: dataframe df_beta_type1 needs index and column-header
# st.header("Type I Fitting")
# type1 = df_beta[df_beta["type"] == 'I'] #full dataframe with just type I probes
# st.write("Type1: ", type1)
# type1_array = type1.iloc[:, 1:66].to_numpy() #full dataframe as array for further calculations
# st.write("Type1_array: ", type1_array)
# type1_fit = stats.beta.fit(type1_array[0]) #parameters fitting to found beta distribution
# st.write("Type1_fit: ", type1_fit)
# df_beta_type1 = pd.DataFrame(stats.beta.pdf(type1_array, type1_fit[0], type1_fit[1],
#                                             type1_fit[2], type1_fit[3]), columns=sample_list_meth)
# st.write(df_beta_type1)
# fig_beta_type1 = plot.density_plot(df_beta_type1, "beta type 1 pdf", -5, 10)
# st.write("Beta distribution Type 1")
# st.pyplot(fig_beta_type1)
#
# # TODO: dataframe df_beta_type2 needs index and column-header
# st.header("Type II Fitting")
# type2 = df_beta[df_beta["type"] == 'II']
# type2_array = type2.iloc[:, 1:66].to_numpy()
# type2_fit = stats.beta.fit(type2_array[0])
# df_beta_type2 = pd.DataFrame(stats.beta.pdf(type2_array, type2_fit[0], type2_fit[1],
#                                             type2_fit[2], type2_fit[3]))
# st.write(df_beta_type2)
# fig_beta_type2 = plot.density_plot(df_beta_type2, "beta type 2 pdf", -5, 10)
# st.write("Beta distribution Type 2")
# st.pyplot(fig_beta_type2)
#
# # NAIVE CLASSIFICATION WITH CUTOFF POINT
# st.header("Naive classification")
#
# df_states = norm.assign_probes_to_state(df_beta)
# st.write("Naive Classification of states")
# df_t = df_states.transpose()
# df_t.drop('type', inplace=True)
# st.write(df_t)
#
# # TODO: Does not work
# # Trying to get the sample proportions for initial naive values.
# sample_proportions = pd.DataFrame()
#
# first_time = True
# # returns normalized values (probabilities for states)
# for column in df_t:
#     step = df_t[column].value_counts(normalize=True)
#     step.to_frame(name=step.iloc[0])
#
# # def neg_likelihood(params, data):
# #    return -1*scipy.beta.logpdf(data, loc=params[0], scale=params[1]).sum()
#
# # array = np.array([50, 10])
# # result = optimize.minimize(neg_likelihood, array, args=data)
# # st.write(result)


# # ------------------------------------------------------------------------GET DATA-#
# df_meth = pd.read_csv('resources/short_methylated_w_types.csv', sep='\t')
# df_unmeth = pd.read_csv('resources/short_unmethylated_w_types.csv', sep='\t')
# sample_list_meth = df_meth.columns.values.tolist()[2:]
# sample_list_unmeth = df_meth.columns.values.tolist()[2:]
# df_unmeth[sample_list_unmeth] = df_unmeth[sample_list_unmeth].astype(int)
# df_meth[sample_list_meth] = df_meth[sample_list_meth].apply(lambda x: x + 1)
# df_unmeth[sample_list_unmeth] = df_unmeth[sample_list_unmeth].apply(
# 	lambda x: x + 1)
# # ---------------------------------------------------------------------------------#

