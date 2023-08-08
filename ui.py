import shutil

import pandas as pd
import streamlit as st

import plotting as plot
import prepare_data as prep
import normalization as norm

import user

from type import DataType

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
		st.header("Normalization of Illumina HumanMethylation Beadchips")
		st.subheader(
			"Methylation normalization because technical variability sucks.")

	# TABLE #
	with st.expander("See information."):
		st.header("What this tool does")
		st.write("""
		This little website makes viewing your methylation data fun and 
		easy with colourful graphs. It is especially designed to view and 
		preprocess data from Illumina 450K and 850K Beadchips.
		The most interesting part should be the 
		Beta-Mixture-Quantile-Normalization.
        """
		         )
		st.header("What to keep in mind")
		st.write("""
		In the sidebar you can select which normalizations you want to 
		see, some additional info and if you probably want to download 
		all preprocessed files in one zip file.
		Htting the "Run Normalization" button (obviously) starts 
		calculating what you selected.\n
		In the sample selection underneath you can select how many 
		samples you want to view.
		Depending on how large the dataset is it can take quite a while 
		to process. Small spinners should show you that the program 
		does something and hasn't crashed.\n
		Important: To perform the Beta-Mixture-Quantile-Normalization 
		correctly you need to fit a three-state beta-mixture model to 
		your data. I used the betamix tool by Schröder, Rahmann and 
		stated the workflow in a file named "betamix-walkthrough" that 
		you can find among the other files in my git.
        """
		         )


def make_page():
	with st.sidebar:
		with st.form("sidebar"):
			b_download = st.checkbox("Download everything at once")
			st.write("What normalizations do you want to see?")
			b_mean = st.checkbox('Mean Normalization')
			b_qn = st.checkbox('Quantile Normalization')
			b_minmax = st.checkbox('Min-Max Normalization')
			b_bmiq = st.checkbox('Beta-Mixture-Quantile Normalization')
			#b_qn_bmiq = st.checkbox('Quantile Normalization + '
			#                        'Beta-Mixture-Quantile Normalization')
			b_qn_bmiq = False
			st.info('Additional Info.')
			show_types = st.checkbox('Additional view: split by types')
			show_boxplot = st.checkbox('Additional boxplots for normalizations')

			submitted = st.form_submit_button("Run Normalization")
	df_beta = plot.beta_value()

	df_sample_to_numbers = pd.DataFrame(df_beta.columns.values.tolist(),
	                                    columns=["samples"])
	options = st.multiselect('Select Samples.',
	                         df_beta.columns.values.tolist(), default=[
			"RB_E_001", "RB_E_002", "RB_E_003"], )
	st.write(df_beta[options])
	df_beta = df_beta[options]
	if submitted:
		make_plots(df_beta, b_mean, b_minmax, b_qn, b_bmiq, b_qn_bmiq,
		           show_types, show_boxplot, options, df_sample_to_numbers)
	with st.sidebar:
		if b_download:
			shutil.make_archive("compressed_download", 'zip', "download")

			filename = "compressed_download.zip"
			with open(filename, "rb") as fp:
				btn = st.download_button(
					label="Download normalized values as zip",
					data=fp,
					file_name="compressed_download.zip",
					mime="application/zip",
					help="When this button is clicked, the site will reload.")


def make_plots(df_beta, b_mean, b_minmax, b_qn, b_bmiq, b_qn_bmiq, show_types,
               show_boxplot, options, df_sample_to_numbers):
	"""Cares about all the plots and Dataframe-Views"""

	with st.spinner("Waiting for data..."):
		# ORIGINAL DATA
		plot.containerize_chart(df_beta, "Raw Beta Values")
		user.convert_df(df_beta, 'download/raw-beta-values.csv')
		df_beta_w_types = prep.add_probetypes(df_beta)
		user.convert_df(df_beta_w_types, 'download/raw-beta-values_w_types.csv')
		df_beta_t1, df_beta_t2 = prep.split_types(df_beta_w_types)
		plot.default_plots(df_beta_t1, df_beta_t2, "Beta Values")
		if show_boxplot:
			plot.boxplots(df_beta_t1)

	with st.spinner("Wait for Normalization..."):
		if b_mean:
			st.subheader("Mean Normalization")
			df_mean = plot.default_plots(DataType.MEAN, df_beta)
			user.convert_df(df_mean, 'download/meannorm.csv')
			if show_types:
				df_mean_w_types = prep.add_probetypes(df_mean)
				df_mean_t1, df_mean_t2 = prep.split_types(df_mean_w_types)
				plot.default_plots(df_mean_t1, df_mean_t2, "Mean Normalized")
			if show_boxplot:
				plot.boxplots(df_mean)
	with st.spinner("Wait for Normalization..."):
		if b_minmax:
			st.subheader("Minmax Normalization")
			df_minmax = plot.default_plots(DataType.MINMAX, df_beta)
			user.convert_df(df_minmax, 'download/minmaxnorm.csv')
			if show_boxplot:
				plot.boxplots(df_minmax)
			if show_types:
				df_minmax_w_types = prep.add_probetypes(df_minmax)
				df_minmax_t1, df_minmax_t2 = prep.split_types(df_minmax_w_types)
				plot.default_plots(df_minmax_t1, df_minmax_t2, "Minmax "
				                                               "Normalized")
	with st.spinner("Wait for Normalization..."):
		if b_qn:
			st.subheader("Quantile Normalization")
			df_qn = plot.default_plots(DataType.QN, df_beta)
			user.convert_df(df_qn, 'download/qnorm.csv')
			if show_boxplot:
				plot.boxplots(df_qn)
			if show_types:
				df_qn_w_types = prep.add_probetypes(df_qn)
				df_qn_t1, df_qn_t2 = prep.split_types(df_qn_w_types)
				plot.default_plots(df_qn_t1, df_qn_t2, "Quantile Normalized")
	with st.spinner("Wait for Normalization..."):
		if b_bmiq:
			st.subheader("BMIQ Normalization")
			st.info("The BMIQ Normalization as presented in Teschendorff et "
			        "al. (2013) only proposed a normalization for type 2 "
			        "probes. This leads to a combined histogram with one type "
			        "normalized and one type as raw beta values. To see the "
			        "normalization better, tick the box 'Additional view: "
			        "split by types'")
			df_bmiq = plot.bmiq_plot(df_beta, df_sample_to_numbers)
			user.convert_df(df_bmiq, 'download/bmiq.csv')
			if show_boxplot:
				plot.boxplots(df_bmiq)
			if show_types:
				df_bmiq_w_types = prep.add_probetypes(df_bmiq)
				df_bmiq_t1, df_bmiq_t2 = prep.split_types(df_bmiq_w_types)
				plot.default_plots(df_bmiq_t1, df_bmiq_t2, "BMIQ "
				                                           "Normalized")
	# STILL BUGGY
	# with st.spinner("Wait for Normalization..."):
	# 	if b_qn_bmiq:
	# 		st.subheader("QN.BMIQ Normalization")
	# 		df_qn_beta = norm.quantile_normalization(options, 'Median')
	# 		df_qn_bmiq = plot.bmiq_plot(df_qn_beta, df_sample_to_numbers)
	# 		user.convert_df(df_qn_bmiq, 'download/bmiq.csv')
	# 		if show_boxplot:
	# 			plot.boxplots(df_qn_bmiq)
	# 		if show_types:
	# 			df_qn_bmiq_w_types = prep.add_probetypes(df_qn_bmiq)
	# 			df_qn_bmiq_t1, df_qn_bmiq_t2 = prep.split_types(
	# 				df_qn_bmiq_w_types)
	# 			plot.default_plots(df_qn_bmiq_t1, df_qn_bmiq_t2, "QN.BMIQ "
	# 			                                           "Normalized")