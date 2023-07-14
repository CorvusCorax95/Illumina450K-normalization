import shutil

import streamlit as st

import plotting as plot
import prepare_data as prep

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


def make_page():
	with st.sidebar:
		with st.form("sidebar"):
			b_download = st.checkbox("Do you want the possibility to download "
			                         "everything at once? Clicking the provided "
			                         "button will reload the site entirely.")
			st.write("What normalizations do you want to see?")
			b_mean = st.checkbox('Mean Normalization')
			b_qn = st.checkbox('Quantile Normalization')
			b_minmax = st.checkbox('Min-Max Normalization')
			b_bmiq = st.checkbox('Beta-Mixture-Quantile Normalization')
			b_qn_bmiq = st.checkbox('Quantile Normalization + '
			                        'Beta-Mixture-Quantile Normalization')
			st.info('Decide here if you want to see the Additional Info.')
			show_types = st.checkbox('Additional view split by types')
			show_boxplot = st.checkbox('Boxplots for BMIQ')

			submitted = st.form_submit_button("Run Normalization")

	df_beta = plot.beta_value()
	if submitted:
		make_plots(df_beta, b_mean, b_minmax, b_qn, b_bmiq, b_qn_bmiq,
		           show_types, show_boxplot)

	if b_download:
		shutil.make_archive("compressed_download", 'zip', "download")

		filename = "compressed_download.zip"
		with open(filename, "rb") as fp:
			btn = st.download_button(
				label="Download normalized values as csv",
				data=fp,
				file_name="compressed_download.zip",
				mime="application/zip",
				help="When this button is clicked, the site will reload.")


def make_plots(df_beta, b_mean, b_minmax, b_qn, b_bmiq, b_qn_bmiq, show_types,
               show_boxplot):
	"""Cares about all the plots and Dataframe-Views"""

	with st.spinner("Waiting for data..."):
		# ORIGINAL DATA
		plot.containerize_chart(df_beta, "Raw Beta Values")
		user.convert_df(df_beta, 'download/raw-beta-values.csv')
		df_beta_w_types = prep.add_probetypes(df_beta)
		user.convert_df(df_beta_w_types, 'download/raw-beta-values_w_types.csv')
		df_beta_t1, df_beta_t2 = prep.split_types(df_beta_w_types)
		plot.default_plots(df_beta_t1, df_beta_t2, "Beta Values")

	with st.spinner("Wait for Normalization..."):
		if b_mean:
			st.subheader("Mean Normalization")
			df_mean = plot.default_plots(DataType.MEAN, df_beta)
			user.convert_df(df_mean, 'download/meannorm.csv')
			if show_types:
				df_mean_w_types = prep.add_probetypes(df_mean)
				df_mean_t1, df_mean_t2 = prep.split_types(df_mean_w_types)
				plot.default_plots(df_mean_t1, df_mean_t2, "Mean Normalized")
	with st.spinner("Wait for Normalization..."):
		if b_minmax:
			st.subheader("Minmax Normalization")
			df_minmax = plot.default_plots(DataType.MINMAX, df_beta)
			user.convert_df(df_minmax, 'download/minmaxnorm.csv')
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
			if show_types:
				df_qn_w_types = prep.add_probetypes(df_qn)
				df_qn_t1, df_qn_t2 = prep.split_types(df_qn_w_types)
				plot.default_plots(df_qn_t1, df_qn_t2, "Quantile Normalized")
	with st.spinner("Wait for Normalization..."):
		if b_bmiq:
			st.subheader("BMIQ Normalization")
			df = plot.bmiq_plot(show_boxplot, df_beta)
			user.convert_df(df, 'download/bmiq.csv')
