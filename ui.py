import shutil
import os
import glob

import streamlit as st
import time

import plotting as plot
import user

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
	with st.sidebar:
		st.write("What normalizations do you want to see?")
		b_mean_norm = st.checkbox('Mean Normalization')
		b_qn_norm = st.checkbox('Quantile Normalization')
		b_minmax_norm = st.checkbox('Min-Max Normalization')
		b_bmiq_norm = st.checkbox('Beta-Mixture-Quantile Normalization')
		st.info('Decide here if you want to see the Dataframes.')
		show_beta_df = st.checkbox('Beta Values as Dataframe')
		show_beta_plot = st.checkbox('Beta Values as Plot')
		show_types_df = st.checkbox('Beta Values as Dataframe (Type sorted)')
		show_bmiq_df = st.checkbox('BMIQ-Normalized Dataframe')
		show_boxplot_df = st.checkbox('Boxplots for BMIQ')
		show_qn_df = st.checkbox('Quantile Normalized Dataframe')
		show_mean_df = st.checkbox('Mean Normalized Dataframe')
		show_minmax_df = st.checkbox('Minmax Normalized Dataframe')




	if st.button("Start Normalization", key="start"):

		# ORIGINAL DATA
		df_meth, df_unmeth = plot.logged_plot()
		user.download_df_as_csv(df_meth, 'download/log-meth.csv')
		user.download_df_as_csv(df_unmeth, 'download/log-unmeth.csv')
		# NORMALIZATION
		if b_mean_norm:
			df_meth, df_unmeth = plot.mean_normalized_plots(
				show_mean_df)
			user.download_df_as_csv(df_meth, 'download/meannorm-meth.csv')
			user.download_df_as_csv(df_unmeth, 'download/meannorm-unmeth.csv')
		if b_qn_norm:
			df_meth, df_unmeth = plot.quantile_normalized_plots(show_qn_df)
			user.download_df_as_csv(df_meth, 'download/qnnorm-meth.csv')
			user.download_df_as_csv(df_unmeth, 'download/qnnorm-unmeth.csv')
		if b_minmax_norm:
			df_meth, df_unmeth = plot.min_max_normalized_plots(
				show_minmax_df)
			user.download_df_as_csv(df_meth, 'download/minmaxnorm-meth.csv')
			user.download_df_as_csv(df_unmeth, 'download/minmaxnorm-unmeth.csv')
		if show_beta_plot or b_bmiq_norm:
			df = plot.everything_beta(show_beta_df, show_types_df,
			                          show_bmiq_df,
			                          show_boxplot_df,
			                          b_bmiq_norm, show_beta_plot)
			user.download_df_as_csv(df, 'download/bmiq.csv')

