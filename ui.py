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
	tab_raw, tab_norm = st.tabs(["Original Data", "Normalized Data"])
	with st.sidebar:
		st.write("What normalizations do you want to see?")
		b_mean_norm = st.checkbox('Mean Normalization')
		b_qn_norm = st.checkbox('Quantile Normalization')
		b_minmax_norm = st.checkbox('Min-Max Normalization')
		b_bmiq_norm = st.checkbox('Beta-Mixture-Quantile Normalization')
		st.info('Decide here if you want to see the Dataframes.')
		beta_df = st.checkbox('Beta Values as Dataframe')
		beta_plot = st.checkbox('Beta Values as Plot')
		types_df = st.checkbox('Beta Values as Dataframe (Type sorted)')
		bmiq_df = st.checkbox('BMIQ-Normalized Dataframe')
		boxplot_df = st.checkbox('Boxplots for BMIQ')
		qn_df = st.checkbox('Quantile Normalized Dataframe')
		mean_df = st.checkbox('Mean Normalized Dataframe')
		minmax_df = st.checkbox('Minmax Normalized Dataframe')


	if st.checkbox("Show Logged Data"):
		with tab_raw:
			plot.logged_plot()
	if st.button("Start Normalization", key="start"):
		with tab_norm:
			if b_mean_norm:
				plot.mean_normalized_plots(mean_df)
			if b_qn_norm:
				plot.quantile_normalized_plots(qn_df)
			if b_minmax_norm:
				plot.min_max_normalized_plots(minmax_df)
			if beta_plot:
				plot.beta_value_plots(beta_df, types_df)
			if b_bmiq_norm:
				plot.bmiq_plots(bmiq_df, boxplot_df)
# #TODO: DONE
# plot.beta_value_plots()
