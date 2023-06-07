import shutil

import streamlit as st

import plotting as plot
import user
from streamlit_extras.dataframe_explorer import dataframe_explorer
from streamlit_extras.chart_container import chart_container

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


def make_plots():
	"""Cares about all the plots and Dataframe-Views"""
	with st.sidebar:
		b_download = st.checkbox("Do you want the possibility to download "
		                         "everything at once? Clicking the provided "
		                         "button will reload the site entirely.")
		st.write("What normalizations do you want to see?")
		b_mean_norm = st.checkbox('Mean Normalization')
		b_qn_norm = st.checkbox('Quantile Normalization')
		b_minmax_norm = st.checkbox('Min-Max Normalization')
		b_bmiq_norm = st.checkbox('Beta-Mixture-Quantile Normalization')
		st.info('Decide here if you want to see the Additional Info.')
		show_boxplot_df = st.checkbox('Boxplots for BMIQ')

	if st.button("Start Normalization", key="start"):
		with st.spinner("Wait for Normalization..."):
			# ORIGINAL DATA
			df_meth, df_unmeth = plot.default_plots(DataType.LOG)
			user.convert_df(df_meth, 'download/log-meth.csv')
			user.convert_df(df_unmeth, 'download/log-unmeth.csv')

			# NORMALIZATION
			if b_mean_norm:
				df_meth, df_unmeth = plot.default_plots(DataType.MEAN)
				user.convert_df(df_meth, 'download/meannorm-meth.csv')
				user.convert_df(df_unmeth, 'download/meannorm-unmeth.csv')
			if b_minmax_norm:
				df_meth, df_unmeth = plot.default_plots(DataType.MINMAX)
				user.convert_df(df_meth, 'download/minmaxnorm-meth.csv')
				user.convert_df(df_unmeth, 'download/minmaxnorm-unmeth.csv')
			if b_qn_norm:
				df_meth, df_unmeth = plot.quantile_normalized_plots()
				user.convert_df(df_meth, 'download/qnnorm-meth.csv')
				user.convert_df(df_unmeth, 'download/qnnorm-unmeth.csv')
			if b_bmiq_norm:
				df = plot.everything_beta(show_boxplot_df)
				user.convert_df(df, 'download/bmiq.csv')

		st.success("Normalization completed!")
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
