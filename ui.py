import shutil

import pandas as pd
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
		Download everything at once: Get an additional button that downloads 
		a zip-File with all calculated dataframes as .csv
		Normalize non-logged intensities: This means that the initial 
		normalization is performed on non-logged 
		data. After the normalization the log is calulated before 
		calculating the beta values.
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
			b_log = st.checkbox("Normalize logarithmized intensities.")
			st.write("What normalizations do you want to see?")
			b_mean = st.checkbox('Mean Normalization')
			b_qn = st.checkbox('Quantile Normalization')
			b_minmax = st.checkbox('Min-Max Normalization')
			b_qn_bmiq = st.checkbox('Beta-Mixture-Quantile Normalization')
			st.info('Additional Info.')
			show_types = st.checkbox('Additional view split by types')
			show_boxplot = st.checkbox('Additional boxplots for '
			                           'normalizations')
			show_dens = st.checkbox('Additional density graphs')
			submitted = st.form_submit_button("Run Normalization")

	if b_log:
		df_meth, df_unmeth = prep.log_data()
	else:
		df_meth, df_unmeth = prep.get_dataframe(True)


	df_sample_to_numbers = pd.DataFrame(df_meth.columns.values.tolist(),
	                                    columns=["samples"])
	multiselect = df_meth.columns.values.tolist()
	options = st.multiselect('Select Samples.',
	                         multiselect, default=[
			"RB_E_001", "RB_E_002", "RB_E_003"], )
	col_left, col_right = st.columns(2)
	with col_left:
		st.write(df_meth[options])
	with col_right:
		st.write(df_unmeth[options])
	df_meth = df_meth[options]
	df_unmeth = df_unmeth[options]
	if submitted:
		make_plots(df_meth, df_unmeth, b_mean, b_minmax, b_qn, b_qn_bmiq,
		           show_types, show_boxplot, show_dens, df_sample_to_numbers)
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


def make_plots(df_meth, df_unmeth, b_mean, b_minmax, b_qn, b_bmiq,
               show_types, show_boxplot, show_dens, df_sample_to_numbers):
	"""Cares about all the plots and Dataframe-Views"""

	with st.spinner("Waiting for data..."):
		# ORIGINAL DATA
		plot.meth_plots(df_meth, df_unmeth, "Raw Light Intensities")
		user.convert_df(df_meth, 'download/raw-meth-values.csv')
		user.convert_df(df_unmeth, 'download/raw-unmeth-values.csv')
		df_beta = prep.beta_value(df_meth, df_unmeth, 100)
		plot.default_plots(DataType.BETA, df_beta)
		user.convert_df(df_beta, 'download/raw-beta-values.csv')
		if show_boxplot:
			plot.boxplots(df_beta, "Raw Beta Values")
		if show_types:
			df_meth = prep.add_probetypes(df_meth)
			df_unmeth = prep.add_probetypes(df_unmeth)
			df_meth_t1, df_meth_t2 = prep.split_types(df_meth)
			df_unmeth_t1, df_unmeth_t2 = prep.split_types(df_unmeth)
			plot.default_plots(df_meth_t1, df_meth_t2, "Methylated Values")
			plot.default_plots(df_unmeth_t1, df_unmeth_t2,
			                   "Unmethylated Values")
		if show_dens:
			st.pyplot(plot.density(df_beta, "Density of Beta Values"))

	with st.spinner("Wait for Normalization..."):
		if b_mean:
			st.subheader("Mean Normalization")
			st.subheader("_Methylated File_")
			df_meth_mean = plot.default_plots(DataType.MEAN, df_meth)
			st.subheader("_Unmethylated File_")
			df_unmeth_mean = plot.default_plots(DataType.MEAN, df_unmeth)
			user.convert_df(df_meth_mean, 'download/meannorm_meth.csv')
			user.convert_df(df_unmeth_mean, 'download/meannorm_unmeth.csv')
			df_mean_beta = plot.beta_value(df_meth_mean, df_unmeth_mean)
			plot.containerize_chart(df_mean_beta, "Mean Normalization (Beta "
			                                      "Values)")
			user.convert_df(df_mean_beta, 'download/mean-beta-values.csv')
			if show_boxplot:
				plot.boxplots(df_mean_beta, "Mean Normalization")
			if show_types:
				df_mean_beta = prep.add_probetypes(df_mean_beta)
				df_mean_t1, df_mean_t2 = prep.split_types(df_mean_beta)
				plot.default_plots(df_mean_t1, df_mean_t2,
				                   "Mean Normalized Beta Values")
			if show_dens:
				st.pyplot(plot.density(df_mean_beta, "Density of Mean "
				                                     "Normalized Beta Values"))
	with st.spinner("Wait for Normalization..."):
		if b_minmax:
			st.subheader("Minmax Normalization")
			st.subheader("_Methylated File_")
			df_meth_minmax = plot.default_plots(DataType.MINMAX, df_meth)
			st.subheader("_Unmethylated File_")
			df_unmeth_minmax = plot.default_plots(DataType.MINMAX, df_unmeth)
			user.convert_df(df_meth_minmax, 'download/minmaxnorm_meth.csv')
			user.convert_df(df_unmeth_minmax, 'download/minmaxnorm_unmeth.csv')
			df_minmax_beta = plot.beta_value(df_meth_minmax, df_unmeth_minmax)
			plot.containerize_chart(df_minmax_beta, "Minmax Normalization ("
			                                        "Beta "
			                                        "Values)")
			user.convert_df(df_minmax_beta, 'download/minmax-beta-values.csv')
			if show_boxplot:
				plot.boxplots(df_minmax_beta, "Minmax Normalization")
			if show_types:
				df_minmax_beta = prep.add_probetypes(df_minmax_beta)
				df_minmax_t1, df_minmax_t2 = prep.split_types(df_minmax_beta)
				plot.default_plots(df_minmax_t1, df_minmax_t2, "Minmax "
				                                               "Normalized Beta Values")
			if show_dens:
				st.pyplot(plot.density(df_minmax_beta, "Density of Minmax "
				                                       "Normalized Beta Values"))
	with st.spinner("Wait for Normalization..."):
		if b_qn:
			st.subheader("Quantile Normalization")
			st.subheader("_Methylated File_")
			df_meth_qn = plot.default_plots(DataType.QN, df_meth)
			st.subheader("_Unmethylated File_")
			df_unmeth_qn = plot.default_plots(DataType.QN, df_unmeth)
			user.convert_df(df_meth_qn, 'download/qnorm_meth.csv')
			user.convert_df(df_unmeth_qn, 'download/qnorm_unmeth.csv')
			df_qn_beta = plot.beta_value(df_meth_qn, df_unmeth_qn)
			plot.containerize_chart(df_qn_beta,
			                        "Quantile Normalization (Beta Values)")
			user.convert_df(df_qn_beta, 'download/qn-beta-values.csv')
			if show_boxplot:
				plot.boxplots(df_qn_beta, "Quantile Normalization")
			if show_types:
				df_qn_w_types = prep.add_probetypes(df_qn_beta)
				df_qn_t1, df_qn_t2 = prep.split_types(df_qn_w_types)
				plot.default_plots(df_qn_t1, df_qn_t2, "Quantile Normalized "
				                                       "Beta Values")
			if show_dens:
				st.pyplot(plot.density(df_qn_beta, "Density of Quantile "
				                                   "Normalized Beta Values"))
	with st.spinner("Wait for Normalization..."):
		if b_bmiq:
			st.subheader("BMIQ Normalization")
			st.info("The BMIQ Normalization as presented in Teschendorff et "
			        "al. (2013) only proposed a normalization for type 2 "
			        "probe values. To see the normalization better, tick the "
			        "box 'Additional view: split by types'")
			df_bmiq = plot.bmiq_plot(df_meth, df_sample_to_numbers)
			user.convert_df(df_bmiq, 'download/bmiq.csv')
			if show_dens:
				st.pyplot(plot.density(df_bmiq, "Density of BMIQ Normalized "
				                                "Beta Values"))
			if show_boxplot:
				plot.boxplots(df_bmiq, "BMIQ Normalization")
			if show_types:
				df_bmiq_w_types = prep.add_probetypes(df_bmiq)
				df_bmiq_t1, df_bmiq_t2 = prep.split_types(df_bmiq_w_types)
				plot.default_plots(df_bmiq_t1, df_bmiq_t2, "BMIQ "
				                                           "Normalized Beta Values")
