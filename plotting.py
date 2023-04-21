import matplotlib.pyplot as plt
import pandas as pd

import streamlit as st

import normalization as norm
import prepare_data as prep


def _density_plot(df, title, x1, x2):
	"""Makes a density-plot from the Dataframe, a title and set boundaries
	for the x axis"""
	plt.style.use('dark_background')
	df.plot.density(linewidth=1, figsize=(20, 10), xlim=(x1, x2))
	plt.legend([])
	plt.title(title)


def _boxplot(df, title):
	"""Makes a density-plot from the Dataframe, a title and set boundaries
	for the x axis"""
	plt.style.use('dark_background')
	plt.boxplot(df)
	plt.title(title)


def logged_plot():
	col_left, col_right = st.columns(2)
	df_log_meth, df_log_unmeth = prep.log_data()
	with col_left:
		st.pyplot(_density_plot(df_log_meth, "Logged plot - methylated", 5, 17))
	with col_right:
		st.pyplot(
			_density_plot(df_log_unmeth, "Logged plot - unmethylated", 5, 17))


def mean_normalized_plots(show_df):
	df_log_meth, df_log_unmeth = prep.log_data()
	# MEAN NORMALIZED
	col_left, col_right = st.columns(2)
	# MIN MAX NORMALIZED
	with col_left:
		st.subheader("Methylated Values")
		df_meannorm_meth = norm.mean_normalization(df_log_meth)
		st.pyplot(
			_density_plot(df_meannorm_meth, "Mean Normalized - methylated", -4,
			              3))
		if show_df:
			st.write(df_meannorm_meth)
	with col_right:
		st.subheader("Unmethylated Values")
		# MEAN NORMALIZED
		df_meannorm_unmeth = norm.mean_normalization(df_log_unmeth)
		st.pyplot(
			_density_plot(df_meannorm_unmeth, "Mean Normalized - unmethylated",
			              -4,
			              3))
		if show_df:
			st.write(df_meannorm_unmeth)
	return df_meannorm_meth, df_meannorm_unmeth


def min_max_normalized_plots(show_df):
	df_log_meth, df_log_unmeth = prep.log_data()
	col_left, col_right = st.columns(2)
	# MIN MAX NORMALIZED
	with col_left:
		st.subheader("Methylated Values")
		df_minmax_meth = norm.min_max_normalization(df_log_meth)
		st.pyplot(
			_density_plot(df_minmax_meth, "Min-Max Normalized - methylated",
			              -0.2,
			              1.25))
		if show_df:
			st.write(df_minmax_meth)
	## MIN MAX NORMALIZED
	with col_right:
		st.subheader("Unmethylated Values")
		df_minmax_unmeth = norm.min_max_normalization(df_log_unmeth)
		st.pyplot(
			_density_plot(df_minmax_unmeth, "Min-Max Normalized - unmethylated",
			              -0.2, 1.25))
		if show_df:
			st.write(df_minmax_unmeth)


def quantile_normalized_plots(show_df):
	df_log_meth, df_log_unmeth = prep.log_data()
	# QUANTILE NORMALIZED
	col_left, col_right = st.columns(2)
	with col_left:
		st.subheader("Methylated Values")
		# last column as Median column
		df_log_meth['Median'] = df_log_meth.median(axis=1)
		reference_options = df_log_meth.columns.values.tolist()[2:]
		reference_meth = st.selectbox('Which reference do you want to use?',
		                              reference_options,
		                              (len(reference_options) - 1), key=0)
		df_qn_meth = norm.quantile_normaliziation(df_log_meth, reference_meth)
		if show_df:
			st.write(df_qn_meth)
		st.pyplot(_density_plot(df_qn_meth,
		                        "Quantile Normalized plot (Median) - methylated",
		                        5,
		                        17))
	with col_right:
		st.subheader("Unmethylated Values")

		## QUANTILE NORMALIZED
		# last column as Median column
		df_log_unmeth['Median'] = df_log_unmeth.median(axis=1)
		reference_options = df_log_unmeth.columns.values.tolist()[2:]
		reference_unmeth = st.selectbox('Which reference do you want to use?',
		                                reference_options,
		                                (len(reference_options) - 1), key=1)
		df_qn_unmeth = norm.quantile_normaliziation(df_log_unmeth,
		                                            reference_unmeth)
		if show_df:
			st.write(df_qn_unmeth)
		st.pyplot(_density_plot(df_qn_unmeth,
		                        "Quantile Normalized plot (Median) - unmethylated",
		                        5, 17))


@st.cache_data
def everything_beta(beta, types, bmiq, boxplot, bmiq_norm, beta_plot):
	df_log_meth, df_log_unmeth = prep.log_data()
	df_beta = prep.beta_value(df_log_meth, df_log_unmeth, 100)
	st.write("Beta Values (Methylation Values)")
	df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
	print(df_beta_t1)
	if beta_plot: #fine
		_beta_value_plots(beta, types, df_beta, df_beta_t1, df_beta_t2)
	if bmiq_norm:
		_bmiq_plots(bmiq, boxplot, df_beta_t2)


@st.cache_data
def _beta_value_plots(beta, types, df_beta, df_beta_t1, df_beta_t2):
	col_left, col_right = st.columns(2)
	if beta:
		st.write(df_beta)
	with col_left:
		if types:
			st.write(df_beta_t1)
		st.pyplot(_density_plot(df_beta_t1, "Beta Values Type 1",
		                        -0.1, 0.2))
	with col_right:
		if types:
			st.write(df_beta_t2)
		st.pyplot(_density_plot(df_beta_t2, "Beta Values type 2",
		                        -0.1, 0.2))


@st.cache_data
def _bmiq_plots(bmiq, boxplot, df_beta_t2):
	df_bmiq = norm.bmiq()
	if bmiq:
		st.write(df_bmiq)
	st.pyplot(_density_plot(df_bmiq, "BMIQ Values", -0.5, 2.25))
	if boxplot:
		bmiq_boxplots(df_beta_t2, df_bmiq)


def bmiq_boxplots(df_beta_t2, df_bmiq):
	st.pyplot(_boxplot(df_beta_t2.loc[:, df_beta_t2.columns != 'type'],
	                   "Beta Value Boxplot"))

	st.pyplot(_boxplot(df_bmiq, "BMIQ Normalized Boxplot"))


def plot_random_dataframe(df):
	from st_aggrid import AgGrid
	st.write(AgGrid(df))
