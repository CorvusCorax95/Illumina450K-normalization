import matplotlib.pyplot as plt

import streamlit as st

from streamlit_extras.chart_container import chart_container

import normalization as norm
import prepare_data as prep
from type import DataType


def _density_plot(df, title):
	"""Makes a density-plot from the Dataframe, a title and set boundaries
	for the x axis"""
	plt.style.use('dark_background')
	df.plot.density(linewidth=1, figsize=(20, 10))
	plt.legend([])
	plt.title(title)


def _boxplot(df, title):
	"""Makes a density-plot from the Dataframe, a title and set boundaries
	for the x axis"""
	plt.style.use('dark_background')
	plt.boxplot(df)
	plt.title(title)


def _containerize_chart(df, name):
	with chart_container(df):
		st.write(name)
		st.pyplot(_density_plot(df, name))


def _switch(data_type, df_log_meth, df_log_unmeth):
	if data_type == DataType.LOG:
		return df_log_meth, df_log_unmeth, "Logged Plots"
	elif data_type == DataType.MEAN:
		df_norm_meth = norm.mean_normalization(df_log_meth)
		df_norm_unmeth = norm.mean_normalization(df_log_unmeth)
		title = "Mean Normalization"
	elif data_type == DataType.MINMAX:
		df_norm_meth = norm.min_max_normalization(df_log_meth)
		df_norm_unmeth = norm.min_max_normalization(df_log_unmeth)
		title = "Min-Max Normalization"
	elif data_type == DataType.QN:
		df_norm_meth = norm.quantile_normaliziation(df_log_meth)
		df_norm_unmeth = norm.quantile_normaliziation(df_log_unmeth)
		title = "Quantile Normalization"
	return df_norm_meth, df_norm_unmeth, title


def default_plots(data_type):
	df_meth, df_unmeth = prep.log_data()
	df_norm_meth, df_norm_unmeth, title = _switch(data_type, df_meth, df_unmeth)

	col_left, col_right = st.columns(2)
	with col_left:
		_containerize_chart(df_norm_meth, str(title + " - methylated"))
	with col_right:
		_containerize_chart(df_norm_unmeth, str(title + " - unmethylated"))
	return df_norm_meth, df_norm_unmeth


def quantile_normalized_plots():
	st.header("Quantile Normalization")
	df_log_meth, df_log_unmeth = prep.log_data()
	# QUANTILE NORMALIZED
	col_left, col_right = st.columns(2)
	with col_left:
		df_log_meth['Median'] = df_log_meth.median(axis=1)
		reference_options = df_log_meth.columns.values.tolist()[2:]
		reference_meth = st.selectbox('Which reference do you want to use?',
		                              reference_options,
		                              (len(reference_options) - 1), key=0)
		df_qn_meth = norm.quantile_normaliziation(df_log_meth, reference_meth)
		_containerize_chart(df_qn_meth, "Quantile Normalized plot - methylated")
	with col_right:
		## QUANTILE NORMALIZED
		# last column as Median column
		df_log_unmeth['Median'] = df_log_unmeth.median(axis=1)
		reference_options = df_log_unmeth.columns.values.tolist()[2:]
		reference_unmeth = st.selectbox('Which reference do you want to use?',
		                                reference_options,
		                                (len(reference_options) - 1),
		                                key="reference")
		df_qn_unmeth = norm.quantile_normaliziation(df_log_unmeth,
		                                            reference_unmeth)
		_containerize_chart(df_qn_meth, "Quantile Normalized plot - "
		                                "unmethylated")

	return df_qn_meth, df_qn_unmeth


def everything_beta(boxplot):
	df_log_meth, df_log_unmeth = prep.log_data()
	df_beta = prep.beta_value(df_log_meth, df_log_unmeth, 100)
	st.write("Beta Values (Methylation Values)")
	df_beta_t1, df_beta_t2 = prep.split_types(df_beta)

	_beta_value_plots(df_beta, df_beta_t1, df_beta_t2, boxplot)
	st.write("BMIQ Normalization")

	return df_beta


def _beta_value_plots(df_beta, df_beta_t1, df_beta_t2, boxplot):
	_containerize_chart(df_beta, "Beta Values")

	col_left, col_right = st.columns(2)
	with col_left:
		_containerize_chart(df_beta_t1, "Beta Values Type 1")
	with col_right:
		_containerize_chart(df_beta_t2, "Beta Values Type 2")

	df_bmiq = norm.bmiq()
	_containerize_chart(df_bmiq, "BMIQ Values")
	if boxplot:
		_boxplots(df_beta_t2, df_bmiq)


def _boxplots(df_beta_t2, df_bmiq):
	st.pyplot(_boxplot(df_beta_t2.loc[:, df_beta_t2.columns != 'type'],
	                   "Beta Value Boxplot"))
	st.pyplot(_boxplot(df_bmiq, "BMIQ Normalized Boxplot"))


def plot_random_dataframe(df):
	from st_aggrid import AgGrid
	st.write(AgGrid(df))
