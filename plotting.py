import matplotlib.pyplot as plt
import mpld3
import streamlit.components.v1 as components

import streamlit as st

from streamlit_extras.chart_container import chart_container
from multipledispatch import dispatch

import normalization as norm
import prepare_data as prep
from type import DataType


# can take dataframes with type column
def _density_plot(df, title):
	"""Makes a density-plot from the Dataframe, a title and set boundaries
	for the x axis."""
	plt.style.use('dark_background')
	df.plot.density(linewidth=1, figsize=(20, 10))
	plt.title(title)
	plt.legend(loc="upper right", ncols=2)
	plt.grid(True)


def _boxplot(df, title):
	"""Makes a boxplot from the Dataframe, a title and set boundaries
	for the x axis"""
	plt.style.use('dark_background')
	plt.boxplot(df)
	plt.title(title)
	plt.grid(True)


# can take dataframes with type column
def containerize_chart(df, name):
	"""Container makes usage of streamlit-extras possible."""
	with chart_container(df):
		st.write(name)
		st.pyplot(_density_plot(df, name))
		st.write("(Mean, STD): ", prep.output_measures(df))


def _switch(data_type, df_beta):
	"""Switch-case for differentiating between the different types of
	normalizations without encountering too much code duplication. Uses a
	custom enum for the type of Normalization provided (DataType)."""
	if data_type == DataType.BETA:
		return df_beta
	elif data_type == DataType.MEAN:
		df_beta_norm = norm.mean_normalization(df_beta)
		title = "Mean Normalization"
	elif data_type == DataType.MINMAX:
		df_beta_norm = norm.min_max_normalization(df_beta)
		title = "Min-Max Normalization"
	elif data_type == DataType.QN:
		df_beta_norm = norm.quantile_normalization('Median')
		title = "Quantile Normalization"
	return df_beta_norm, title


# return df to save as csv
@dispatch(DataType, object)
def default_plots(data_type, df_beta):
	"""Splits page into two columns to show methylated and unmethylated case
	side-by-side."""

	df_beta_norm, title = _switch(data_type, df_beta)
	col_left, col_right = st.columns(2)
	with col_left:
		containerize_chart(df_beta, "Raw Beta Values")
	with col_right:
		containerize_chart(df_beta_norm, title)
	return df_beta_norm


@dispatch(object, object, str)
def default_plots(df_t1, df_t2, title):
	"""Splits page into two columns to show methylated and unmethylated case
	side-by-side."""

	col_left, col_right = st.columns(2)
	with col_left:
		containerize_chart(df_t1, str(title + " Type 1 Probes"))
	with col_right:
		containerize_chart(df_t2, str(title + " Type 2 Probes"))


def quantile_normalized_plots():
	"""
	1. Set median to be the last column (acts as default reference).
	2. Build a selectbox to set your reference.
	3. Make normalization.
	TODO: Maybe build a reference dataframe?
	TODO: reload does not work.
	"""
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
		df_qn_meth = norm.quantile_normalization(df_log_meth, reference_meth)
		containerize_chart(df_qn_meth, "Quantile Normalized plot - methylated")
	with col_right:
		## QUANTILE NORMALIZED
		# last column as Median column
		df_log_unmeth['Median'] = df_log_unmeth.median(axis=1)
		reference_options = df_log_unmeth.columns.values.tolist()[2:]
		reference_unmeth = st.selectbox('Which reference do you want to use?',
		                                reference_options,
		                                (len(reference_options) - 1),
		                                key="reference")
		df_qn_unmeth = norm.quantile_normalization(df_log_unmeth,
		                                           reference_unmeth)
		containerize_chart(df_qn_meth, "Quantile Normalized plot - "
		                               "unmethylated")

	return df_qn_meth, df_qn_unmeth


@st.cache_data
def beta_value():
	df_meth, df_unmeth = prep.get_dataframe(True)
	df_beta = prep.beta_value(df_meth, df_unmeth, 100)
	return df_beta


def bmiq_plot(show_boxplot, df):
	"""Provides necessary dataframes and calls all corresponding methods."""
	df_beta = prep.add_probetypes(df)
	print('here')
	df_bmiq = norm.bmiq(df_beta)
	containerize_chart(df_bmiq, "BMIQ Values")
	if show_boxplot:
		_boxplots(df_bmiq)
	return df_bmiq


def _boxplots(df_bmiq):
	"""Builds evaliuation plots."""
	st.pyplot(_boxplot(df_bmiq, "BMIQ Normalized Boxplot"))
