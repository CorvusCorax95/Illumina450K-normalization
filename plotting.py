import matplotlib.pyplot as plt
import streamlit as st
import seaborn as sns

from streamlit_extras.chart_container import chart_container
from multipledispatch import dispatch

import normalization as norm
import prepare_data as prep
from type import DataType


# can take dataframes with type column
def _hist_plot(df):
	"""Makes a density-plot from the Dataframe, a title and set boundaries
	for the x axis."""
	plt.style.use('dark_background')
	for sample in df.columns.values.tolist():
		sns.histplot(df[sample], kde=False, bins=20,
		             label=sample, element='poly', fill=False)
	plt.legend(loc="upper right", ncols=2)
	plt.grid(True, alpha=0.35)


def density(df, title):
	sns.kdeplot(df)
	plt.title(title)
	plt.legend(loc="upper right", ncols=2)
	plt.grid(True, alpha=0.35)

# can take dataframes with type column
def containerize_chart(df, name):
	"""Container makes usage of streamlit-extras possible."""
	if "type" in df.columns.values.tolist():
		df.drop(['type'], axis=1)
	if "Median" in df.columns.values.tolist():
		df.drop(['Median'], axis=1)
	with chart_container(df):
		st.write(name)
		st.pyplot(_hist_plot(df))
		st.write("(Mean, Variance, STD): ", prep.output_measures(df))



def _switch(data_type, df):
	"""Switch-case for differentiating between the different types of
	normalizations without encountering too much code duplication. Uses a
	custom enum for the type of Normalization provided (DataType)."""
	if 'type' in df.columns.values.tolist():
		df = df.drop(['type'], axis=1)
	if data_type == DataType.BETA:
		return df, "Beta Values"
	elif data_type == DataType.MEAN:
		df_norm = norm.mean_normalization(df)
		title = "Mean Normalization"
	elif data_type == DataType.MINMAX:
		df_norm = norm.min_max_normalization(df)
		title = "Min-Max Normalization"
	elif data_type == DataType.QN:
		sample_list = df.columns.values.tolist()
		df_norm = norm.quantile_normalization(df, 'Median')
		title = "Quantile Normalization"
	return df_norm, title


# return df to save as csv
@dispatch(DataType, object)
def default_plots(data_type, df):
	"""Splits page into two columns to show methylated and unmethylated case
	side-by-side."""
	df_norm, title = _switch(data_type, df)
	containerize_chart(df_norm, title)
	return df_norm


@dispatch(object, object, str)
def default_plots(df_t1, df_t2, title):
	"""Splits page into two columns to show methylated and unmethylated case
	side-by-side."""
	col_left, col_right = st.columns(2)
	with col_left:
		containerize_chart(df_t1, str(title + " Type 1 Probes"))
	with col_right:
		containerize_chart(df_t2, str(title + " Type 2 Probes"))

def meth_plots(df_meth, df_unmeth, title):
	"""Splits page into two columns to show methylated and unmethylated case
	side-by-side."""
	col_left, col_right = st.columns(2)
	with col_left:
		containerize_chart(df_meth, str(title + "(Methylated)"))
	with col_right:
		containerize_chart(df_unmeth, str(title + "(Unmethylated)"))


@st.cache_data
def beta_value():
	df_meth, df_unmeth = prep.get_dataframe(True)
	df_beta = prep.beta_value(df_meth, df_unmeth, 100)
	return df_beta

@dispatch(object, object)
@st.cache_data
def beta_value(df_meth, df_unmeth):
	df_beta = prep.beta_value(df_meth, df_unmeth, 100)
	return df_beta


def bmiq_plot(df_meth, df_sample_to_numbers):
	"""Provides necessary dataframes and calls all corresponding methods."""
	df_bmiq = norm.bmiq(df_meth, df_sample_to_numbers)
	containerize_chart(df_bmiq, "BMIQ Values")
	return df_bmiq


def boxplots(df):
	"""Builds evaluation plots."""
	columns = df.columns.values.tolist()
	with chart_container(df[columns]):
		st.write("Boxplot")
		st.pyplot(
			_boxplot_df(df[columns]))


def _boxplot_df(df):
	"""Makes a boxplot from the Dataframe, a title and set boundaries
	for the x axis"""
	boxprops = dict(color='lightblue', linewidth=5)
	flierprops = dict(marker='o', markerfacecolor='firebrick')
	medianprops = dict(linewidth=5, color='yellow')
	meanprops = dict(marker='D', markeredgecolor='green',
	                 markerfacecolor='green',
	                 linewidth=5)
	labels = df.columns.values.tolist()
	plt.boxplot(df, labels=labels, showfliers=True, showmeans=True,
	            showcaps=True, showbox=True, flierprops=flierprops,
	            boxprops=boxprops, meanprops=meanprops, medianprops=medianprops)
	plt.style.use('dark_background')
	plt.grid(True, alpha=0.35)
