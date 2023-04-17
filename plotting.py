import matplotlib.pyplot as plt
import pandas as pd

import streamlit as st

import normalization as norm
import prepare_data as prep


# TODO: alle prep.function in die normalization verschieben

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
	# plt.legend([])
	plt.title(title)


def methylated_plots():
	df_log_meth, df_log_unmeth = prep.log_data()
	st.pyplot(_density_plot(df_log_meth, "Logged plot - methylated", 5, 17))

	# MEAN NORMALIZED
	df_meannorm_meth = norm.mean_normalization(df_log_meth)
	st.pyplot(
		_density_plot(df_meannorm_meth, "Mean Normalized - methylated", -4, 3))
	# MIN MAX NORMALIZED
	df_minmax_meth = norm.min_max_normalization(df_log_meth)
	st.pyplot(
		_density_plot(df_minmax_meth, "Min-Max Normalized - methylated", -0.2,
		              1.25))
	# QUANTILE NORMALIZED
	# last column as Median column
	df_log_meth['Median'] = df_log_meth.median(axis=1)
	reference_options = df_log_meth.columns.values.tolist()[2:]
	reference_meth = st.selectbox('Which reference do you want to use?',
	                              reference_options,
	                              (len(reference_options) - 1), key=0)
	df_qn_meth = norm.quantile_normaliziation(df_log_meth, reference_meth)
	st.pyplot(_density_plot(df_qn_meth,
	                        "Quantile Normalized plot (Median) - methylated", 5,
	                        17))


def unmethylated_plots():
	df_log_meth, df_log_unmeth = prep.log_data()
	st.pyplot(_density_plot(df_log_unmeth, "Logged plot - unmethylated", 5, 17))
	# MEAN NORMALIZED
	df_meannorm_unmeth = norm.mean_normalization(df_log_unmeth)
	st.pyplot(
		_density_plot(df_meannorm_unmeth, "Mean Normalized - unmethylated", -4,
		              3))
	## MIN MAX NORMALIZED
	df_minmax_unmeth = norm.min_max_normalization(df_log_unmeth)
	st.pyplot(
		_density_plot(df_minmax_unmeth, "Min-Max Normalized - unmethylated",
		              -0.2, 1.25))
	## QUANTILE NORMALIZED
	# last column as Median column
	df_log_unmeth['Median'] = df_log_unmeth.median(axis=1)
	reference_options = df_log_unmeth.columns.values.tolist()[2:]
	reference_unmeth = st.selectbox('Which reference do you want to use?',
	                                reference_options,
	                                (len(reference_options) - 1), key=1)
	df_qn_unmeth = norm.quantile_normaliziation(df_log_unmeth, reference_unmeth)
	st.pyplot(_density_plot(df_qn_unmeth,
	                        "Quantile Normalized plot (Median) - unmethylated",
	                        5, 17))


def beta_value_plots():
	st.header("Beta Values")
	df_beta = pd.read_csv('df_beta.csv', sep='\t')
	st.write("Beta Values (Methylation Values)")
	df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
	st.pyplot(_density_plot(df_beta_t1.loc[1:], "Beta Values Type 1",
	                        -0.5, 2.25))
	st.pyplot(_density_plot(df_beta_t2.loc[1:], "Beta Values type 2",
	                        -0.5, 2.25))
	df_bmiq = norm.bmiq()
	st.header("BMIQ Values")
	st.pyplot(_density_plot(df_bmiq, "BMIQ Values", -0.5, 2.25))
	st.write(df_beta_t2.loc[:, df_beta_t2.columns != 'type'])
	st.pyplot(_boxplot(df_beta_t2.loc[:, df_beta_t2.columns != 'type'],
	                   "Beta Value Boxplot"))
	st.write(df_bmiq)
	st.pyplot(_boxplot(df_bmiq, "BMIQ Normalized Boxplot"))


# st.pyplot(_boxplot(df_bmiq, "Title"))
# plot_random_dataframe(df_beta)


def plot_random_dataframe(df):
	from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, \
		DataReturnMode
	st.write(AgGrid(df))
