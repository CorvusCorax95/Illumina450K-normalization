import matplotlib.pyplot as plt

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
	df_meth, df_unmeth = prep.get_values_as_dataframe_w_types()
	del df_meth["type"]
	del df_unmeth["type"]
	df_beta = prep.beta_value(df_meth, df_unmeth, 100)
	prep.df_to_h5(df_beta)
	st.write("Made data.h5.")
	# st.header("Beta Values")
	# st.write(df_beta)
	# st.write("Beta Values (Methylation Values)")
	# st.pyplot(_density_plot(df_beta, "Beta Values", -0.5, 1.5))
