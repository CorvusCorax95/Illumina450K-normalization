import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

import normalization as norm

# WORKS

sns.set()
# ------------------------------------------------------------------#
light_yellow = (1, 0.949, 0.639, 1)  # Type 1 Methylated
dark_yellow = (0.91, 0.635, 0, 1)  # Type 1 Methylated
hex_yellow = "#e8a200"

light_blue = (0.557, 0.6, 1, 1)  # Type 1 Unmethylated
dark_blue = (0, 0.314, 0.541, 1)  # Type 1 Unmethylated
hex_blue = "#00508a"

light_green = (0.62, 0.839, 0.502, 1)  # Type 2 Methylated
dark_green = (0.012, 0.369, 0.067, 1)  # Type 2 Methylated
hex_green = "#035e00"

light_red = (1, 0.682, 0.643, 1)  # Type 2 Unmethylated
dark_red = (0.541, 0.043, 0, 1)  # Type 2 Unmethylated
hex_red = "#8a0b00"

light_purple = (0.678, 0.596, 0.839, 1)  # Type 2 Unmethylated
dark_purple = (0.212, 0.047, 0.541, 1)  # Type 2 Unmethylated
hex_purple = "#360c8a"

df_beta = pd.read_csv('beta_values.csv', sep='\t', index_col=0)

df_bmiq = pd.read_csv('bmiq_norm.csv', sep='\t', index_col=0)
df_bmiq_t1 = pd.read_csv('bmiq_norm_t1.csv', sep='\t', index_col=0)
df_bmiq_t2 = pd.read_csv('bmiq_norm_t2.csv', sep='\t', index_col=0)

df_beta_t1 = pd.read_csv('beta_values_type1.csv', sep='\t', index_col=0)
df_beta_t2 = pd.read_csv('beta_values_type2.csv', sep='\t', index_col=0)

df_mean_beta_t1 = norm.mean_normalization(df_beta_t1)
df_mean_beta_t2 = norm.mean_normalization(df_beta_t2)

df_minmax_beta_t1 = norm.min_max_normalization(df_beta_t1)
df_minmax_beta_t2 = norm.min_max_normalization(df_beta_t2)

df_beta_t1['Median'] = df_beta_t1.median(axis=1)
df_beta_t2['Median'] = df_beta_t2.median(axis=1)

df_qn_beta_t1 = norm.quantile_normaliziation(df_beta_t1, 'Median')
df_qn_beta_t2 = norm.quantile_normaliziation(df_beta_t2, 'Median')

fig, ax = plt.subplots()

fig.suptitle("Value Ranges of Normalizations")


df_beta_t1['Median'] = df_beta_t1.median(axis=1)
df_beta_t2['Median'] = df_beta_t2.median(axis=1)
df_bmiq['Median'] = df_bmiq.median(axis=1)
df_mean_beta_t1['Median'] = df_mean_beta_t1.median(axis=1)
df_mean_beta_t2['Median'] = df_mean_beta_t2.median(axis=1)
df_minmax_beta_t1['Median'] = df_minmax_beta_t1.median(axis=1)
df_minmax_beta_t2['Median'] = df_minmax_beta_t2.median(axis=1)
df_qn_beta_t1['Median'] = df_qn_beta_t1.median(axis=1)
df_qn_beta_t2['Median'] = df_qn_beta_t2.median(axis=1)


plot = sns.boxplot(data=[
	df_beta_t1['Median'],
	df_beta_t2['Median'],
	df_mean_beta_t1['Median'],
	df_mean_beta_t2['Median'],
	df_minmax_beta_t1['Median'],
	df_minmax_beta_t2['Median'],
	df_qn_beta_t1['Median'],
	df_qn_beta_t2['Median'],
	df_bmiq['Median']])

xtick_loc = plot.get_xticks()
plot.set_xticks(ticks=xtick_loc, labels=[
	"Beta T1",
	"Beta T2",
	"Mean T1",
	"Mean T2",
	"Minmax T1",
	"Minmax T2",
	"QN T1",
	"QN T2",
	"BMIQ"])

plt.show()
