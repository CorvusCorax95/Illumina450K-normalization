import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

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

df_meth = pd.read_csv('log_meth.csv', sep='\t', index_col=0)
df_unmeth = pd.read_csv('log_unmeth.csv', sep='\t', index_col=0)
df_beta = pd.read_csv('beta_values.csv', sep='\t', index_col=0)
df_t1 = pd.read_csv('beta_values_type1.csv', sep='\t', index_col=0)
df_t2 = pd.read_csv('beta_values_type2.csv', sep='\t', index_col=0)

df_bmiq = pd.read_csv('bmiq_norm.csv', sep='\t', index_col=0)
df_bmiq_t1 = pd.read_csv('bmiq_norm_t1.csv', sep='\t', index_col=0)
df_bmiq_t2 = pd.read_csv('bmiq_norm_t2.csv', sep='\t', index_col=0)

df_mean_m = pd.read_csv('mean_norm_m.csv', sep='\t', index_col=0)
df_mean_u = pd.read_csv('mean_norm_u.csv', sep='\t', index_col=0)
df_mean_m_t1 = pd.read_csv('mean_norm_m_t1.csv', sep='\t', index_col=0)
df_mean_u_t1 = pd.read_csv('mean_norm_u_t1.csv', sep='\t', index_col=0)
df_mean_m_t2 = pd.read_csv('mean_norm_m_t2.csv', sep='\t', index_col=0)
df_mean_u_t2 = pd.read_csv('mean_norm_u_t2.csv', sep='\t', index_col=0)

df_minmax_m = pd.read_csv('minmax_norm_m.csv', sep='\t', index_col=0)
df_minmax_u = pd.read_csv('minmax_norm_u.csv', sep='\t', index_col=0)
df_minmax_m_t1 = pd.read_csv('minmax_norm_m_t1.csv', sep='\t', index_col=0)
df_minmax_u_t1 = pd.read_csv('minmax_norm_u_t1.csv', sep='\t', index_col=0)
df_minmax_m_t2 = pd.read_csv('minmax_norm_m_t2.csv', sep='\t', index_col=0)
df_minmax_u_t2 = pd.read_csv('minmax_norm_u_t2.csv', sep='\t', index_col=0)

df_qn_u = pd.read_csv('qn_norm_u.csv', sep='\t', index_col=0)
df_qn_m = pd.read_csv('qn_norm_m.csv', sep='\t', index_col=0)
df_qn_u_t1 = pd.read_csv('qn_norm_u_t1.csv', sep='\t', index_col=0)
df_qn_m_t1 = pd.read_csv('qn_norm_m_t1.csv', sep='\t', index_col=0)
df_qn_u_t2 = pd.read_csv('qn_norm_u_t2.csv', sep='\t', index_col=0)
df_qn_m_t2 = pd.read_csv('qn_norm_m_t2.csv', sep='\t', index_col=0)

fig, ax = plt.subplots()

fig.suptitle("Value Ranges")

df_mean_u['Median'] = df_mean_u.median(axis=1)
df_mean_m['Median'] = df_mean_m.median(axis=1)
df_qn_u['Median'] = df_qn_u.median(axis=1)
df_qn_m['Median'] = df_qn_m.median(axis=1)
df_minmax_u['Median'] = df_minmax_u.median(axis=1)
df_minmax_m['Median'] = df_minmax_m.median(axis=1)
df_meth['Median'] = df_meth.median(axis=1)
df_unmeth['Median'] = df_unmeth.median(axis=1)

df_bmiq['Median'] = df_bmiq.median(axis=1)

plot = sns.boxplot(data=[
	df_meth['Median'],
	df_mean_m['Median'],
	df_minmax_m['Median'],
	df_qn_m['Median'],
	df_unmeth['Median'],
	df_mean_u['Median'],
	df_minmax_u['Median'],
	df_qn_u['Median'],
	df_bmiq['Median']])

xtick_loc = plot.get_xticks()
plot.set_xticks(ticks=xtick_loc, labels=[
	"Raw (m)",
	"mean (m)",
	"min-max (m)",
	"QN (m)",
	"raw (u)",
	"mean (u)",
	"min-max (u)",
	"QN (u)",
	"BMIQ"])

plt.show()
