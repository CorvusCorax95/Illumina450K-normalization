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

fig, ax = plt.subplots(2, 2)

fig.suptitle("Compare normalizations")

ax[0, 0].set_title("A) BMIQ Normalized Density", loc='left')
ax[0, 0].set_xlabel("beta-values")

ax[1, 0].set_title("B) Quantile Normalization Density", loc='left')
ax[1, 0].set_xlabel("beta-values")

ax[0, 1].set_title("C) Mean Normalization Density", loc='left')
ax[0, 1].set_xlabel("beta-values")

ax[1, 1].set_title("D) Min-Max Normalization Density", loc='left')
ax[1, 1].set_xlabel("beta-values")

sample_list = df_bmiq.columns.values.tolist()[1:]

del df_bmiq['type']
del df_mean_u['type']
del df_mean_m['type']
del df_minmax_u['type']
del df_minmax_m['type']

df_bmiq['Median'] = df_bmiq.median(axis=1)

df_beta_mean = pd.read_csv("mean_norm_beta.csv", sep='\t', index_col=0)
df_beta_minmax = pd.read_csv("minmax_norm_beta.csv", sep='\t', index_col=0)
df_beta_qn = pd.read_csv("qn_norm_beta.csv", sep='\t', index_col=0)
df_beta_bmiq = pd.read_csv("bmiq_norm.csv", sep='\t', index_col=0)

sample_list = df_beta_bmiq.columns.values.tolist()[1:]
for sample in sample_list:
	sns.kdeplot(df_beta_mean[sample], color=light_yellow, ax=ax[0, 1])
	sns.kdeplot(df_beta_minmax[sample], color=light_green, ax=ax[1, 1])
	sns.kdeplot(df_beta_bmiq[sample], color=light_purple, ax=ax[0, 0])

sample_list = df_beta_qn.columns.values.tolist()[1:]
for sample in sample_list:
	sns.kdeplot(df_beta_qn[sample], color=light_red, ax=ax[1, 0])

df_beta_mean['Median'] = df_beta_mean.median(axis=1)
df_beta_minmax['Median'] = df_beta_minmax.median(axis=1)
df_beta_qn['Median'] = df_beta_qn.median(axis=1)
df_beta_bmiq['Median'] = df_beta_qn.median(axis=1)

sns.kdeplot(df_beta_mean['Median'], color=dark_yellow, ax=ax[0, 1], legend=True)

sns.kdeplot(df_beta_minmax['Median'], color=dark_green, ax=ax[1, 1],
            legend=True)

sns.kdeplot(df_beta_qn['Median'], color=dark_red, ax=ax[1, 0], legend=True)

sns.kdeplot(df_bmiq['Median'], color=dark_purple, ax=ax[0, 0], legend=True)

plt.show()
