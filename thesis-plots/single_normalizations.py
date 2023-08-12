import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import prepare_data as prep

sns.set(font_scale=1.7)
# ------------------------------------------------------------------#
light_yellow = (1, 0.949, 0.639, 1)  # Type 1 Methylated
dark_yellow = (0.91, 0.635, 0, 1)  # Type 1 Methylated
hex_yellow = "#e8a200"

light_red = (1, 0.682, 0.643, 1)  # Type 2 Unmethylated
dark_red = (0.541, 0.043, 0, 1)  # Type 2 Unmethylated
hex_red = "#8a0b00"

light_purple = (0.678, 0.596, 0.839, 1)  # beta
dark_purple = (0.212, 0.047, 0.541, 1)  # beta
hex_purple = "#360c8a"

df_beta_types = pd.read_csv(
	'compressed_download(7)/raw-beta-values_w_types.csv', index_col=0)
df_mean = pd.read_csv(
	'compressed_download(7)/meannorm.csv', index_col=0)
df_mean_types = prep.add_probetypes(df_mean)

df_minmax = pd.read_csv(
	'compressed_download(7)/minmaxnorm.csv', index_col=0)
df_minmax_types = prep.add_probetypes(df_minmax)

df_qn = pd.read_csv(
	'compressed_download(7)/qnorm.csv', index_col=0)
df_qn_types = prep.add_probetypes(df_qn)

df_bmiq = pd.read_csv(
	'compressed_download(7)/bmiq.csv', index_col=0)
df_bmiq_types = prep.add_probetypes(df_qn)

fig, ax = plt.subplots(2, sharey=True)

df = df_bmiq_types
df_1, df_2 = prep.split_types(df)
df_t1 = df_1.copy()
df_t2 = df_2.copy()

ax[0].set_title("Combined Quantile Normalization")
ax[1].set_title("Split Quantile Normalization")

for sample in df.columns.values.tolist():
	if sample == 'type':
		continue
	sns.kdeplot(data=df[sample], ax=ax[0], color=light_purple)
	sns.kdeplot(data=df_t1[sample], ax=ax[1], color=light_yellow)
	sns.kdeplot(data=df_t2[sample], ax=ax[1], color=light_red)

del df['type']

df_med = df.copy()

df_med['Median'] = df.median(axis=1)
df_t1['Median'] = df_t1.median(axis=1)
df_t2['Median'] = df_t2.median(axis=1)

sns.kdeplot(data=df_med['Median'], ax=ax[0], color=dark_purple).set(xlabel="beta-value")
sns.kdeplot(data=df_t1['Median'], ax=ax[1], color=dark_yellow).set(xlabel="beta-value")
sns.kdeplot(data=df_t2['Median'], ax=ax[1], color=dark_red).set(
	xlabel="beta-value")



plt.show()
