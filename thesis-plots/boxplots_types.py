import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import prepare_data as prep

sns.set(font_scale=1.6)

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
df_bmiq_types = prep.add_probetypes(df_bmiq)

sns.set_style('whitegrid')

y = 'RB_E_007'
x = 'type'
sns.set_palette('Spectral')
fig, ax = plt.subplots(2, 3)


ax[0][0].set_title('Beta Values')

ax[0][1].set_title('Mean Normalization')

ax[1][0].set_title('Minmax Normalization')

ax[1][1].set_title('Quantile Normalization')

ax[0][2].set_title('BMIQ')

ax[1][2].set_axis_off()

sns.boxplot(data=df_beta_types, x=y, y=x, ax=ax[0][0])

sns.boxplot(data=df_mean_types, x=y, y=x, ax=ax[0][1])

sns.boxplot(data=df_minmax_types, x=y, y=x, ax=ax[1][0])

sns.boxplot(data=df_qn_types, x=y, y=x, ax=ax[1][1])

sns.boxplot(data=df_bmiq_types, x=y, y=x, ax=ax[0][2])

fig2, ax = plt.subplots(1)

ax.set_title('Beta values')

d = {'Raw beta Values': df_beta_types[y],
     'Mean normalized'    : df_mean_types[y],
     'Minmax Normalized'  : df_minmax_types[y],
     'Quantile Normalized': df_qn_types[y],
     'BMIQ Normalized': df_bmiq_types[
	     y]}
df = pd.DataFrame(data=d)

sns.boxplot(data=df)

plt.show()
