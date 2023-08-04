import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

import normalization as norm

sns.set(font_scale=2)
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

light_purple = (0.678, 0.596, 0.839, 1)  # beta
dark_purple = (0.212, 0.047, 0.541, 1)  # beta
hex_purple = "#360c8a"


df_beta = pd.read_csv('beta_values.csv', sep='\t', index_col=0)

df_beta_t1 = pd.read_csv('beta_values_type1.csv', sep='\t', index_col=0)
df_beta_t2 = pd.read_csv('beta_values_type2.csv', sep='\t', index_col=0)

df_mean_beta_t1 = norm.mean_normalization(df_beta_t1)
df_mean_beta_t2 = norm.mean_normalization(df_beta_t2)


fig, ax = plt.subplots(2, sharex=True)

fig.suptitle("Mean Normalizations")

ax[0].set_title("A) Type 1 beta-value densities", loc='left')
ax[0].set_xlabel("beta-values")

ax[1].set_title("B) Type 2 beta-value densities", loc='left')
ax[1].set_xlabel("beta-values")

sample_list = df_mean_beta_t1.columns.values.tolist()[1:]

for sample in sample_list:
	sns.kdeplot(df_mean_beta[sample], color=light_blue, ax=ax[1])
	sns.kdeplot(df_mean_beta_t1[sample], color=light_blue, ax=ax[1])
	sns.kdeplot(df_mean_beta_t2[sample], color=light_purple, ax=ax[1])

df_mean_beta_t1['Median'] = df_mean_beta_t1.median(axis=1)
df_mean_beta_t2['Median'] = df_mean_beta_t2.median(axis=1)

sns.kdeplot(df_mean_beta_t1['Median'], color=dark_blue, ax=ax[1])
sns.kdeplot(df_mean_beta_t2['Median'], color=dark_purple, ax=ax[1])

plt.show()