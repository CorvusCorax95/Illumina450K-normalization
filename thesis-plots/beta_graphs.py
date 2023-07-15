import pandas as pd

import prepare_data as prep
import matplotlib.pyplot as plt
import seaborn as sns

### SEABORN ###

# makes two side-by-side log value densities
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

light_purple = (0.678, 0.596, 0.839, 1)  # beta
dark_purple = (0.212, 0.047, 0.541, 1)  # beta
hex_purple = "#360c8a"

# ------------------------------------------------------------------#
# PREP WORK
df_beta = pd.read_csv('beta_values.csv', sep='\t', index_col=0)
df_t1 = pd.read_csv('beta_values_type1.csv', sep='\t', index_col=0)
df_t2 = pd.read_csv('beta_values_type2.csv', sep='\t', index_col=0)
df_beta_norm = pd.read_csv('bmiq_norm.csv', sep='\t', index_col=0)
df_beta_t1 = df_t1.copy()
df_beta_t2 = df_t2.copy()


del df_beta['type']
del df_beta_norm['type']

df_beta_t1['Median'] = df_beta_t1.median(axis=1)
df_beta_t2['Median'] = df_beta_t2.median(axis=1)
df_beta['Median'] = df_beta.median(axis=1)
df_beta_norm['Median'] = df_beta_norm.median(axis=1)
#------------------------------------------------------------------#
fig_beta, ax = plt.subplots(2)
fig_beta.suptitle("Beta Value Density Plots")
ax[0].set_title("A) Beta Values", loc='left')
ax[0].set(xlabel="beta value")
ax[1].set_title("B) Beta Values Type 1 & 2", loc="left")
ax[1].set(xlabel="beta value")
sample_list = df_beta.columns.values.tolist()[1:]
for sample in sample_list:
	sns.kdeplot(df_beta[sample], color=light_blue, ax=ax[0])
	sns.kdeplot(df_beta_t1[sample], color=light_red, ax=ax[1])
	sns.kdeplot(df_beta_t2[sample], color=light_yellow, ax=ax[1])

sns.kdeplot(df_beta_t1['Median'], color=dark_red, ax=ax[1])
sns.kdeplot(df_beta_t2['Median'], color=dark_yellow, ax=ax[1])
sns.kdeplot(df_beta['Median'], color=dark_blue, ax=ax[0])

lst = ('beta-values', 'beta-values type 1', 'beta-values type 2')
legend = fig_beta.legend(labels=lst, frameon=False,
                                  loc='lower center', ncol=5)

handles = legend.legend_handles

colors = [hex_purple, hex_red, hex_yellow]

for i, handle in enumerate(handles):
	handle.set_color(colors[i])


plt.show()