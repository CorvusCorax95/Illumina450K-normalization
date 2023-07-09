import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

df_meth_wo_types = pd.read_csv('log_meth_wo_types.csv', sep='\t',
                               index_col='probe')
df_unmeth_wo_types = pd.read_csv('log_unmeth_wo_types.csv', sep='\t',
                                 index_col='probe')
df_meth = pd.read_csv('log_meth.csv', sep='\t', index_col='probe')
df_meth = pd.read_csv('log_unmeth.csv', sep='\t', index_col='probe')

df_t1_m = pd.read_csv('log_meth_type1.csv', sep='\t', index_col='probe')
df_t2_m = pd.read_csv('log_meth_type2.csv', sep='\t', index_col='probe')

df_t1_u = pd.read_csv('log_unmeth_type1.csv', sep='\t', index_col='probe')
df_t2_u = pd.read_csv('log_unmeth_type2.csv', sep='\t', index_col='probe')

df_t1_m['Median'] = df_t1_m.median(axis=1)
df_t2_m['Median'] = df_t2_m.median(axis=1)
df_t1_u['Median'] = df_t1_u.median(axis=1)
df_t2_u['Median'] = df_t2_u.median(axis=1)

# show the plot
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

fig_log_densities, ax = plt.subplots(2, 3)
fig_log_densities.suptitle("logarithmic density functions")

ax[0, 0].set_title("A) Methylated density plots", loc='left')
ax[0, 0].set_xlabel("log(2)-brightness")

ax[1, 0].set_title("B) Unmethylated density plots", loc='left')
ax[1, 0].set_xlabel("log(2)-brightness")

ax[0, 1].set_title("C) Type 1 density plots", loc='left')
ax[0, 1].set_xlabel("log(2)-brightness")

ax[1, 1].set_title("D) Type 2 density plots", loc='left')
ax[1, 1].set_xlabel("log(2)-brightness")

ax[0, 2].set_title("E) Boxplots", loc='left')

ax[1, 2].set_title("F) Median density plots", loc='left')

sample_list = df_t1_m.columns.values.tolist()[1:]
for sample in sample_list:
	# METHYLATED
	sns.kdeplot(df_t1_m[sample], color=light_yellow, ax=ax[0, 0])
	sns.kdeplot(df_t2_m[sample], color=light_green, ax=ax[0, 0])
	# UNMETHYLATED
	sns.kdeplot(df_t1_u[sample], color=light_blue, ax=ax[1, 0])
	sns.kdeplot(df_t2_u[sample], color=light_red, ax=ax[1, 0])

# METHYLATED MEDIAN
sns.kdeplot(df_t1_m['Median'], color=dark_yellow, ax=ax[0, 0], legend=True)
sns.kdeplot(df_t2_m['Median'], color=dark_green, ax=ax[0, 0], legend=True)
# UNMETHYLATED MEDIAN
sns.kdeplot(df_t1_u['Median'], color=dark_blue, ax=ax[1, 0], legend=True)
sns.kdeplot(df_t2_u['Median'], color=dark_red, ax=ax[1, 0], legend=True)

for sample in sample_list:
	# TYPE 1
	sns.kdeplot(df_t1_m[sample], color=light_yellow, ax=ax[0, 1])
	sns.kdeplot(df_t1_u[sample], color=light_blue, ax=ax[0, 1])
	# TYPE 2
	sns.kdeplot(df_t2_m[sample], color=light_green, ax=ax[1, 1])
	sns.kdeplot(df_t2_u[sample], color=light_red, ax=ax[1, 1])

# TYPE 1 MEDIAN
sns.kdeplot(df_t1_m['Median'], color=dark_yellow, ax=ax[0, 1], legend=True)
sns.kdeplot(df_t1_u['Median'], color=dark_blue, ax=ax[0, 1], legend=True)
# TYPE 2 MEDIAN
sns.kdeplot(df_t2_u['Median'], color=dark_red, ax=ax[1, 1], legend=True)
sns.kdeplot(df_t2_m['Median'], color=dark_green, ax=ax[1, 1], legend=True)

# TYPE 1 MEDIAN
sns.kdeplot(df_t1_m['Median'], color=dark_yellow, ax=ax[1, 2], legend=True)
sns.kdeplot(df_t1_u['Median'], color=dark_blue, ax=ax[1, 2], legend=True)
# TYPE 2 MEDIAN
sns.kdeplot(df_t2_u['Median'], color=dark_red, ax=ax[1, 2], legend=True)
sns.kdeplot(df_t2_m['Median'], color=dark_green, ax=ax[1, 2], legend=True)

# ------------------------------------------------------------------#

sns.boxplot(ax=ax[0, 2], data=[df_t1_u['Median'], df_t1_m['Median'], df_t2_m[
	'Median'], df_t2_u['Median']])

lst = ('Unmethylated Type 1', 'Methylated Type 1', 'Unmethylated Type 2',
       'Methylated Type 2')
legend = fig_log_densities.legend(labels=lst, frameon=False,
                                  loc='lower center', ncol=4)

handles = legend.legend_handles

colors = [hex_blue, hex_yellow, hex_red, hex_green]

for i, handle in enumerate(handles):
	handle.set_color(colors[i])

cap = "text"

plt.show()
