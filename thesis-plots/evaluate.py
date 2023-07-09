import pandas as pd

import normalization as norm

df_beta = pd.read_csv('beta_values.csv', sep='\t', index_col=0)

df_beta_t1 = pd.read_csv('beta_values_type1.csv', sep='\t', index_col=0)
df_beta_t2 = pd.read_csv('beta_values_type2.csv', sep='\t', index_col=0)

df_minmax_m = pd.read_csv('minmax_norm_m.csv', sep='\t', index_col=0)
df_minmax_u = pd.read_csv('minmax_norm_u.csv', sep='\t', index_col=0)
df_minmax_u_t1 = pd.read_csv('minmax_norm_u_t1.csv', sep='\t', index_col=0)
df_minmax_u_t2 = pd.read_csv('minmax_norm_u_t2.csv', sep='\t', index_col=0)
df_minmax_m_t1 = pd.read_csv('minmax_norm_m_t1.csv', sep='\t', index_col=0)
df_minmax_m_t2 = pd.read_csv('minmax_norm_m_t2.csv', sep='\t', index_col=0)

df_mean_m = pd.read_csv('mean_norm_m.csv', sep='\t', index_col=0)
df_mean_u = pd.read_csv('mean_norm_u.csv', sep='\t', index_col=0)
df_mean_u_t1 = pd.read_csv('mean_norm_u_t1.csv', sep='\t', index_col=0)
df_mean_u_t2 = pd.read_csv('mean_norm_u_t2.csv', sep='\t', index_col=0)
df_mean_m_t1 = pd.read_csv('mean_norm_m_t1.csv', sep='\t', index_col=0)
df_mean_m_t2 = pd.read_csv('mean_norm_m_t2.csv', sep='\t', index_col=0)

df_qn_m = pd.read_csv('qn_norm_m.csv', sep='\t', index_col=0)
df_qn_u = pd.read_csv('qn_norm_u.csv', sep='\t', index_col=0)
df_qn_u_t1 = pd.read_csv('qn_norm_u_t1.csv', sep='\t', index_col=0)
df_qn_u_t2 = pd.read_csv('qn_norm_u_t2.csv', sep='\t', index_col=0)
df_qn_m_t1 = pd.read_csv('qn_norm_m_t1.csv', sep='\t', index_col=0)
df_qn_m_t2 = pd.read_csv('qn_norm_m_t2.csv', sep='\t', index_col=0)

df_beta_t1['Median'] = df_beta_t1.median(axis=1)
df_beta_t2['Median'] = df_beta_t2.median(axis=1)

df_minmax_beta_t1 = norm.min_max_normalization(df_beta_t1)
df_minmax_beta_t2 = norm.min_max_normalization(df_beta_t2)
df_mean_beta_t1 = norm.mean_normalization(df_beta_t1)
df_mean_beta_t2 = norm.mean_normalization(df_beta_t2)
df_qn_beta_t1 = norm.quantile_normaliziation(df_beta_t1, 'Median')
df_qn_beta_t2 = norm.quantile_normaliziation(df_beta_t2, 'Median')


lst = [df_beta, df_mean_u, df_mean_m, df_minmax_u, df_minmax_m, df_qn_u, df_qn_m]


for df in lst:
	df['mean'] = df.mean(axis=1)
	df['median'] = df.median(axis=1)
	df['variance'] = df.var(axis=1)

	count_I = 0
	count_II = 0
probe_list = df_beta.index.to_list()
for df in lst:
	count_I = 0
	count_II = 0
	for probe in probe_list:
		if df['type'] == 'I' and count_I < 5:
			print(df['RB_E_001']['type'])
			count_I = count_I + 1
		if df['RB_E_001']['type'] == 'II' and count_II < 5:
			print(df['RB_E_001']['type'])
			count_II = count_II + 1