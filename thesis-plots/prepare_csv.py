import prepare_data as prep
import normalization as norm

#-------------------------------------
# DEFAULT IS WITH TYPES
# EXCEPT FOR SPLITTED DATAFRAMES
#-------------------------------------

# BASE DATAFRAMES

	# Log Values
df_meth, df_unmeth = prep.get_dataframe(True)
log_m, log_u = prep.log_data()

log_m = prep.add_probetypes(log_m)
log_u = prep.add_probetypes(log_u)

log_m_t1, log_m_t2 = prep.split_types(log_m)
log_u_t1, log_u_t2 = prep.split_types(log_u)

	# Beta Values
beta_values = prep.beta_value(df_meth, df_unmeth, 100)
beta_values = prep.add_probetypes(beta_values)
beta_t1, beta_t2 = prep.split_types(beta_values)

# SAVE BASE DATAFRAMES

log_m.to_csv('log_meth.csv', sep='\t')
log_u.to_csv('log_unmeth.csv', sep='\t')
log_m_t1.to_csv('log_meth_type1.csv', sep='\t')
log_u_t1.to_csv('log_unmeth_type1.csv', sep='\t')
log_m_t2.to_csv('log_meth_type2.csv', sep='\t')
log_u_t2.to_csv('log_unmeth_type2.csv', sep='\t')

beta_values.to_csv('beta_values.csv', sep='\t')
beta_t1.to_csv('beta_values_type1.csv', sep='\t')
beta_t2.to_csv('beta_values_type2.csv', sep='\t')

del log_u['type']
del log_m['type']
del beta_values['type']


# NORMALIZED DATAFRAMES

	# mean
mean_norm_beta = norm.mean_normalization(beta_values)

	# minmax

minmax_norm_beta = norm.min_max_normalization(beta_values)

	# qn

qn_m = log_m
qn_m['Median'] = qn_m.median(axis=1)

qn_u = log_u
qn_u['Median'] = qn_m.median(axis=1)

qn_norm_u = norm.quantile_normalization(qn_u, 'Median')
qn_norm_m = norm.quantile_normalization(qn_m, 'Median')
qn_norm_beta = prep.beta_value(qn_norm_m, qn_norm_u, 100)

	# bmiq
bmiq_norm = norm.bmiq(beta_values)
bmiq_norm.index.name = "probe"
bmiq_norm = prep.add_probetypes(bmiq_norm)
#bmiq_norm_t1, bmiq_norm_t2 = prep.split_types(bmiq_norm)


mean_norm_beta.to_csv('mean_norm_beta.csv', sep='\t')

minmax_norm_beta.to_csv('minmax_norm_beta.csv', sep='\t')

qn_norm_beta.to_csv('qn_norm_beta.csv', sep='\t')

bmiq_norm.to_csv('bmiq_norm.csv', sep='\t')