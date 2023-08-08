
import prepare_data as prep


# Ansatz funktioniert nicht
# df_meth, df_unmeth = prep.log_data()
#
# df_meth['Median'] = df_meth.median(axis=1)
# df_unmeth['Median'] = df_unmeth.median(axis=1)
#
# df_qn_meth = qn.quantile_normalize(df_meth, target=df_meth['Median'])
# df_qn_unmeth = qn.quantile_normalize(df_unmeth, target=df_unmeth['Median'])
# df_beta = prep.beta_value(df_qn_meth, df_qn_unmeth, 100)
df_meth, df_unmeth = prep.get_dataframe(False)
df_beta = prep.beta_value(df_meth, df_unmeth, 100)

df_beta.to_csv('df_beta.csv', sep='\t')
df_beta = prep.add_probetypes(df_beta)
df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
prep.df_to_h5(df_beta_t1, "type1_probes")
prep.df_to_h5(df_beta_t2, "type2_probes")