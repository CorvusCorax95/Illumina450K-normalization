import plotting as plot
import prepare_data as prep

df_beta = plot.beta_value()
df_beta.to_csv('df_beta.csv', sep='\t')
df_beta = prep.add_probetypes(df_beta)
df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
prep.df_to_h5(df_beta_t1, "type1_probes")
prep.df_to_h5(df_beta_t2, "type2_probes")