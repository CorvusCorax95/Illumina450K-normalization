import pandas as pd
from prepare_data import log_data
from prepare_data import beta_value
from prepare_data import m_value

##------------------------------------------------------------------------GET DATA-##
df_meth = pd.read_csv('resources/short_methylated_w_types.csv', sep='\t')
df_unmeth = pd.read_csv('resources/short_unmethylated_w_types.csv', sep='\t')
sample_list_meth = df_meth.columns.values.tolist()[2:]
sample_list_unmeth = df_meth.columns.values.tolist()[2:]
df_unmeth[sample_list_unmeth] = df_unmeth[sample_list_unmeth].astype(int)
df_meth[sample_list_meth] = df_meth[sample_list_meth].apply(lambda x: x + 1)
df_unmeth[sample_list_unmeth] = df_unmeth[sample_list_unmeth].apply(lambda x: x + 1)
df_log_meth = log_data(df_meth, 'logged_meth.csv')
df_log_unmeth = log_data(df_unmeth, 'logged_unmeth.csv')
##---------------------------------------------------------------------------------##

beta_value(df_meth, df_unmeth, 100)