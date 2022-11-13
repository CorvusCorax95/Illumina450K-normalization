import pandas as pd
from prepare_data import log_data

from plotting import density_plot
from normalization import mean_normalization
from normalization import min_max_normalization
from prepare_data import output_measures
from normalization import quantile_normaliziation

### DATA PREPARATION ###

# from prepare_data import prepare_data
# prepare_data()

# Main holds the dataframes!
#
#

### RAW PLOTTING ###

df_meth = pd.read_csv('resources/short_methylated_w_types.csv', sep='\t')
df_unmeth = pd.read_csv('resources/short_unmethylated_w_types.csv', sep='\t')
df_log_meth = log_data(df_meth, 'logged_meth')
df_log_unmeth = log_data(df_unmeth, 'logged_unmeth')


## RAW
density_plot(df_meth, "Raw plot-methylated")
density_plot(df_unmeth, "Raw plot-unmethylated")
output_measures(df_meth, "Raw Methylated")
output_measures(df_meth, "Raw Unmethylated")

## LOGGED
density_plot(df_log_meth, "Logged plot-methylated")
density_plot(df_log_unmeth, "Logged plot-unmethylated")
output_measures(df_log_meth, "Logged Methylated")
output_measures(df_log_unmeth, "Logged Unmethylated")

## MEAN NORMALIZATION
df_meannorm_meth = mean_normalization(df_log_meth, "mean normalized - methylated")
df_meannorm_unmeth = mean_normalization(df_log_unmeth, "mean normalized - unmethylated")
output_measures(df_meannorm_meth, "Mean Normalized Methylated")
output_measures(df_meannorm_unmeth, "Mean Normalized Unmethylated")

## MIN MAX NORMALIZATION
df_minmax_meth = min_max_normalization(df_log_meth, "min-max normalized - methylated")
df_minmax_unmeth = min_max_normalization(df_log_unmeth, "min-max normalized - unmethylated")
output_measures(df_minmax_meth, "Min-Max Normalized Methylated")
output_measures(df_minmax_unmeth, "Min-Max Normalized Unmethylated")

## QUANTILE NORMALIZATION -> noch kaputt
df_qn_meth = quantile_normaliziation(df_log_meth, "quantile normalized - methylated")
df_qn_unmeth = quantile_normaliziation(df_log_unmeth, "quantile normalized - unmethylated")
output_measures(df_qn_meth, "Quantile Normalized Methylated")
output_measures(df_qn_unmeth, "Quantile Normalized Unmethylated")

#print(df_log_meth.head())
#print(df_qn_meth.head())




print("All done")
