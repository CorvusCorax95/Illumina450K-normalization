import pandas as pd
from plotting import density_plot

import qnorm as qn


def mean_normalization(df, title):
    sample_list = df.columns.values.tolist()[2:]
    norm_df = (df[sample_list] - df[sample_list].mean()) / df[sample_list].std()

    density_plot(norm_df, title)
    return norm_df

def min_max_normalization(df, title):
    sample_list = df.columns.values.tolist()[2:]
    mm_df = (df[sample_list] - df[sample_list].min()) / (df[sample_list].max() - df[sample_list].min())

    density_plot(mm_df, title)

    return mm_df

def quantile_normaliziation(df, title):
    #sample_list = df.columns.values.tolist()[2:]
    df_qn = qn.quantile_normalize(df, axis=1)
    density_plot(df_qn, title)
    return df_qn