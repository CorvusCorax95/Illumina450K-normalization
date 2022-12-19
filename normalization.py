import numpy as np
import pandas as pd
import qnorm as qn
import prepare_data


def mean_normalization(df):
    sample_list = df.columns.values.tolist()[2:]
    norm_df = (df[sample_list] - df[sample_list].mean()) / df[sample_list].std()

    return norm_df


def min_max_normalization(df):
    sample_list = df.columns.values.tolist()[2:]
    mm_df = (df[sample_list] - df[sample_list].min()) / (df[sample_list].max() - df[sample_list].min())

    return mm_df


def quantile_normaliziation(df, reference):

    df_qn = qn.quantile_normalize(df, target=df[reference])

    return df_qn

def beta_mixture_normalization(df):
    # 1. Fitting 3-state-beta-mixture models to type I and type II probes separately.
    print("hello")

    # 2. Type-II-probes with U-state

    # 3. Type-II-probes with M-state

    # 4. Type-II-probes with H-state



# idea: copy df, replacing values in df according to methylation state with U, M, H

def assign_probes_to_state(df):

    sample_list = df.columns.values.tolist()[2:]

    df = df.applymap(lambda x: set_states(x) if type(x) == float else x)
    print(df)
    return df

def set_states(x):
    # unmethylated
    if float(x) <= 0.3:
        return 'M'
    # fully-methylated
    elif float(x) >= 0.6:
        return 'U'
    # hemi-methylated
    else:
        return 'H'