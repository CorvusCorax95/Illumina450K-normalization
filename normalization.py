import numpy as np
import qnorm as qn
import scipy.stats as scipy
import scipy.optimize as optimize


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

# TODO: implement!
def beta_mixture_normalization(df):
    # 1. Fitting 3-state-beta-mixture models to type I and type II probes separately.
    print("hello")

    data = df.iloc[:1]
    array = np.array([50, 10])
    result = optimize.minimize(neglikelihood, array, args=data)
    return result
    # 2. Type-II-probes with U-state

    # 3. Type-II-probes with M-state

    # 4. Type-II-probes with H-state



# idea: copy df, replacing values in df according to methylation state with U, M, H

# assign_probes_to_state replaces every value with a given state from set_states
# used to fit a beta distribution on it

def assign_probes_to_state(df):
    df = df.applymap(lambda x: set_states(x) if type(x) == float else x)
    return df

# cutoffs from Schröder/Rahmann 2017
def set_states(x):

    # unmethylated
    if float(x) <= 0.25:
        return 'M'
    # fully-methylated
    elif float(x) >= 0.75:
        return 'U'
    # hemi-methylated
    else:
        return 'H'

    ## I fit 3 different distribution: One beta distribuition for each.

def likelihood(params, data):
    return scipy.norm.logpdf(data, loc=params[0],scale=params[1]).sum()


def neglikelihood(params,data):
    return -1*likelihood(params,data)