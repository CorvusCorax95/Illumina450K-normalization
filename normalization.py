import numpy as np
import pandas as pd
import qnorm as qn
import scipy.stats as stats
import scipy.optimize as optimize

import prepare_data as prep


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
# want to use the estimate.py from betamix my Schröder, Rahmann
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
    return stats.norm.logpdf(data, loc=params[0], scale=params[1]).sum()


def neglikelihood(params,data):
    return -1*likelihood(params,data)


def bmiq():
    '''STEP 1'''
    '''Fitting of a three-state (unmethylated-U, hemimethylated-H, fully
    methylated-M) beta mixture model to the type1 and type2 probes
    separately. For sake of convenience we refer to intermediate allelic
    methylation as hemimethylation even though hemimethylation
    is most often used in the context of strand-specific methylation.
    -> Realized with betamix (Schroeder, Rahmann)'''
    df_meth, df_unmeth = prep.get_values_as_dataframe_w_types()

    '''Prepping for betamix'''
    #df_beta = prep.beta_value(df_meth, df_unmeth, 100)
    #df_beta.to_csv('df_beta.csv', sep='\t')
    df_beta = pd.read_csv('df_beta.csv', sep='\t', index_col=0)
    df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
    # prep.df_to_h5(df_beta_t1, "type1_probes")
    # prep.df_to_h5(df_beta_t2, "type2_probes")
    probe_list_t1 = df_meth.loc[df_meth["type"] == "I"].index.values.tolist()
    probe_list_t2 = df_meth.loc[df_meth["type"] == "II"].index.values.tolist()
    df_t1_parameters = prep.get_parameters(
        "C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
        "\\betamix-results\\type1_probes-est.h5")
    df_t2_parameters = prep.get_parameters(
        "C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
        "\\betamix-results\\type2_probes-est.h5")
    df_t1_weights = prep.get_weights(
        "C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
        "\\betamix-results\\type1_probes-est.h5", probe_list_t1)
    df_t2_weights = prep.get_weights(
        "C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
        "\\betamix-results\\type2_probes-est.h5", probe_list_t2)
    df_w_class_t1 = prep.get_classes(
        "C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
        "\\betamix-results\\type1_probes-eval.h5", probe_list_t1)
    df_w_class_t2 = prep.get_classes(
        "C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
        "\\betamix-results\\type2_probes-eval.h5", probe_list_t2)

    '''Prep Step 2 & 3'''
    ''' Calculate m2U, then let U2L(U2R) = set of U2 probes with 
    beta-values smaller(larger) than mIIU'''
    ''' Calculate m2M, then let M2L(M2R) = set of M2 probes with 
    beta-values smaller(larger) than m2M'''

    m2U = df_t2_parameters['a']['U'] / (df_t2_parameters['a'][
                                            'U']+df_t2_parameters['b']['U'])
    m2M = df_t2_parameters['a']['M'] / (df_t2_parameters['a'][
                                            'M']+df_t2_parameters['b']['M'])

    '''STEP 2'''
    '''for type2 probes with U-state: transform their probabilities of 
    belonging to the U-state to quantiles using the inverse of the cumulative
    beta-distribution with beta parameters (aU2, bU2)'''
    '''Transform: Inverse Transform sampling?'''
    u2l_list = []
    u2r_list= []
    df_bmiq_U = None
    sample_list = df_beta_t2.columns.values.tolist()[2:]
    for sample in sample_list:
        q_u_list = []
        for probe in probe_list_t2:
            value = df_beta_t2.loc[probe][sample]
            if value <= m2U:
                u2l_list.append(probe)
            else:
                u2r_list.append(probe)
            rb = np.array(df_beta_t2[sample][probe])
            '''p: probability of probe belonging to the U state'''
            p_u = stats.beta.cdf(rb, df_t2_parameters['a']['U'],
                         df_t2_parameters['b']['U'])
            #print("Sample: ", sample, "Probe: ", probe, "pU: ", p_u)
            q_u = stats.beta.ppf(p_u, df_t1_parameters['a']['U'],
                         df_t1_parameters['b']['U'])
            #print("p_U: ", p_u, "q_U: ", p_u)
            q_u_list.append(q_u)
        if df_bmiq_U is None:
            df_bmiq_U = pd.DataFrame(q_u_list, index=probe_list_t2, columns=[
                sample])
        else:
            df_bmiq_U[sample] = q_u_list

    '''STEP 3'''
    '''for type2 probes with M-state: transform their probabilities of 
    belonging to the M-state to quantiles using the inverse of the cumulative
    beta-distribution with beta parameters (aM2, bM2)'''
    m2l_list = []
    m2r_list= []
    df_bmiq_M = None
    for sample in sample_list:
        q_m_list = []
        for probe in probe_list_t2:
            value = df_beta_t2.loc[probe][sample]
            if value <= m2M:
                m2l_list.append(probe)
            else:
                m2r_list.append(probe)
            rb = np.array(df_beta_t2[sample][probe])
            '''p: probability of probe belonging to the U state'''
            p_m = stats.beta.cdf(rb, df_t2_parameters['a']['M'],
                           df_t2_parameters['b']['M'])
            q_m = stats.beta.ppf(p_m, df_t1_parameters['a']['M'],
                           df_t1_parameters['b']['M'])
            q_m_list.append(q_m)
        if df_bmiq_M is None:
            df_bmiq_M = pd.DataFrame(q_m_list, index=probe_list_t2, columns=[
                sample])
        else:
            df_bmiq_M[sample] = q_m_list

    # TODO: Find out why we get primarily 0 and 1 in methylated data
    return df_bmiq_U, df_bmiq_M

    '''STEP 4'''
    '''for type2 probes with H-state: perform a dilation (scale) 
    dransformation to "fit" the data into the "gap" with endpoints defined by
    max(eta2U) and min(eta2M)'''


