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
    for probe in probe_list_t2:
        value = df_beta_t2.loc[probe][1]
        if value <= m2U:
            u2l_list.append(probe)
        else:
            u2r_list.append(probe)
    rbe001 = np.array(df_beta_t2['RB_E_001'])
    '''p: probability of probe belonging to the U state'''
    p = stats.beta.cdf(rbe001, df_t2_parameters['a']['U'],
                       df_t2_parameters['b']['U'])
    q = stats.beta.ppf(p, df_t1_parameters['a']['U'],
                       df_t1_parameters['b']['U'])
    print(q)
    eta2U = q[0]
    df_bmiq2 = pd.DataFrame(q, columns=['RB_E_001'], index=probe_list_t2)
    return df_bmiq2
    '''STEP 3'''
    '''for type2 probes with M-state: transform their probabilities of 
    belonging to the M-state to quantiles using the inverse of the cumulative
    beta-distribution with beta parameters (aM2, bM2)'''
    m2l_list = []
    m2r_list= []
    for probe in probe_list_t2:
        value = df_beta_t2.loc[probe][1]
        if value <= m2M:
            m2l_list.append(probe)
        else:
            m2r_list.append(probe)

    eta2M = 0 #normalized values of the type2 U-probes
    '''STEP 4'''
    '''for type2 probes with H-state: perform a dilation (scale) 
    dransformation to "fit" the data into the "gap" with endpoints defined by
    max(eta2U) and min(eta2M)'''


# # TODO: UNDER CONSTRUCTION
# # BMIQ
#
#
# # trying to make data work in betamix
# df_beta.to_hdf('data.h5', key='df', mode='w')
#
# # FITTING BETA DISTRIBUTION TO OUR DATA
# m = df_beta.iloc[:, 1:66]  # Beta values in Array without classes
# # transform into numpy array (for fitting)
# naive_array = m.to_numpy()
#
# # fit: Return estimates of shape (if applicable), location, and scale parameters from data. The default
# # estimation method is Maximum Likelihood Estimation (MLE), but Method of Moments (MM) is also available.
# beta_params = stats.beta.fit(naive_array[0])
# # returns probability density function with given parameters (here, parameters from our array)
# st.write("Beta Distribution PDF")
# df_beta_pdf = pd.DataFrame(stats.beta.pdf(naive_array, beta_params[0], beta_params[1],
#                                           beta_params[2], beta_params[3]))
# fig_beta = plot.density_plot(df_beta_pdf, "beta pdf", -5, 10)
# st.pyplot(fig_beta)
#
# # Teschendorff
# # Fitting of a three-state (unmethylated-U, hemimethylated-H, fully
# # methylated-M) beta mixture model to the type1 and type2 probes
# # separately.
#
# # TODO: dataframe df_beta_type1 needs index and column-header
# st.header("Type I Fitting")
# type1 = df_beta[df_beta["type"] == 'I'] #full dataframe with just type I probes
# st.write("Type1: ", type1)
# type1_array = type1.iloc[:, 1:66].to_numpy() #full dataframe as array for further calculations
# st.write("Type1_array: ", type1_array)
# type1_fit = stats.beta.fit(type1_array[0]) #parameters fitting to found beta distribution
# st.write("Type1_fit: ", type1_fit)
# df_beta_type1 = pd.DataFrame(stats.beta.pdf(type1_array, type1_fit[0], type1_fit[1],
#                                             type1_fit[2], type1_fit[3]), columns=sample_list_meth)
# st.write(df_beta_type1)
# fig_beta_type1 = plot.density_plot(df_beta_type1, "beta type 1 pdf", -5, 10)
# st.write("Beta distribution Type 1")
# st.pyplot(fig_beta_type1)
#
# # TODO: dataframe df_beta_type2 needs index and column-header
# st.header("Type II Fitting")
# type2 = df_beta[df_beta["type"] == 'II']
# type2_array = type2.iloc[:, 1:66].to_numpy()
# type2_fit = stats.beta.fit(type2_array[0])
# df_beta_type2 = pd.DataFrame(stats.beta.pdf(type2_array, type2_fit[0], type2_fit[1],
#                                             type2_fit[2], type2_fit[3]))
# st.write(df_beta_type2)
# fig_beta_type2 = plot.density_plot(df_beta_type2, "beta type 2 pdf", -5, 10)
# st.write("Beta distribution Type 2")
# st.pyplot(fig_beta_type2)
#
# # NAIVE CLASSIFICATION WITH CUTOFF POINT
# st.header("Naive classification")
#
# df_states = norm.assign_probes_to_state(df_beta)
# st.write("Naive Classification of states")
# df_t = df_states.transpose()
# df_t.drop('type', inplace=True)
# st.write(df_t)
#
# # TODO: Does not work
# # Trying to get the sample proportions for initial naive values.
# sample_proportions = pd.DataFrame()
#
# first_time = True
# # returns normalized values (probabilities for states)
# for column in df_t:
#     step = df_t[column].value_counts(normalize=True)
#     step.to_frame(name=step.iloc[0])
#
# # def neg_likelihood(params, data):
# #    return -1*scipy.beta.logpdf(data, loc=params[0], scale=params[1]).sum()
#
# # array = np.array([50, 10])
# # result = optimize.minimize(neg_likelihood, array, args=data)
# # st.write(result)

