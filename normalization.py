import numpy as np
import qnorm as qn
import scipy.stats as scipy
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
    return scipy.norm.logpdf(data, loc=params[0],scale=params[1]).sum()


def neglikelihood(params,data):
    return -1*likelihood(params,data)


def bmiq():
    df_meth, df_unmeth = prep.get_values_as_dataframe_w_types()
    # df_beta = prep.beta_value(df_meth, df_unmeth, 100)
    '''Prepping for betamix'''
    # df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
    # prep.df_to_h5(df_beta_t1, "type1_probes")
    # prep.df_to_h5(df_beta_t2, "type2_probes")
    df = prep.betamix_estimates_to_df(
        "C:\\Users\\lisar\\Documents\\University\\Illumina450K-normalization"
        "\\betamix-results\\type1_probes-est.h5")
    list = prep.get_est_parameters(df, 1)
    print(list)

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

