import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as scipy

import normalization
import plotting
from prepare_data import log_data

from plotting import density_plot
from normalization import mean_normalization
from normalization import min_max_normalization
from prepare_data import output_measures
from normalization import quantile_normaliziation
from normalization import assign_probes_to_state
from prepare_data import beta_value
from prepare_data import m_value
import scipy.optimize as optimize

import streamlit as st


### DATA PREPARATION ###

# from prepare_data import prepare_data
# prepare_data()

# Main holds the dataframes!
#
def make_header():
    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.set_page_config(page_title="450K normalization", page_icon=":chart_with_upwards_trend:", layout="wide")

    ## HEADER ##
    with st.container():  # for wrapping contents
        st.header("Normalization of Illumina HumanMethylation 450K Beadchips")
        st.subheader("Methylation normalization because technical variability sucks.")

    ## TABLE ##
    with st.container():
        left_column, right_column = st.columns(2)
        with left_column:
            st.header("Who I am")
            st.write("""
            Lisa
            - Bioinformatics Student (Bachelor of Science)
            - struggling my way through university
            """
                     )
        with right_column:
            st.header("My University")
            st.write("""
            Universität des Saarlandes
            - Center of Bioinformatics
            - Chair of Algorithmic Bioinformatics
            - university in a land just used for size comparisons
            """
                     )


def make_plots():
    with st.container():
        left_column, right_column = st.columns(2)
        with left_column:
            # TODO: DONE
            #plotting.methylated_plots(df_log_meth)

            ## BETA MIXTURE NORMALIZED
            df_beta = beta_value(df_meth, df_unmeth, 100)
            fig_beta = density_plot(df_beta, "Beta Values", -0.5, 1.5)
            st.pyplot(fig_beta)

            # TODO: UNDER CONSTRUCTION
            #single_col = df_beta.iloc[:, 1]
            #array = single_col.to_numpy()
            #st.write(df_beta)
            #st.write(df_beta.iloc[:, 1])

        with right_column:
            # TODO: DONE
            #plotting.unmethylated_plots(df_log_unmeth)

            ## BETA MIXTURE NORMALIZED
            df_m = m_value(df_meth, df_unmeth)
            fig_m = density_plot(df_m, "M Values", -25, -5)
            st.pyplot(fig_m)

            # TODO: UNDER CONSTRUCTION
            #df_states = assign_probes_to_state(df_beta)
            #st.write(df_states)
            #st.write(df_beta)
            #st.write(df_beta.iloc[:, 1])

        #def neglikelihood(params, data):
        #    return -1*scipy.beta.logpdf(data, loc=params[0], scale=params[1]).sum()

        #array = np.array([50, 10])
        #result = optimize.minimize(neglikelihood, array, args=data)
        #st.write(result)


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


make_header()
make_plots()
