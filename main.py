import pandas as pd
from prepare_data import log_data

from plotting import density_plot
from normalization import mean_normalization
from normalization import min_max_normalization
from prepare_data import output_measures
from normalization import quantile_normaliziation
import plotly.express as px

import streamlit as st

### DATA PREPARATION ###

# from prepare_data import prepare_data
# prepare_data()

# Main holds the dataframes!
#



def do_stuff():

    sample_list = df_meth.columns.values.tolist()[2:]
    df_meth[sample_list] = df_meth[sample_list].apply(lambda x: x + 1)
    df_unmeth[sample_list] = df_unmeth[sample_list].apply(lambda x: x + 1)
    df_log_meth = log_data(df_meth, 'logged_meth.csv')
    df_log_unmeth = log_data(df_unmeth, 'logged_unmeth.csv')


    ## RAW
    #density_plot(df_meth, "Raw plot-methylated")
    #density_plot(df_unmeth, "Raw plot-unmethylated")
    #output_measures(df_meth, "Raw Methylated")
    #output_measures(df_meth, "Raw Unmethylated")

    ## LOGGED
    # density_plot(df_log_meth, "Logged plot-methylated")
    # density_plot(df_log_unmeth, "Logged plot-unmethylated")
    # output_measures(df_log_meth, "Logged Methylated")
    # output_measures(df_log_unmeth, "Logged Unmethylated")
    #
    # ## MEAN NORMALIZATION
    # df_meannorm_meth = mean_normalization(df_log_meth, "mean normalized - methylated")
    # df_meannorm_unmeth = mean_normalization(df_log_unmeth, "mean normalized - unmethylated")
    # output_measures(df_meannorm_meth, "Mean Normalized Methylated")
    # output_measures(df_meannorm_unmeth, "Mean Normalized Unmethylated")
    #
    # ## MIN MAX NORMALIZATION
    # df_minmax_meth = min_max_normalization(df_log_meth, "min-max normalized - methylated")
    # df_minmax_unmeth = min_max_normalization(df_log_unmeth, "min-max normalized - unmethylated")
    # output_measures(df_minmax_meth, "Min-Max Normalized Methylated")
    # output_measures(df_minmax_unmeth, "Min-Max Normalized Unmethylated")
    #
    # ## QUANTILE NORMALIZATION -> noch kaputt
    # df_qn_meth = quantile_normaliziation(df_log_meth, "quantile normalized - methylated")
    # df_qn_unmeth = quantile_normaliziation(df_log_unmeth, "quantile normalized - unmethylated")
    # output_measures(df_qn_meth, "Quantile Normalized Methylated")
    # output_measures(df_qn_unmeth, "Quantile Normalized Unmethylated")

    print("All done")


def build_website():
    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.set_page_config(page_title="450K normalization", page_icon=":chart_with_upwards_trend:", layout="wide")

    ## HEADER ##
    with st.container(): #for wrapping contents
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
    #with st.container():
    #    st.write(df_meth)
    #    st.write(df_unmeth)

    with st.container():
        st.write("Methylated Plots")
        left_column, right_column = st.columns(2)
        with left_column:
            fig_dens_meth = density_plot(df_log_meth, "Logged plot - methylated")
            st.pyplot(fig_dens_meth)

        with right_column:
            df_meannorm_meth = mean_normalization(df_log_meth)
            fig_meannorm_meth = density_plot(df_meannorm_meth, "Mean Normalized - methylated")
            st.pyplot(fig_meannorm_meth)


    with st.container():
        st.write("Unmethylated Plots")
        left_column, right_column = st.columns(2)
        with left_column:
            fig_dens_unmeth = density_plot(df_log_unmeth, "Logged plot - unmethylated")
            st.pyplot(fig_dens_unmeth)

        with right_column:
            df_meannorm_unmeth = mean_normalization(df_log_unmeth)
            fig_meannorm_unmeth = density_plot(df_meannorm_unmeth, "Mean Normalized - unmethylated")
            st.pyplot(fig_meannorm_unmeth)


#do_stuff()

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


build_website()