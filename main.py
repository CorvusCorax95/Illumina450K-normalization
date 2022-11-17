import matplotlib.pyplot as plt
import pandas as pd
from prepare_data import log_data

from plotting import density_plot
from normalization import mean_normalization
from normalization import min_max_normalization
from prepare_data import output_measures
from normalization import quantile_normaliziation

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
            st.write("Methylated Plots")
            fig_dens_meth = density_plot(df_log_meth, "Logged plot - methylated", 5, 17)
            st.pyplot(fig_dens_meth)
            ## MEAN NORMALIZED
            df_meannorm_meth = mean_normalization(df_log_meth)
            fig_meannorm_meth = density_plot(df_meannorm_meth, "Mean Normalized - methylated", -4, 3)
            st.pyplot(fig_meannorm_meth)
            ## MIN MAX NORMALIZED
            df_minmax_meth = min_max_normalization(df_log_meth)
            fig_mm_meth = density_plot(df_minmax_meth, "Min-Max Normalized - methylated", -0.2, 1.25)
            st.pyplot(fig_mm_meth)
            ## QUANTILE NORMALIZED
            df_qn_meth = quantile_normaliziation(df_log_meth)
            fig_qn_meth = density_plot(df_qn_meth, "Quantile Normalized plot (Median) - methylated", 5, 17)
            st.pyplot(fig_qn_meth)


        with right_column:
            st.write("Unmethylated Plots")
            fig_dens_unmeth = density_plot(df_log_unmeth, "Logged plot - unmethylated", 5, 17)
            st.pyplot(fig_dens_unmeth)
            ## MEAN NORMALIZED
            df_meannorm_unmeth = mean_normalization(df_log_unmeth)
            fig_meannorm_unmeth = density_plot(df_meannorm_unmeth, "Mean Normalized - unmethylated", -4, 3)
            st.pyplot(fig_meannorm_unmeth)
            ## MIN MAX NORMALIZED
            df_minmax_unmeth = min_max_normalization(df_log_unmeth)
            fig_mm_unmeth = density_plot(df_minmax_unmeth, "Min-Max Normalized - unmethylated", -0.2, 1.25)
            st.pyplot(fig_mm_unmeth)
            ## QUANTILE NORMALIZED
            df_qn_unmeth = quantile_normaliziation(df_log_unmeth)
            fig_qn_unmeth = density_plot(df_qn_unmeth, "Quantile Normalized plot (Median) - unmethylated", 5, 17)
            st.pyplot(fig_qn_unmeth)

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