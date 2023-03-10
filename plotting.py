
import matplotlib.pyplot as plt

import streamlit as st

from normalization import *


def density_plot(df, title, x1, x2):
    plt.style.use('dark_background')
    df.plot.density(linewidth=1, figsize=(20, 10), xlim=(x1, x2))
    plt.legend([])
    plt.title(title)

def methylated_plots(df_log_meth):
    st.write("Methylated Plots")
    fig_dens_meth = density_plot(df_log_meth, "Logged plot - methylated", 5, 17)
    st.pyplot(fig_dens_meth)
    # MEAN NORMALIZED
    df_meannorm_meth = mean_normalization(df_log_meth)
    fig_meannorm_meth = density_plot(df_meannorm_meth, "Mean Normalized - methylated", -4, 3)
    st.pyplot(fig_meannorm_meth)
    ## MIN MAX NORMALIZED
    df_minmax_meth = min_max_normalization(df_log_meth)
    fig_mm_meth = density_plot(df_minmax_meth, "Min-Max Normalized - methylated", -0.2, 1.25)
    st.pyplot(fig_mm_meth)
    ## QUANTILE NORMALIZED
    # last column as Median column
    df_log_meth['Median'] = df_log_meth.median(axis=1)
    reference_options = df_log_meth.columns.values.tolist()[2:]
    reference_meth = st.selectbox('Which reference do you want to use?',  reference_options, (len(reference_options)-1), key=0)
    df_qn_meth = quantile_normaliziation(df_log_meth, reference_meth)
    fig_qn_meth = density_plot(df_qn_meth, "Quantile Normalized plot (Median) - methylated", 5, 17)
    st.pyplot(fig_qn_meth)

def unmethylated_plots(df_log_unmeth):
    st.write("Unmethylated Plots")
    fig_dens_unmeth = density_plot(df_log_unmeth, "Logged plot - unmethylated", 5, 17)
    st.pyplot(fig_dens_unmeth)
    # MEAN NORMALIZED
    df_meannorm_unmeth = mean_normalization(df_log_unmeth)
    fig_meannorm_unmeth = density_plot(df_meannorm_unmeth, "Mean Normalized - unmethylated", -4, 3)
    st.pyplot(fig_meannorm_unmeth)
    ## MIN MAX NORMALIZED
    df_minmax_unmeth = min_max_normalization(df_log_unmeth)
    fig_mm_unmeth = density_plot(df_minmax_unmeth, "Min-Max Normalized - unmethylated", -0.2, 1.25)
    st.pyplot(fig_mm_unmeth)
    ## QUANTILE NORMALIZED
    # last column as Median column
    df_log_unmeth['Median'] = df_log_unmeth.median(axis=1)
    reference_options = df_log_unmeth.columns.values.tolist()[2:]
    reference_unmeth = st.selectbox('Which reference do you want to use?', reference_options, (len(reference_options)-1), key=1)
    df_qn_unmeth = quantile_normaliziation(df_log_unmeth, reference_unmeth)
    fig_qn_unmeth = density_plot(df_qn_unmeth, "Quantile Normalized plot (Median) - unmethylated", 5, 17)
    st.pyplot(fig_qn_unmeth)

def beta_value_plots(df_beta):
    # showing off the beta values
    fig_beta = density_plot(df_beta, "Beta Values", -0.5, 1.5)
    st.header("Beta Values")
    st.write(df_beta)
    st.write("Beta Values (Methylation Values)")
    st.pyplot(fig_beta)