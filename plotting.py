import pandas as pd
import matplotlib.pyplot as plt
from prepare_data import log_data
from prepare_data import make_dataframe


# Kann ich beim plotten nur die "inneren" 95% nehmen um Ausreißer zu ignorieren?

def density_plot(df, title, x1, x2):
    plt.style.use('dark_background')
    df.plot.density(linewidth=1, figsize=(20, 10), xlim=(x1, x2))
    plt.legend([])
    plt.title(title)