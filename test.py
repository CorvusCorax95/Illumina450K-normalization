import pandas as pd
import prepare_data as prep

##------------------------------------------------------------------------GET DATA-##
df_meth = pd.read_csv('resources/short_methylated_w_types.csv', sep='\t')
df_unmeth = pd.read_csv('resources/short_unmethylated_w_types.csv', sep='\t')
sample_list_meth = df_meth.columns.values.tolist()[2:]
sample_list_unmeth = df_meth.columns.values.tolist()[2:]
df_unmeth[sample_list_unmeth] = df_unmeth[sample_list_unmeth].astype(int)
df_meth[sample_list_meth] = df_meth[sample_list_meth].apply(lambda x: x + 1)
df_unmeth[sample_list_unmeth] = df_unmeth[sample_list_unmeth].apply(lambda x: x + 1)
df_log_meth = prep.log_data(df_meth, 'logged_meth.csv')
df_log_unmeth = prep.log_data(df_unmeth, 'logged_unmeth.csv')
df_beta = prep.beta_value(df_meth, df_unmeth, 100)
##---------------------------------------------------------------------------------##

import numpy as np

from scipy.stats import beta

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

a, b = 2.31, 0.627

mean, var, skew, kurt = beta.stats(a, b, moments='mvsk')

x = np.linspace(beta.ppf(0.01, a, b), beta.ppf(0.99, a, b), 100)

print(df_meth)
print(df_beta)
print(x)

ax.plot(x, beta.pdf(x, a, b),

        'r-', lw=5, alpha=0.6, label='beta pdf')

rv = beta(a, b)

ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

vals = beta.ppf([0.001, 0.5, 0.999], a, b)

np.allclose([0.001, 0.5, 0.999], beta.cdf(vals, a, b))
True

r = beta.rvs(a, b, size=1000)

ax.hist(r, density=True, bins='auto', histtype='stepfilled', alpha=0.2)

ax.set_xlim([x[0], x[-1]])

ax.legend(loc='best', frameon=False)

plt.show()
