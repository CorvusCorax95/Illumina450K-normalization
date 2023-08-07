import numpy as np
import pandas as pd
from h5py import h5

import prepare_data as prep




def df_to_h5(df, filename):
	"""As preparation to use the betamix-tool by Schröder, Rahmann the
	dataframe is converted to a fitting h5 file."""
	# delete existing file first by hand!
	sample_list = df.columns.values.tolist()[1:]
	df.reset_index(drop=True)
	filename = filename + ".h5"
	hf = h5.File(filename, "w")
	group = hf.create_group("data")
	i = 0
	for sample in sample_list:
		array = np.array(df[sample], dtype=np.float64)
		#dset = group.create_dataset(str(i), data=np.sort(array))
		i = i + 1
	hf.close()


# betamix commands
# python estimate.py -F -t 1E-5 mytestfile.h5 --resultpath mytestfile-est.h5
# python evaluate.py -F -i snps -f pdf mytestfile.h5 mytestfile-est.h5 mytestfile-eval.h5

def get_parameters(path, n):
	'''After using the betamix-tool we need to get the esimated parameters
	from the h5 file provided.'''
	num = str(n)
	f = h5.File(path, "r")
	grp = f['estimation']
	subgrp = grp[num]
	data = subgrp["ab"]
	return pd.DataFrame(data, index=["U", "H", "M"], columns=["a", "b"])


def get_classes(path, probe_list, n):
	'''gets the class to probe and returns dataframe'''
	'''0: unmethylated, 1: hemimethylated, 2: methylated'''
	num = str(n)
	f = h5.File(path, "r")
	grp = f['evaluation']
	subgrp = grp[num]
	subsubgrp = subgrp['mix']
	data = subsubgrp['classes']
	return pd.DataFrame(data, index=probe_list, columns=["sample"])

def bmiq():
	'''STEP 1'''
	'''Fitting of a three-state (unmethylated-U, hemimethylated-H, fully
	methylated-M) beta mixture model to the type1 and type2 probes
	separately. For sake of convenience we refer to intermediate allelic
	methylation as hemimethylation even though hemimethylation
	is most often used in the context of strand-specific methylation.
	-> Realized with betamix (Schroeder, Rahmann)'''
	df_meth, df_unmeth = prep.get_dataframe(True)

	'''Prepping for betamix'''
	df_beta = prep.beta_value(df_meth, df_unmeth, 100)
	df_beta.to_csv('df_beta.csv', sep='\t')
	df_beta_t1, df_beta_t2 = prep.split_types(df_beta)
	sample_list = df_beta_t2.columns.values.tolist()
	prep.df_to_h5(df_beta_t1, "type1_probes")
	prep.df_to_h5(df_beta_t2, "type2_probes")
	'''list with all names of probes with type 2'''
	'''[...]'''