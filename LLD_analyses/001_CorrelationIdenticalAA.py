#This script will check for peptides with identical sequence and correlate the presence/absence profiles.
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from numpy import cov
from numpy import std
from scipy.stats import pearsonr
import matplotlib.pyplot as plt 



#Need to define a pickle with the 0/1 presence absence matrix (rows = samples, columns = peptides)
Pickle_matrix = ""
Enrichment_data = pd.read_pickle(Pickle_matrix)
#Need to define a file with annotation per peptide, including sequence in pythonic position 3 (4rth column)
Annotation_file = ""


def get_correlation(X, Y):
	'Compute Pearson correlation from scipy between two vectors'
	corr, _ = pearsonr(X, Y) 
	return(corr)
def Run_only_shigella():
	'This is an example of using only'
	#Get only peptides belonign to Shigella
	Shigella = ["agilent_173099","agilent_141247", "agilent_120811", "agilent_211821", "agilent_235763", "agilent_221918", "agilent_113080", "agilent_230549","agilent_124313","agilent_71908","agilent_231662","agilent_221258","agilent_24468"]
	ED = Enrichment_data[Shigella]
	CORR = ED.corr()
	Select = pd.np.triu(pd.np.ones(CORR.shape), k=1).astype(bool)
	CORR = CORR.where(Select).stack().reset_index()
	M  =CORR.iloc[:,2].mean()
	Max = CORR.iloc[:,2].max()
	Min = CORR.iloc[:,2].min()
	print(M, Max, Min)
	exit()

dic_protein = {}
#Correlate identical proteins
#Need to read annotation file to identificate which peptides have identical 
with open(Annotation_file) as F:
	for line in F:
		l = line.rstrip().split(",")
		if l[0]== "order": continue
		Seq = l[3]
		if "(" in Seq:
			if Seq.split()[0] not in dic_protein:  dic_protein[Seq.split()[0]] = []
			dic_protein[Seq.split()[0]].append(l[0])
		else:	continue

#Iterate thorugh the different peptide groups and compute correlation			
Means = []
Total = None
for id_seq in dic_protein:
	Prots = dic_protein[id_seq]
	if len(Prots) < 2: continue
	DF = Enrichment_data[Enrichment_data.columns.intersection(Prots)]
	if DF.shape[1] < 2: continue 
	
	CORR = DF.corr()
	Select = pd.np.triu(pd.np.ones(CORR.shape), k=1).astype(bool)
	CORR = CORR.where(Select).stack().reset_index()
	M  =CORR.iloc[:,2].mean()
	Means.append(M)
	if isinstance(Total, pd.DataFrame):
		Total = Total.append(CORR)
	else: Total = CORR
	
Total.to_csv("Results/Corr_sameSeq.tsv", sep="\t")
fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
# We can set the number of bins with the `bins` kwarg
axs.hist(Means, bins=20)
plt.savefig("Results/Mean_corr_sameSeq.png")


