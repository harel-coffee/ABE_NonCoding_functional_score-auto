

"""
Program to get average editing frequency from CrispEsso.

"""

import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import string
import glob
import numpy as np

def revcomp(seq):
	try: ## python2
		tab = string.maketrans(b"ACTG", b"TGAC")
	except:  ## python3
		tab = bytes.maketrans(b"ACTG", b"TGAC")
	return seq.translate(tab)[::-1]

def get_ref_alt(input_gRNA,table_header_list):
	table_header=[x.split(".")[0] for x in table_header_list]
	table_header = "".join(table_header)
	if input_gRNA == table_header:
		return "A","G"
	elif revcomp(input_gRNA) == table_header:
		return "T","C"
	else:
		print ("something is wrong, Exist!")
		exit()
	

def parse_df(x):
	freq_table = "/".join(x.split("/")[:-1])+"/Quantification_window_nucleotide_percentage_table.txt"
	gRNA = x.split("_")[-1].replace(".txt","")
	df = pd.read_csv(freq_table,sep="\t",index_col=0)
	df[df < 0.01] = 0
	my_freq = []
	flag_type = []
	ref,alt = get_ref_alt(gRNA,df.columns.tolist())
	for c in df.columns:
		base = c.split(".")[0]
		freq = 0
		if base == ref:
			freq = df.at[alt,c]
		my_freq.append(freq)
	if ref == "T":
		return my_freq[::-1]
	return my_freq
		
			
			
		

	

inputs = pd.read_csv("input.list",header=None)
inputs = inputs[0].tolist()


df = pd.DataFrame([parse_df(i) for i in inputs])


df.index = inputs
# print (df.head())
df = df.replace(0, np.NaN)
df.to_csv("average_model.csv")

avg = pd.DataFrame(df.mean())
avg[1] = avg.index.tolist()
print (avg)

plt.figure(figsize=(5,2))

sns.barplot(x=1,y=0,data=avg)
# plt.xticks(rotation=90)
plt.ylabel("Editing frequency")
plt.xlabel('A position')
plt.rcParams.update({'font.size': 10})
plt.savefig("editing_frequecy_barplot.pdf", bbox_inches='tight')



