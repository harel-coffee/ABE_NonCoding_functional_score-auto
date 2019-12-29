#!/usr/bin/env python

import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import string
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import scipy
import seaborn as sns
import glob



df = pd.read_csv("Editable_A_scores.combined.scores.csv",index_col=0)
df['chr'] = [x[:-1].split(":")[0] for x in df.index]
df['start'] = [int(x[:-1].split(":")[-1].split("-")[0]) for x in df.index]
df['end'] = [int(x[:-1].split(":")[-1].split("-")[1]) for x in df.index]
df['name'] = [x[:-1] for x in df.index]
df.index = df.name
print (df.head())




names = ['CADD',"DeepSEA",'HbFBase']
for n in names:  
	df['%s_rank'%(n)] = -df[n].rank(ascending=False)

print (df.head())



df = df.sort_values(['chr','start'])
print (df.head())
df = df.reset_index(drop=True)



def to_matrix(df,n,top=500):
	c=n
	size=100
	tmp = df.copy()
	tmp = tmp.sort_values(n,ascending=False).head(n=top)
	my_index_list = tmp.index.tolist()
	out = []
	for i in my_index_list:
		line = []
		current_chr = df.at[i,"chr"]
		for j in range(-size,size+1):
			try:
				chr = df.at[i+j,"chr"]
			except:
				chr = "None"
			if chr == current_chr:
				value = df.at[i+j,"%s_rank"%(c)]
			else:
				value = np.nan
			line.append(value)
		out.append(line)
	out_df = pd.DataFrame(out)
	# print (out_df.head())
	sel_cols = out_df.columns.tolist()
	out_df.index = my_index_list
	out_df['chr'] = tmp['chr']
	out_df['start'] = tmp['start']
	out_df['end'] = tmp['end']
	out_df['name'] = tmp['name']
	out_df['value'] = "."
	out_df['strand'] = "."
	# print (df["%s_rank"%(c)].mean())
	out_df = out_df.fillna(df["%s_rank"%(c)].mean())
	out_df[['chr','start','end','name','value','strand']+sel_cols].to_csv("%s.computeMatrix.bed"%(c),header=False,index=False,sep="\t")
	# print (out_df.head())
	return out_df[sel_cols]	
color_dict ={}
color_dict['HbFBase']='red'
color_dict['CADD']='green'

color_dict['DeepSEA']='blue'
fig, ax = plt.subplots()
for n in names:
	print (n)
	result_df = to_matrix(df,n)
	mean_line = pd.DataFrame(result_df.mean())
	test = pd.melt(result_df)
	# print (test.head())
	sns.lineplot(x="variable", y="value", data=test,c=color_dict[n],ax=ax,label=n)
ax.set_xticklabels(list(range(-100,101,25)))
ax.set_yticklabels(['']+list(range(6000,999,-1000))+[1])
plt.ylabel("Average rank")
plt.xlabel("Downstream / Upstream neighbors")   
plt.xlim(0,201)
plt.savefig("Score_ranks_comparison_top500.pdf", bbox_inches='tight')

