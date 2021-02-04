# -*- coding: utf-8 -*-

# exec(open("barplot_NG_revision.py").read())
import pandas as pd 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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
import scipy
import numpy as np

color_dict = {}
# color_dict['RF:all']="#fa2111"
# color_dict['RF:TFBS']="#fa8911"
# color_dict['RF:Epi']="#fc42a5"
# color_dict['DeepSEA']="#3e5df7"
# color_dict['DeepSEA']="#213fff"
# color_dict['CADD']="#00bd3c"
# color_dict['CADD']="#f58a67"


color_dict['RF:all']="#E64B35"
# color_dict['RF:TFBS']="#4DBBD5"
# color_dict['RF:Epi']="#00A087"
color_dict['DeepSEA']="#3C5488"
color_dict['CADD']="#F39B7F"
files = glob.glob("*auROC*iter*.csv")
sns.set(style="whitegrid")
def parse_csv(x):
	df = pd.read_csv(x)
	df = pd.melt(df)
	df['resolution'] = "$\pm$"+x.split("_")[2]+"bp"
	return df

df_list = [parse_csv(x) for x in files]
df = pd.concat(df_list)
# plt.figure(figsize=(12,5))
# ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["$\pm$0bp","$\pm$5bp","$\pm$50bp","$\pm$100bp","$\pm$500bp"],hue_order=['RF:all','RF:TFBS','RF:Epi','CADD','DeepSEA'])
# ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["$\pm$0bp","$\pm$5bp","$\pm$50bp","$\pm$100bp","$\pm$500bp"],hue_order=['RF:all','CADD','DeepSEA'])
# ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.mean,data=df,order=["$\pm$0bp","$\pm$5bp","$\pm$50bp","$\pm$100bp","$\pm$500bp"][::-1],hue_order=['RF:all','DeepSEA','CADD'],capsize=0.1)
print (df.shape)
f, (ax2, ax1) = plt.subplots(ncols=1, nrows=2, sharex=True,gridspec_kw={'height_ratios': [3.5, 1]})
# ax = sns.tsplot(time=x, data=y, ax=ax1)
# ax = sns.tsplot(time=x, data=y, ax=ax2)

# ax1.set_xlim(0, 6.5)
# ax2.set_xlim(13.5, 20)

ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,data=df,order=["$\pm$0bp","$\pm$5bp","$\pm$50bp","$\pm$100bp","$\pm$500bp"][::-1],hue_order=['RF:all','DeepSEA','CADD'],capsize=0.1,ax=ax1)
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,data=df,order=["$\pm$0bp","$\pm$5bp","$\pm$50bp","$\pm$100bp","$\pm$500bp"][::-1],hue_order=['RF:all','DeepSEA','CADD'],capsize=0.1,ax=ax2)
ax1.set_ylim(0, 0.1)
ax2.set_ylim(0.4, 0.75)



# hide the spines between ax and ax2
ax2.spines['bottom'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.xaxis.tick_top()
ax2.tick_params(labeltop='off')  # don't put tick labels at the top
ax1.xaxis.tick_bottom()

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax1.transAxes)  # switch to the bottom axes
ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


# plt.legend(loc='upper right',title="")
ax2.legend(loc='upper right',title="")
ax1.get_legend().remove()
ax2.set_xlabel("")
# ax1.set_xticks([0,0.1])
ax1.yaxis.set_ticks(np.arange(0, 0.1001, 0.1))
ax1.set_ylabel("")
ax2.set_ylabel("AUROC")
# plt.ylim(0,0.72)
plt.savefig("barplot_auROC_mean_only_all_NG.pdf", bbox_inches='tight')