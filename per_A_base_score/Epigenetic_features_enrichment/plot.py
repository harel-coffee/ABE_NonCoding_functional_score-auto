
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
from decimal import Decimal

pos = pd.read_csv("pos.bed._yli11_2020-01-06_84353d2d3842.averageBW.csv",index_col=0)
neg = pd.read_csv("neg.bed._yli11_2020-01-06_84353d2d3842.averageBW.csv",index_col=0)

def plot_violin(pos,neg,feature_name):
	
	sns.set_style("whitegrid")
	plt.figure()
	top_n = [np.log2(x+1) for x in pos]
	bot_n = [np.log2(x+1) for x in neg]
	plot_df = pd.DataFrame([top_n,bot_n]).T
	plot_df.columns = ['High',"Low"]
	print (plot_df.describe())
	plot_df = pd.melt(plot_df)
	color_dict={}
	color_dict['High'] = "red"
	color_dict['Low'] = "blue"
	sns.violinplot(x="variable",y='value',data=plot_df,palette =color_dict,linewidth=3,width=0.7,cut=3)


	y=plot_df['value'].max()+1
	h=0.5
	# print (scipy.stats.mannwhitneyu(pos,neg).pvalue)
	plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c="black")
	plt.text(0.5, y+h+0.05, "Mann-Whitney U test: %.2E" % scipy.stats.mannwhitneyu(top_n,bot_n).pvalue, ha='center', va='bottom', color="black")
	plt.ylim(-0.5,plot_df['value'].max()+2)
	plt.xticks([0,1],['High HbF score','Low HbF score'])
	plt.xlabel("HbFBase scores")
	plt.ylabel(feature_name)
	plt.savefig("%s-HbFBase-high-low.pdf"%(feature_name), bbox_inches='tight')

for c in pos.columns:
	plot_violin(pos[c].tolist(),neg[c].tolist(),c)

















