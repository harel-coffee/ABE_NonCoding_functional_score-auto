
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

cadd = pd.read_csv("DeepSEA_CADD_GERP/gRNA_all_A.CADD.vcf",sep="\t")
deepsea = pd.read_csv("DeepSEA_CADD_GERP/DeepSEA.out.funsig",index_col=0)
gerp =pd.read_csv("DeepSEA_CADD_GERP/gRNA_all_GERP.tsv",sep="\t")


gerp.index = gerp['chrom']+":"+gerp['start'].astype(str)+"-"+gerp['end'].astype(str)
deepsea['name'] = deepsea['chr']+":"+(deepsea['pos']-1).astype(str)+"-"+deepsea['pos'].astype(str)
deepsea.index = deepsea['name']
cadd['name'] = "chr"+cadd['#Chrom'].astype(str)+":"+(cadd['Pos']-1).astype(str)+"-"+cadd['Pos'].astype(str)
cadd.index = cadd['name']


df = pd.read_csv("Editable_A_scores.tsv",sep="\t",index_col=0)
df['CADD'] = cadd['PHRED']
df['DeepSEA'] = deepsea['Functional significance score']
df['GERP'] = gerp['gerp_bp_score']
df['DeepSEA'] = df['DeepSEA'].apply(lambda x:-np.log10(x))
df.index = df.coord
df[['CADD','DeepSEA','GERP','HbFBase']].to_csv("Editable_A_scores.combined.scores.csv")
df = pd.read_csv("Editable_A_scores.combined.scores.csv",index_col=0)

# deepsea violin

from decimal import Decimal
sns.set_style("whitegrid")
plt.figure()
top_n = df[df['HbFBase']>=50]['DeepSEA'].tolist()
bot_n = df[df['HbFBase']==0]['DeepSEA'].tolist()
plot_df = pd.DataFrame([top_n,bot_n]).T
plot_df.columns = ['High',"Low"]
print (plot_df.describe())
plot_df = pd.melt(plot_df)
color_dict={}
color_dict['High'] = "#213fff"
color_dict['Low'] = "#6e899c"
sns.violinplot(x="variable",y='value',data=plot_df,palette =color_dict,linewidth=3,width=0.7,cut=3)

import matplotlib.pyplot as plt
y=5.2
h=0.3
print (scipy.stats.mannwhitneyu(top_n,bot_n).pvalue)
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c="black")
plt.text(0.5, y+h+0.05, "Mann-Whitney U test: %.2E" % scipy.stats.mannwhitneyu(top_n,bot_n).pvalue, ha='center', va='bottom', color="black")
plt.ylim(-0.5,6)
plt.xticks([0,1],['High HbF score','Low HbF score'])
plt.xlabel("HbFBase scores")
plt.ylabel("DeepSEA scores (-log10)")
plt.savefig("DeepSEA-HbFBase-high-low.pdf", bbox_inches='tight')


# CADD violin

from decimal import Decimal
sns.set_style("whitegrid")
plt.figure()
top_n = df[df['HbFBase']>=50]['CADD'].tolist()
bot_n = df[df['HbFBase']==0]['CADD'].tolist()
print (scipy.stats.mannwhitneyu(top_n,bot_n).pvalue)

plot_df = pd.DataFrame([top_n,bot_n]).T
plot_df.columns = ['High',"Low"]
# print (plot_df.describe())
plot_df = pd.melt(plot_df)
color_dict={}
color_dict['High'] = "#00bd3c"
color_dict['Low'] = "#7d827e"
# plt.figure(figsize=(7,4))
sns.violinplot(x="variable",y='value',data=plot_df,palette =color_dict,linewidth=3,width=0.7,cut=3)
import matplotlib.pyplot as plt
y=38
h=2
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c="black")
plt.text(0.5, y+h+0.05, "Mann-Whitney U test: %.2E" % scipy.stats.mannwhitneyu(top_n,bot_n).pvalue, ha='center', va='bottom', color="black")
plt.ylim(-10,45)
# plt.xlim(-1,2)
# my_list = []
# for k in color_dict:
# plt.legend(handles=my_list)
plt.xticks([0,1],['High HbF score','Low HbF score'])
plt.xlabel("HbFBase scores")
plt.ylabel("CADD scores")
plt.savefig("CADD-HbFBase-high-low.pdf", bbox_inches='tight')


# GERP violin

from decimal import Decimal
sns.set_style("whitegrid")
plt.figure()
top_n = df[df['HbFBase']>=50]['GERP'].tolist()
bot_n = df[df['HbFBase']==0]['GERP'].tolist()
plot_df = pd.DataFrame([top_n,bot_n]).T
plot_df.columns = ['High',"Low"]
print (plot_df.describe())
plot_df = pd.melt(plot_df)
color_dict={}
color_dict['High'] = "#b569e0"
color_dict['Low'] = "#4d2d5e"
sns.violinplot(x="variable",y='value',data=plot_df,palette =color_dict,linewidth=3,width=0.7,cut=0)
import matplotlib.pyplot as plt
y=7
h=0.3
# print (scipy.stats.ttest_ind(top_n,bot_n).pvalue)
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c="black")
plt.text(0.5, y+h+0.05, "Mann-Whitney U test: %.2E" % scipy.stats.mannwhitneyu(top_n,bot_n).pvalue, ha='center', va='bottom', color="black")
plt.ylim(-12,9)
plt.xticks([0,1],['High HbF score','Low HbF score'])
plt.xlabel("HbFBase scores")
plt.ylabel("GERP")
plt.savefig("GERP-HbFBase-high-low.pdf", bbox_inches='tight')



