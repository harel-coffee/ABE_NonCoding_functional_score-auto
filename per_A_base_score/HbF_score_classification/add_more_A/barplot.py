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
color_dict['RF:all']="#fa2111"
color_dict['RF:TFBS']="#fa8911"
color_dict['RF:Epi']="#fc42a5"
color_dict['DeepSEA']="#213fff"
color_dict['CADD']="#00bd3c"
files = glob.glob("*auROC*iter*.csv")

def parse_csv(x):
	df = pd.read_csv(x)
	df = pd.melt(df)
	df['resolution'] = "±"+x.split("_")[2]+"bp"
	return df

df_list = [parse_csv(x) for x in files]
df = pd.concat(df_list)
plt.figure(figsize=(12,5))
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["±0bp","±5bp","±50bp","±100bp","±500bp"],hue_order=['RF:all','RF:TFBS','RF:Epi','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,0.75)
plt.savefig("barplot_auROC_mean.pdf", bbox_inches='tight')