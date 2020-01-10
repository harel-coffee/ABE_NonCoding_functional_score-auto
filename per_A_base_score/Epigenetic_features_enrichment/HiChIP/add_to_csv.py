import pandas as pd
import glob
for i in ['pos','neg']:
	df1 = pd.read_csv(glob.glob("%s*BW.csv"%(i))[0],index_col=0)
	df2 = pd.read_csv(glob.glob("%s*degree.csv"%(i))[0],index_col=0)
	print (df1.head())
	
	df2.index = df2["0"]+":"+df2["1"].astype(str)+"-"+df2["2"].astype(str)
	print (df2.head())
	df1['#HiChIP_degree'] = df2["6"]
	df1.to_csv("%s_all_bw_plus_HiChIP.csv"%(i))
