import pandas as pd
df = pd.read_csv("locus-level-FDR.csv")
df.head()
df[0] = [x.split(":")[0] for x in df.l]
df[1] = [x.split(":")[1].split("-")[0] for x in df.l]
df[2] = [x.split(":")[1].split("-")[1] for x in df.l]
df.head()
for j in ['log_EBM_FDR']:
	df[[0,1,2,j]].to_csv("%s.bdg"%(j),index=False,header=False,sep="\t")
