import pandas as pd
df = pd.read_csv("Editable_A_scores.combined.scores.csv")
df.head()
df[0] = [x[:-1].split(":")[0] for x in df.coord]
df[1] = [x[:-1].split(":")[1].split("-")[0] for x in df.coord]
df[2] = [x[:-1].split(":")[1].split("-")[1] for x in df.coord]
df.head()
def adjust_for_zero(x):
	if x == 0:
		return -0.1
	else:
		return x
df['HbFBase'] = [adjust_for_zero(x) for x in df['HbFBase']]
for j in ['HbFBase','DeepSEA','CADD']:
	df[[0,1,2,j]].to_csv("%s.bdg"%(j),index=False,header=False,sep="\t")
