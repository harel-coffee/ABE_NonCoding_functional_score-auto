
import pandas as pd
import networkx as nx

def to_bed(x,out):
	x[1] = x[1].astype(int)
	x[2] = x[2].astype(int)
	x.to_csv(out,sep="\t",header=False,index=False)

def read_bedpe(f):
	df = pd.read_csv(f,sep="\t",header=None)
	df[1] = df[1].astype(int)
	df[2] = df[2].astype(int)	
	df[4] = df[4].astype(int)
	df[5] = df[5].astype(int)	
	df = df.dropna()
	return df
	
def read_bed(f):
	df = pd.read_csv(f,sep="\t",header=None)
	df[1] = df[1].astype(int)
	df[2] = df[2].astype(int)	
	df = df.dropna()
	return df
	
def mango2bed(df):
	# print (df[[0,1,2]].head())
	# print (df[[3,4,5]].head())
	tmp = df[[0,1,2]].copy()
	tmp2 = df[[3,4,5]].copy()
	tmp2.columns = tmp.columns
	
	tmp = pd.concat([tmp,tmp2],axis=0)
	# print (tmp.head())
	tmp['name'] = tmp[0]+"-"+tmp[1].astype(str)+"-"+tmp[2].astype(str)
	tmp = tmp.drop_duplicates('name')
	to_bed(tmp,"mango.bed")
	
mango_file = "Hudep2_D0_H3K27AC_HiChIP_FS.interactions.all.mango"
import sys
input_query_bed = sys.argv[1]
df = read_bedpe(mango_file)
df['source'] = df[0]+"-"+df[1].astype(str)+"-"+df[2].astype(str)
df['target'] = df[3]+"-"+df[4].astype(str)+"-"+df[5].astype(str)
print (df.head())
g=nx.from_pandas_edgelist(df, 'source', 'target')
my_dict = g.degree
# print (my_dict)
df1 = [[x[0],x[1]] for x in g.degree]
df1 = pd.DataFrame(df1)
df1.columns=['Genomic coordinates','Number of contacts']
df1 = df1.set_index('Genomic coordinates')
print (df1.head())
mango2bed(df)

## overlap query bed with mango bed
import os
os.system("bedtools intersect -a %s -b mango.bed -wao > intersect.bed"%(input_query_bed))

df2 = read_bed("intersect.bed")
df2[df2.columns[-2]] = df2[df2.columns[-2]].map(df1['Number of contacts'].to_dict())
df2 = df2.fillna(0)
df2.to_csv("%s_degree.csv"%(input_query_bed))