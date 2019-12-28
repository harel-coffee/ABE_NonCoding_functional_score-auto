import pandas as pd
import numpy as np
from EmpiricalBrownsMethod import *
from scipy.stats import pearsonr
df = pd.read_csv("gRNA_all_A.bed",sep="\t",header=None)
print (df.head())
df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)+df[5].astype(str)
A_by_gRNA = pd.DataFrame(df.groupby('name')[3].agg(', '.join))
A_by_gRNA.columns = ['gRNAs']
abe = pd.read_csv("ABE_high_vs_low_mageck_RRA_results.sgrna_summary.txt",sep="\t")
abe.index = [x.split("-")[-1] for x in abe.sgrna]
myFDR = abe['FDR'].to_dict()
df[3] = df[3].map(myFDR)
df[3] = df[3].astype(str)
df[4] = df[4].astype(str)
A_by_FDR = pd.DataFrame(df.groupby('name')[3].agg(', '.join))
A_by_gRNA['FDRs'] = A_by_FDR[3]
A_by_pos = pd.DataFrame(df.groupby('name')[4].agg(', '.join))
A_by_gRNA['pos'] = A_by_pos[4]
A_by_gRNA.head()
A_by_gRNA['coord'] = A_by_gRNA.index.tolist()
A_by_gRNA.index = [x[:-1] for x in A_by_gRNA.index]
print (A_by_gRNA.head())

def row_apply(r):
    gRNA_list = r["gRNAs"].split(", ")
    FDR_list = r['FDRs'].split(", ")
    pos_list = r['pos'].split(", ")
    new_gRNA_list = []
    new_FDR_list = []
    new_pos_list = []
    for i in range(len(pos_list)):
        if 0<=int(pos_list[i])<=12: ## from 0 to 12 positions are observed editing events
            new_gRNA_list.append(gRNA_list[i])
            new_FDR_list.append(FDR_list[i])            
            new_pos_list.append(pos_list[i])   
    return [", ".join(new_gRNA_list),", ".join(new_FDR_list),", ".join(new_pos_list),len(new_gRNA_list)]

A_by_gRNA[['filtered_gRNAs','filtered_FDRs',"filtered_pos","length"]] = A_by_gRNA.apply(row_apply, axis=1).apply(pd.Series)

A_by_gRNA = A_by_gRNA[A_by_gRNA['length']>0]

freq = pd.read_csv("editing_frequency.list",header=None)
freq = freq[0].tolist()


def new_FDR(x):
    FDRs= [float(y) for y in x['filtered_FDRs'].split(", ")]
    new_FDRs = []
    pos_list = [int(y) for y in x['filtered_pos'].split(", ")]
    for i in range(len(pos_list)):
        pos = pos_list[i]
        FDR = FDRs[i]+1e-300
        new_FDRs.append(min(FDR/freq[pos],1))
    return ", ".join([str(f) for f in new_FDRs])    
A_by_gRNA['new_FDR'] = A_by_gRNA.apply(new_FDR,axis=1)
ABE_count = pd.read_csv("ABE_RRA_results.normalized.txt",sep="\t")
ABE_count.index = [x.split("-")[-1] for x in ABE_count.sgRNA]
ABE_count = ABE_count.drop(['sgRNA','Gene'],axis=1)
ABE_count.columns = [0,1,2,3]
print (ABE_count.head())

# EmpiricalBrownsMethod(DataMatrix, Pvalues, extra_info = True)
def combine_p_value(x):
    """row apply for A_by_gRNA"""
    ## get count matrix
    if len(x['filtered_gRNAs'].split(", "))==1:
        return float(x['new_FDR'])
    myCount = get_count_matrix(x['filtered_gRNAs'].split(", ")).values
    ## get p-value vector
    p_value_vector = np.array([float(y) for y in x['new_FDR'].split(", ")])
    result = EmpiricalBrownsMethod(myCount, p_value_vector, extra_info = True)
    return result[0]
    
def get_count_matrix(gRNA_list):
    return ABE_count.loc[gRNA_list]
A_by_gRNA['EBM_FDR']=A_by_gRNA.apply(combine_p_value,axis=1)    
A_by_gRNA['HbFBase'] = [-np.log10(x) for x in A_by_gRNA['EBM_FDR']]
A_by_gRNA.to_csv("Editable_A_scores.tsv",sep="\t")



