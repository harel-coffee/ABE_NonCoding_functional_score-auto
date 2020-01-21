import pandas as pd
df = pd.read_csv("all_A_features.csv",index_col=0)
df2 = pd.read_csv("editable_A_all_bw_plus_HiChIP.csv",index_col=0)
df.index = [x[:-1] for x in df.index]
df = pd.concat([df,df2],axis=1)
print (df.isnull().any().any())
df = df.drop(['e53295538398-Hudep2_D0_GATA1.all','386b0c66251a-1686607_H2_ABEmax_RNAseq.all','deaf983df0e4-1788283_H2_ABEmax.all','GERP'],axis=1)
df.to_csv("ML_data_complete.csv")
# import readline
# readline.write_history_file("combine_data.py")
