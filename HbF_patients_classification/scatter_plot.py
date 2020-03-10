import glob
import matplotlib
import pandas as pd
matplotlib.use('agg')
import seaborn as sns
import numpy as np
import scipy
import glob
import sys
import matplotlib.pyplot as plt
import os
import logging
df = pd.read_csv("snp_data.tsv",sep="\t")
df.head()
common_list = []
novel_list = []
for c in df.columns:
	if c == "HbF":
		continue
	if "chr" in c:
		novel_list.append(c)
	else:
		common_list.append(c)
common_list
df[common_list].sum()
df['common'] = df[common_list].sum(axis=1)
df['novel'] = df[novel_list].sum(axis=1)
df.head()
# plt.figure()
# sns.scatterplot(data=df,x="common",y="HbF",y_bins)
# plt.savefig('common.pdf', bbox_inches='tight')

from decimal import Decimal

def plot_correlation(df):

	plt.figure()
	sns.regplot(x=df['common'],y=df['HbF'],x_bins =20)
	r,p = scipy.stats.pearsonr(df["common"],df["HbF"])
	plt.xlabel("number of common variants")
	plt.ylabel("HbF")
	plt.text(df['common'].quantile(.05), df['HbF'].quantile(.85), 'r=%f'%(r))
	plt.text(df['common'].quantile(.05), df['HbF'].quantile(.8), 'p=%.2E'%(Decimal(p)))

	plt.savefig("common.pdf", bbox_inches='tight')
plot_correlation(df)
plt.close()
def plot_correlation2(df):

	plt.figure()
	sns.regplot(x=df['novel'],y=df['HbF'],x_bins =20)
	r,p = scipy.stats.pearsonr(df["novel"],df["HbF"])
	plt.xlabel("number of novel variants")
	plt.ylabel("HbF")
	plt.text(df['novel'].quantile(.05), 14, 'r=%f'%(r))
	plt.text(df['novel'].quantile(.05), 13.5, 'p=%.2E'%(Decimal(p)))

	plt.savefig("novel.pdf", bbox_inches='tight')
plot_correlation2(df)
plt.close()
df['all'] = df['novel']+df['common']
def plot_correlation3(df):

	plt.figure()
	sns.regplot(x=df['all'],y=df['HbF'],x_bins =20)
	r,p = scipy.stats.pearsonr(df["all"],df["HbF"])
	plt.xlabel("number of novel variants")
	plt.ylabel("HbF")
	plt.text(df['all'].quantile(.05), df['HbF'].quantile(.85), 'r=%f'%(r))
	plt.text(df['all'].quantile(.05), df['HbF'].quantile(.8), 'p=%.2E'%(Decimal(p)))

	plt.savefig("all.pdf", bbox_inches='tight')
plot_correlation3(df)