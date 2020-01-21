#! /home/yli11/.conda/envs/py2/bin/python
"""
use cross validation to plot mean ROC curve, show std

ref:

https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html#sphx-glr-auto-examples-model-selection-plot-roc-crossval-py

Note that you have to tune the parameters yourself

"""
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from scipy import interp
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
# import xgboost as xgb
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
from sklearn.model_selection import KFold,StratifiedKFold
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression,RidgeClassifier,SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB 
from sklearn.ensemble import RandomForestClassifier
# from mlxtend.classifier import StackingCVClassifier
import os
import warnings
from sklearn.metrics import roc_curve,roc_auc_score,average_precision_score,precision_recall_curve
from sklearn import metrics
from sklearn.datasets import load_iris
# from mlxtend.classifier import StackingCVClassifier
# from mlxtend.feature_selection import ColumnSelector
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')
from sklearn.exceptions import ConvergenceWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=ConvergenceWarning)
from sklearn.ensemble import RandomForestRegressor,GradientBoostingRegressor,RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics.scorer import make_scorer
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
import uuid
from sklearn.base import TransformerMixin
from sklearn.datasets import make_regression
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor,GradientBoostingClassifier
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import LinearRegression, Ridge
import scipy
import numpy as np
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR,LinearSVC
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge,Lars,BayesianRidge
from copy import deepcopy as dp
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier,RadiusNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.gaussian_process import GaussianProcessClassifier
# from xgboost import XGBClassifier
from scipy import interp
from joblib import Parallel, delayed
def read_fasta(f):
	my_dict = {}
	for r in SeqIO.parse(f, "fasta"):
		my_dict[r.id] = str(r.seq).upper()
	return my_dict  
	

def sklearn_RF(par=False):

	est = RandomForestClassifier(n_estimators=1000,random_state=0,warm_start=False,n_jobs=-1)
	if par:
		est = RandomForestClassifier(**par)
	myDict = {}
	return est, myDict


def plot_top_features(reg,X,y,output):
	
	current_feature_df = pd.DataFrame()
	current_feature_df['features'] = X.columns.tolist()

	reg.fit(X,y)
	try:
		current_feature_df['score'] = list(reg.feature_importances_) 
	except:
		try:
			current_feature_df['score'] = list(reg.coef_) 
		except:
			current_feature_df['score'] = list(reg.coef_[0]) 
				
	current_feature_df = current_feature_df.sort_values('score',ascending=False)
	

	plt.figure(figsize=(len(current_feature_df['features']*2),8))
	sns.barplot(x=current_feature_df['features'],y=current_feature_df['score'] )
	plt.xticks(rotation=90)
	plt.xlabel("")
	plt.ylabel("Feature importance")
	plt.savefig("%s_feature_importance.pdf"%(output), bbox_inches='tight')  

def get_top_features(reg,X,y,output):
	
	current_feature_df = pd.DataFrame()
	current_feature_df['features'] = X.columns.tolist()

	reg.fit(X,y)
	try:
		current_feature_df['score'] = list(reg.feature_importances_) 
	except:
		try:
			current_feature_df['score'] = list(reg.coef_) 
		except:
			current_feature_df['score'] = list(reg.coef_[0]) 
				
	current_feature_df = current_feature_df.sort_values('score',ascending=False)
	return current_feature_df['features'].tolist()



def to_true_label(pos,neg):
	pos_fasta = read_fasta(pos)
	neg_fasta = read_fasta(neg)
	df = pd.DataFrame()
	df['true'] = [1]*len(pos_fasta.keys())+[-1]*len(neg_fasta.keys())
	df.index = pos_fasta.keys()+neg_fasta.keys()
	return df

	
def gkm_SVM_fit_transform(pos_train,neg_train,pos_test,neg_test,d=3):
	"""From HemTools"""
	
	addon_string = str(uuid.uuid4()).split("-")[-1]
	train_out = addon_string
	classify_out = addon_string + ".classify.out"
	train_command = "gkmtrain %s %s %s >/dev/null 2>&1"%(pos_train,neg_train, train_out)
	
	classify_input = addon_string + ".classify_input.fa"
	combine_pos_neg_command = "cat %s %s > %s"%(pos_test,neg_test,classify_input)
	classify_command = "gkmpredict %s %s.model.txt %s >/dev/null 2>&1"%(classify_input,train_out,classify_out)
	## run
	os.system(train_command)
	os.system(combine_pos_neg_command)
	os.system(classify_command)
	
	y_pred = pd.read_csv(classify_out,sep="\t",header=None,index_col=0)
	y_pred.columns = ['pred']
	os.system("rm %s*"%(addon_string))  
	return y_pred
	
def write_fasta(file_name,myDict):
	out = open(file_name,"wt")
	for k in myDict:
		out.write(">"+k+"\n")
		out.write(myDict[k]+"\n")
	out.close()
def get_fasta_file(X_train,y_train,X_test,y_test):
	addon_string = str(uuid.uuid4()).split("-")[-1]
	pos_train_file = addon_string+"_pos_train_file.fa"
	neg_train_file = addon_string+"_neg_train_file.fa"
	pos_test_file = addon_string+"_pos_test_file.fa"
	neg_test_file = addon_string+"_neg_test_file.fa"


	write_fasta(pos_train_file,X_train.loc[y_train[y_train==1].index.tolist()].to_dict())
	write_fasta(pos_test_file,X_test.loc[y_test[y_test==1].index.tolist()].to_dict())
	write_fasta(neg_train_file,X_train.loc[y_train[y_train==0].index.tolist()].to_dict())
	write_fasta(neg_test_file,X_test.loc[y_test[y_test==0].index.tolist()].to_dict())

	return pos_train_file,neg_train_file,pos_test_file,neg_test_file,addon_string

def find_neighbor(x,y,distance=50):
	chr_x = x.split(":")[0]
	start_x = int(x.split(":")[-1].split("-")[0]) 
	chr_y = y.split(":")[0]
	start_y = int(y.split(":")[-1].split("-")[0]) 
	# print (chr_x,chr_y,start_x,start_y)
	if chr_x != chr_y:
		return False
	if abs(start_y-start_x)<=distance:
		return True
	return False
	

def leave_one_chr_out_CV(y):
	"""Y is a list of genomic coord"""
	tmp = pd.DataFrame(y)
	tmp[1] = [x.split(":")[0] for x in y]
	# print (tmp.head())
	out = []
	for s,d in tmp.groupby(1):
		out.append([tmp[tmp[1]!=s][0].tolist(),tmp[tmp[1]==s][0].tolist()])
	return out

	
def random_sample_balance(y,random_state):
	pos_size = y[y==1].shape[0]
	# print ("pos_size",pos_size)
	sub_sample_index_list = y[y==1].index.tolist()+y[y==0].sample(n=pos_size,random_state=random_state).index.tolist()

	return leave_one_chr_out_CV(y.loc[sub_sample_index_list].index.tolist())

def piecewise_roc(x,y):

	plot_df = pd.DataFrame()
	plot_df['x'] = list(x)
	plot_df['y'] = list(y)
	plot_df = plot_df.sort_values('x')
	mean_x = np.linspace(0, 1, 1000)
	mean_y = interp(mean_x, plot_df['x'], plot_df['y'])
	mean_y[0] = 0.0
	mean_y[-1] = 1
	return mean_x,mean_y

def piecewise_prc(x,y):

	plot_df = pd.DataFrame()
	plot_df['x'] = list(x)
	plot_df['y'] = list(y)
	plot_df = plot_df.sort_values('x')
	mean_x = np.linspace(0, 1, 1000)
	mean_y = interp(mean_x, plot_df['x'], plot_df['y'])
	mean_y[0] = 1
	mean_y[-1] = 0.5
	return mean_x,mean_y


def simple_CV_evaluation(model,X,y,k,random_state):
	
	my_pred=[]
	my_true=[]

	auPRC_list = []
	auROC_list = []
	auROC_x_list = []
	auROC_y_list = []
	auPRC_x_list = []
	auPRC_y_list = []
	for train_index, test_index in random_sample_balance(y,random_state):
		# print ("train size: %s test size: %s"%(len(train_index),len(test_index)))
		# print ("train_index:",np.random.choice(train_index,3))
		# print ("test_index:",np.random.choice(test_index,3))
		# print ("total",len(train_index)+len(test_index))
		X_train, X_test = X.loc[train_index], X.loc[test_index]
		y_train, y_test = y.loc[train_index], y.loc[test_index]

		if "RF" in k:
			current_model = dp(model)
			current_model.fit(X_train,y_train)
			pred_y = current_model.predict_proba(X_test)
			pred_y = [x[1] for x in pred_y]
		elif k == "LS-GKM":

			pos_train_file,neg_train_file,pos_test_file,neg_test_file,addon_string = get_fasta_file(X_train,y_train,X_test,y_test)

			pred_y = gkm_SVM_fit_transform(pos_train_file,neg_train_file,pos_test_file,neg_test_file)
			pred_y = pred_y.loc[y_test.index.tolist()]['pred'].tolist()
			os.system("rm %s*"%(addon_string))

		else:

			pred_y = X_test.tolist()

		y_test = y_test.tolist()
		my_pred += pred_y
		my_true += y_test
		try:
			auROC = roc_auc_score(y_test,pred_y)
			auPRC = average_precision_score(y_test,pred_y)
			# print ("model %s auPRC: %s. auROC: %s"%(k,auPRC,auROC))
			auPRC_list.append(auPRC)
			auROC_list.append(auROC)
			x_predict,y_predict,_ = roc_curve(y_test,pred_y)
			x_predict,y_predict = piecewise_roc(x_predict,y_predict)
			
			
			auROC_x_list+=list(x_predict)
			auROC_y_list+=list(y_predict)		
			y_predict,x_predict,_ = precision_recall_curve(y_test,pred_y)
			x_predict,y_predict = piecewise_prc(x_predict,y_predict)
			auPRC_x_list+=list(x_predict)
			auPRC_y_list+=list(y_predict)			
		except:
			print ("y_test pos :",y_test[y_test==1].shape)
			print ("y_test neg :",y_test[y_test==0].shape)		
			print ("y_train pos :",y_train[y_train==1].shape)
			print ("y_train neg :",y_train[y_train==0].shape)		
		
			pass

	df = pd.DataFrame()
	
	df['true']=my_true
	df['pred']=my_pred
	df['label'] = k
	return df,auROC_list,auPRC_list,[auROC_x_list,auROC_y_list],[auPRC_x_list,auPRC_y_list]

def plot_auROC_multi(df,color_dict,output):
	sns.set_style("white")
	plt.figure()
	
	
	for s,d in df.groupby('label'):
		plot_df = pd.DataFrame()
		x_predict,y_predict,_ = roc_curve(d['true'],d['pred'])
		auc = roc_auc_score(d['true'],d['pred'])
		# print (auc)
		plot_df['x'] = x_predict
		plot_df['y'] = y_predict
		plot_df.to_csv("plot_df_auROC_%s_x_y_resolution_%s.csv"%(s,output),index=False)
		sns.lineplot(data=plot_df,x="x",y="y",ci=0,label="%s AUC:%.2f"%(s,auc),color=color_dict[s])
	plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('False positive rate')
	plt.ylabel('True positive rate')	
	plt.title('ROC curve')
	plt.legend(loc='lower right',title="")
	handles, labels = plt.gca().get_legend_handles_labels()
	# order = [5,4,3,2,0,1]
	order = [4,3,2,0,1]
	plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])	

	plt.savefig("%s_auROC.pdf"%(output), bbox_inches='tight')
	plt.close()

def plot_average_auROC_multi(x,y,roc_mean,color_dict,output):
	sns.set_style("white")
	plt.figure()
	for k in x:
		plot_df = pd.DataFrame()
		plot_df['x'] = x[k]
		plot_df['y'] = y[k]
		plot_df = plot_df.sort_values('x')
		# print (plot_df.head())
		auc = roc_mean[k]
		print (k,auc)
		# mean_fpr = np.linspace(0, 1, 25)
		# interp_tpr = interp(mean_fpr, plot_df['x'], plot_df['y'])
		# interp_tpr[0] = 0.0
		# interp_tpr[-1] = 1
		# plot_df = pd.DataFrame()
		# plot_df['x'] = mean_fpr
		# plot_df['y'] = interp_tpr
		print (plot_df.head())
		sns.lineplot(data=plot_df,x="x",y="y",ci=0,label="%s AUC:%.2f"%(k,auc),color=color_dict[k])
	plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('False positive rate')
	plt.ylabel('True positive rate')	
	plt.title('ROC curve')
	plt.legend(loc='lower right',title="")
	# handles, labels = plt.gca().get_legend_handles_labels()
	# order = [5,4,3,2,0,1]
	# order = [4,3,2,0,1]
	# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
	plt.savefig("%s_auROC.pdf"%(output), bbox_inches='tight')
	plt.close()


def plot_auPRC_multi(df,color_dict,output):
	sns.set_style("white")
	plt.figure()

	for s,d in df.groupby('label'):
		plot_df = pd.DataFrame()
		y_predict,x_predict,_ = precision_recall_curve(d['true'],d['pred'])
		auc = average_precision_score(d['true'],d['pred'])
		# print (auc)
		plot_df['x'] = x_predict
		plot_df['y'] = y_predict
		plot_df.to_csv("plot_df_auPRC_%s_x_y_resolution_%s.csv"%(s,output),index=False)
		sns.lineplot(data=plot_df,x="x",y="y",ci=0,label="%s AUC:%.2f"%(s,auc),color=color_dict[s])
		# plt.step(x_predict, y_predict, where='post',label="%s AUC:%.2f"%(s,auc),color=color_dict[s])
	# plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('Recall')
	plt.ylabel('Precision') 
	plt.title('Precision-Recall curve')
	plt.legend(loc='best',title="")
	# plt.savefig("auPRC.png")
	plt.savefig("%s_auPRC.pdf"%(output), bbox_inches='tight')
	plt.close()

def plot_average_auPRC_multi(x,y,prc_mean,color_dict,output):
	sns.set_style("white")
	plt.figure()
	for k in x:
		plot_df = pd.DataFrame()
		auc = prc_mean[k]
		# mean_recall = np.linspace(0, 1, 1000)
		# interp_prec = interp(mean_recall, x[k], y[k])
		# interp_prec[0] = 1
		# interp_prec[-1] = 0.5
		# plot_df['x'] = mean_recall
		# plot_df['y'] = interp_prec			
		# print (k,auc)
		plot_df['x'] = x[k]
		plot_df['y'] = y[k]
		print (plot_df.head())
		sns.lineplot(data=plot_df,x="x",y="y",ci=0,label="%s AUC:%.2f"%(k,auc),color=color_dict[k])
	# plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('Recall')
	plt.ylabel('Precision') 
	plt.title('Precision-Recall curve')
	plt.legend(loc='lower right',title="")
	# handles, labels = plt.gca().get_legend_handles_labels()
	# order = [5,4,3,2,0,1]
	# order = [4,3,2,0,1]
	# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])	
	plt.savefig("%s_auPRC.pdf"%(output), bbox_inches='tight')
	plt.close()


def define_high_low(x,cutoff):
	if x>=cutoff:
		return 1
	elif x<=10:
		return 0
	else:
		return -1

def define_low_filter(x,pos_list,distance):
	"""keep 1 remove 0"""
	if x in pos_list:
		return 1
	for p in pos_list:
		if find_neighbor(x,p,distance):
			return 0
	return 1


def boxplot_paired_t_test(df,color_dict,ylabel,output):
	sns.set_style("whitegrid")
	myMin = df.min().min()
	myMax = df.max().max()
	plot_df = pd.melt(df)
	plt.figure()
	ax=sns.boxplot(x="variable",y='value',data=plot_df,palette =color_dict,linewidth=3,order=['RF-motif-epi','LS-GKM','CADD','DeepSEA'])
	for patch in ax.artists:
		r, g, bb, _ = patch.get_facecolor()
		patch.set_facecolor((r, g, bb, .3))
	
	sns.swarmplot(x="variable", y="value", data=plot_df,palette =color_dict,order=['RF-motif-epi','LS-GKM','CADD','DeepSEA'])
	unit=0.01
	
	y=myMax+unit*1
	h=unit*1.5
	plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c="black")
	plt.text(0.5, y+h+unit, "Paired T-test: %.2E" % scipy.stats.ttest_rel(df['RF-motif-epi'],df['LS-GKM']).pvalue, ha='center', va='bottom', color="black")
	
	# y=myMax+unit*6
	# h=unit*1.5
	# plt.plot([0, 0, 2, 2], [y, y+h, y+h, y], lw=1.5, c="black")
	# plt.text(0.5, y+h+unit, "Paired T-test: %.2E" % scipy.stats.ttest_rel(df['RF-motif-epi'],df['CADD']).pvalue, ha='center', va='bottom', color="black")
	print ("RF vs. CADD Paired T-test: %.2E" % scipy.stats.ttest_rel(df['RF-motif-epi'],df['CADD']).pvalue)
	print ("RF vs. DeepSEA Paired T-test: %.2E" % scipy.stats.ttest_rel(df['RF-motif-epi'],df['DeepSEA']).pvalue)
	# y=myMax+unit*11
	# h=unit*1.5
	# plt.plot([0, 0, 3, 3], [y, y+h, y+h, y], lw=1.5, c="black")
	# plt.text(0.5, y+h+unit, "Paired T-test: %.2E" % scipy.stats.ttest_rel(df['RF-motif-epi'],df['DeepSEA']).pvalue, ha='center', va='bottom', color="black")
		
	plt.ylim(myMin-unit*5,myMax+0.1)
	plt.ylabel(ylabel)
	plt.savefig("%s.pdf"%(output), bbox_inches='tight')


#--------------- DEFINE PARAMETERS -----------------#

N_samples = 100
HbF_cutoff = 30
distance = int(sys.argv[1])
# distance = 500


## define colors for line plot and boxplot
color_dict = {}
# color_dict['RF:all']="#fa2111"
# color_dict['RF:TFBS']="#fa8911"
# color_dict['RF:Epi']="#fc42a5"
# color_dict['DeepSEA']="#213fff"
# color_dict['CADD']="#00bd3c"
# https://nanx.me/ggsci/reference/pal_npg.html
color_dict['RF:all']="#E64B35"
color_dict['RF:TFBS']="#4DBBD5"
color_dict['RF:Epi']="#00A087"
color_dict['DeepSEA']="#3C5488"
color_dict['CADD']="#F39B7F"
# color_dict['LS-GKM']="#464d4f"

auROC_dict={}
auPRC_dict={}
roc_x_dict = {}
roc_y_dict = {}
prc_x_dict = {}
prc_y_dict = {}
for k in color_dict:
	auROC_dict[k]=[]
	auPRC_dict[k]=[]
	roc_x_dict[k]=[]
	roc_y_dict[k]=[]
	prc_x_dict[k]=[]
	prc_y_dict[k]=[]


## step 1: get X,Y for each of the 6 types


## define global X Y
df = pd.read_csv("ML_data_complete.csv",index_col=0)
target = "HbFBase"
df['label'] = [define_high_low(x,HbF_cutoff) for x in df[target]]
Y = df['label']
print ("distance",distance)
print ("Negative Size:",Y[Y==0].shape)
print ("Positive Size:",Y[Y==1].shape)

df['filter'] = [define_low_filter(x,Y[Y!=0].index.tolist(),distance) for x in df.index.tolist()]
# df.to_csv("filtered_label_ML_data.csv")
df = df[df['filter']==1]
df = df[df['label']!=-1]
Y = df['label']
print ("After removing positive A neighbors")
print ("Negative Size:",Y[Y==0].shape)
print ("Positive Size:",Y[Y==1].shape)
X = df.drop([target,'label','filter'],axis=1)

## filtering step didn't change HbF score distribution
for s,d in df.groupby("label"):
	print (s)
	print (d[target].describe())

## get X,Y for each of the 6 types

X_DeepSEA = X['DeepSEA']
X_CADD = X['CADD']
X_gkm = X['seq']
not_all_features = ['DeepSEA','CADD','seq']
TFBS_features = []
epi_features = []
for c in X.columns:
	if "motif" in c:
		TFBS_features.append(c)
	else:
		if c in not_all_features:
			continue
		epi_features.append(c)
X_TFBS = X[TFBS_features]
X_epi = X[epi_features]
X_all = X.drop(not_all_features,axis=1)

print ("all features size:",X_all.shape)
print ("epi features size:",X_epi.shape)
print ("TFBS features size:",X_TFBS.shape)


## RF model
model,params=sklearn_RF()



X_dict = {}
X_dict['RF:all']=X_all
X_dict['RF:TFBS']=X_TFBS
X_dict['RF:Epi']=X_epi
X_dict['DeepSEA']=X_DeepSEA
X_dict['CADD']=X_CADD
# X_dict['LS-GKM']=X_gkm 


## RF top features
plot_top_features(dp(model),X_all,Y,"All_features_distance_%s"%(distance))




## step 2: random sample X and Y, and at least 100 samples



def parallel_unit(model,X_dict,Y):
	random_state = np.random.randint(0,999999)
	out_dict = {}
	for k in X_dict:
		ddf,roc_list,prc_list,roc_x_y,prc_x_y = simple_CV_evaluation(dp(model),X_dict[k],Y,k,random_state)   
		out_dict[k]=[ddf,roc_list,prc_list,roc_x_y,prc_x_y]
	return out_dict
out_list = Parallel(n_jobs=-1)(delayed(parallel_unit)(model,X_dict,Y) for i in range(N_samples))




## summarize results
## out_list, from joblib classification
df_list = []

for item in out_list:
	for k in X_dict:
		auROC_dict[k]+=item[k][1]
		auPRC_dict[k]+=item[k][2]
		roc_x_dict[k]+=item[k][3][0]
		roc_y_dict[k]+=item[k][3][1]
		prc_x_dict[k]+=item[k][4][0]
		prc_y_dict[k]+=item[k][4][1]
		df_list.append(item[k][0])
plot_df = pd.concat(df_list)
roc_mean = pd.DataFrame(auROC_dict).mean().to_dict()
# print (roc_mean)
prc_mean = pd.DataFrame(auPRC_dict).mean().to_dict()
plot_average_auROC_multi(roc_x_dict,roc_y_dict,roc_mean,color_dict,"average_auROC_%s"%(distance))
plot_average_auPRC_multi(prc_x_dict,prc_y_dict,prc_mean,color_dict,"average_auPRC_%s"%(distance))

# plot_auROC_multi(plot_df,color_dict,distance)
# plot_auPRC_multi(plot_df,color_dict,distance)
plot_df.to_csv("plot_df_%s_100iter.csv"%(distance))
pd.DataFrame(auROC_dict).to_csv("auROC_dict_%s_100iter.csv"%(distance),index=False)
pd.DataFrame(auPRC_dict).to_csv("auPRC_dict_%s_100iter.csv"%(distance),index=False)


