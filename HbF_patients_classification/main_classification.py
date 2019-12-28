
"""
use cross validation to plot mean ROC curve, show std

ref:

https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html#sphx-glr-auto-examples-model-selection-plot-roc-crossval-py

Note that you have to tune the parameters yourself

"""
from scipy import interp
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import xgboost as xgb
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
import umap
import warnings
from sklearn.metrics import roc_curve,roc_auc_score,average_precision_score,precision_recall_curve
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


def sklearn_RF(par=False):

	est = RandomForestClassifier(n_estimators=1000,random_state=0,warm_start=False,n_jobs=-1,class_weight={1:4,0:1})
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
	

def simple_CV_evaluation(model,params,X,y):

	outer = StratifiedKFold(n_splits=3,shuffle=True)
	my_pred=[]
	my_true=[]

	best_features = X.columns.tolist()

	auPRC_list = []
	auROC_list = []
	for train_index, test_index in outer.split(X,y):
		X_train, X_test = X.iloc[train_index], X.iloc[test_index]
		y_train, y_test = y.iloc[train_index], y.iloc[test_index]
		
		model.fit(X_train[best_features].values,y_train)
		pred_y = model.predict_proba(X_test[best_features].values)
		pred_y = [x[1] for x in pred_y]
		y_test = y_test.tolist()
		auROC = roc_auc_score(y_test,pred_y)
		auPRC = average_precision_score(y_test,pred_y)
		my_pred += pred_y
		my_true += y_test
		print ("auPRC: %s. auROC: %s"%(auPRC,auROC))
		auPRC_list.append(auPRC)
		auROC_list.append(auROC)

	df = pd.DataFrame()
	df['true']=my_true
	df['pred']=my_pred
	return df,auROC_list,auPRC_list

def plot_auROC_multi(df,color_dict):
	sns.set_style("white")
	plt.figure()
	
	
	for s,d in df.groupby('label'):
		plot_df = pd.DataFrame()
		x_predict,y_predict,_ = roc_curve(d['true'],d['pred'])
		auc = roc_auc_score(d['true'],d['pred'])
		print (auc)
		plot_df['x'] = x_predict
		plot_df['y'] = y_predict
		sns.lineplot(data=plot_df,x="x",y="y",ci=0,label="%s AUC:%.2f"%(s,auc),color=color_dict[s])
	plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('False positive rate')
	plt.ylabel('True positive rate')	
	plt.title('ROC curve')
	plt.legend(loc='best',title="")
	# plt.savefig("auROC.png")
	plt.savefig("auROC.pdf", bbox_inches='tight')
	plt.close()

def plot_auPRC_multi(df,color_dict):
	sns.set_style("white")
	plt.figure()

	for s,d in df.groupby('label'):
		plot_df = pd.DataFrame()
		x_predict,y_predict,_ = precision_recall_curve(d['true'],d['pred'])
		auc = average_precision_score(d['true'],d['pred'])
		print (auc)
		plot_df['x'] = x_predict
		plot_df['y'] = y_predict
		sns.lineplot(data=plot_df,x="x",y="y",ci=0,label="%s AUC:%.2f"%(s,auc),color=color_dict[s])
	# plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.xlabel('Recall')
	plt.ylabel('Precision')	
	plt.title('Precision-Recall curve')
	plt.legend(loc='best',title="")
	# plt.savefig("auPRC.png")
	plt.savefig("auPRC.pdf", bbox_inches='tight')
	plt.close()

def define_high_low(x,mu,sigma):
	t = 1
	low = mu-t*sigma
	high = mu+t*sigma
	high2= mu+t*sigma
	# print (low,high)
	if low <= x <= high:
		return 0
	if x > high2:
		return 1
	return -1

def boxplot_paired_t_test(a,b,color_dict,ylabel,output):
	sns.set_style("whitegrid")
	df = pd.DataFrame()
	df['All_variants'] = a
	df['GWAS_only'] = b
	myMin = df.min().min()
	myMax = df.max().max()
	plot_df = pd.melt(df)
	plt.figure()
	ax=sns.boxplot(x="variable",y='value',data=plot_df,palette =color_dict,linewidth=3)
	for patch in ax.artists:
		r, g, bb, _ = patch.get_facecolor()
		patch.set_facecolor((r, g, bb, .3))
	
	sns.swarmplot(x="variable", y="value", data=plot_df,palette =color_dict)
	unit=0.01
	y=myMax+unit*1
	h=unit*2
	plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c="black")
	plt.text(0.5, y+h+unit, "Paired T-test: %.2E" % scipy.stats.ttest_rel(a,b).pvalue, ha='center', va='bottom', color="black")
	plt.ylim(myMin-unit*5,myMax+0.1)
	plt.ylabel(ylabel)
	plt.savefig("%s.pdf"%(output), bbox_inches='tight')

def main():

	color_dict = {}
	color_dict['All_variants']="#fa6525"
	color_dict['GWAS_only']="#00b8a5"

	df = pd.read_csv("snp_data.tsv",sep="\t")
	target = "HbF"
	
	df = df.fillna(0)
	mu = df[target].mean()
	sigma = df[target].std()
	print (mu,sigma)
	df['label'] = [define_high_low(x,mu,sigma) for x in df[target]]
	
	df = df[df['label']>=0] 
	Y = df['label']
	print (Y[Y==0].shape)
	print (Y[Y==1].shape)
	X = df.drop([target,'label'],axis=1)
	sel_columns = []
	for c in X.columns:
		if not "chr" in c:
			sel_columns.append(c)
	X1 = X[sel_columns]

	model,params=sklearn_RF()
	
	auROC_list_a=[]
	auROC_list_b=[]
	auPRC_list_a=[]
	auPRC_list_b=[]
	df_list = []
	
	for i in range(100):
		ddf,a,c = simple_CV_evaluation(model,params,X,Y)
		ddf.to_csv("%s_prediction.csv"%("ALL"),index=False)
		plot_top_features(dp(model),X,Y,"ALL")
		ddf['label']="All_variants"
		auROC_list_a+=a
		auPRC_list_a+=c
		ddf2,b,d = simple_CV_evaluation(model,params,X1,Y)
		ddf2.to_csv("%s_prediction.csv"%("GWAS"),index=False)
		plot_top_features(dp(model),X,Y,"GWAS")
		ddf2['label']="GWAS_only"
		df_list.append(ddf)
		df_list.append(ddf2)
		auROC_list_b+=b
		auPRC_list_b+=d
	
	plot_df = pd.concat(df_list)
	plot_auROC_multi(plot_df,color_dict)
	plot_auPRC_multi(plot_df,color_dict)
	boxplot_paired_t_test(auROC_list_a,auROC_list_b,color_dict,"Area under ROC curve","auROC_boxplot")
	boxplot_paired_t_test(auPRC_list_a,auPRC_list_b,color_dict,"Area under Precision-Recall curve","auPRC_boxplot")

if __name__ == "__main__":
	main()







