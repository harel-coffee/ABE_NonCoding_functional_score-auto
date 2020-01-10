import numpy as np
import os
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
try:
	from StringIO import StringIO ## for Python 2
except ImportError:
	from io import StringIO ## for Python 3
import uuid
from joblib import Parallel, delayed

import argparse
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
from janggu.data import Cover
try:
	from StringIO import StringIO ## for Python 2
except ImportError:
	from io import StringIO ## for Python 3
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
import warnings
from sklearn.metrics import roc_curve,roc_auc_score,average_precision_score,accuracy_score
import warnings
warnings.filterwarnings('ignore')
# warnings.simplefilter(action='ignore', category=FutureWarning)
from sklearn.ensemble import RandomForestRegressor,GradientBoostingRegressor,RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics.scorer import make_scorer
from sklearn.model_selection import train_test_split
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
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge,Lars,BayesianRidge
from copy import deepcopy as dp
"""

Feature extraction (Top motif scores)

1. using janggu get DNA one-hot

3. read meme get motif PWMs in both strands

4. scan motifs get score_list, max(pos_strand,neg_strand)

with tree-based methods, we don't need to do normalization here

5. for each seq, get top N scores from (4) and their footprint score (given their positions), get adjusted score

Dependency
----------

meme (to get motif revcomp)
bedtools (to get fasta sequences for gkm_svm)

python library
--------------
janggu (tensorflow + keras)
biopython
sklearn
joblib


"""

	
def read_fasta(f):
	my_dict = {}
	for r in SeqIO.parse(f, "fasta"):
		my_dict[r.id] = str(r.seq).upper()
	return my_dict  
	
def read_motif(meme_file):
	revcomp_file = "/tmp/"+str(uuid.uuid4())
	os.system("meme-get-motif -rc -all %s > %s"%(meme_file,revcomp_file))
	original_motif_label = "++original++"
	revcomp_motif_label = "--revcomp--"
	
	dict1 = parse_meme(meme_file,label=original_motif_label)
	dict2 = parse_meme(revcomp_file,label=revcomp_motif_label)
	myDict = {}
	for k in dict1:
		motif_name = k.replace(original_motif_label,"")
		myDict[motif_name]=[dict1[k].T.values,dict2[k.replace(original_motif_label,revcomp_motif_label)].T.values]
	return myDict
	
def parse_meme(file,label=""):
	"""function to read meme file to pd.DataFrame"""
	lines = open(file).readlines()
	i = 0
	myDict = {}
	while i < len(lines):
		myList = lines[i].strip().split()
		if len(myList) < 1:
			i = i + 1
			continue
		if myList[0] == "MOTIF":
			if lines[i+1].strip() == "":
				desc = lines[i+2].strip().split()
				flag = True
			else:
				desc = lines[i+1].strip().split()
				flag = False
			try:
				motifLength = int(desc[5])
			except:
				print (desc)
				i = i+1
				continue
			if flag:
				myString = "\n".join(map(lambda x:"\t".join(x.strip().split()),lines[i+3:i+3+motifLength])).replace("  "," ")
				df = pd.read_csv(StringIO(myString), sep="\t",header=None)
				df.columns=['A','C','G','T']
				myDict[myList[1]+label] = df
				if df.shape[0] != motifLength or df.shape[1] !=4:
					print ("something is wrong")
				i = i+3+motifLength
				continue
			else:
				myString = "\n".join(map(lambda x:"\t".join(x.strip().split()),lines[i+2:i+2+motifLength])).replace("  "," ")
				df = pd.read_csv(StringIO(myString), sep="\t",header=None)
				df.columns=['A','C','G','T']
				myDict[myList[1]+label] = df
				i = i+2+motifLength 
				if df.shape[0] != motifLength or df.shape[1] !=4:
					print ("something is wrong")				
				continue
		i = i+1
	return myDict
	
def motif_scan(s,m):
	## s, m are numpy array
	## s.shape = L*4
	## m.shape = 4*W
	L = s.shape[0]
	W = m.shape[1]
	score_list = []
	for i in range(L-W):
		sub = np.matmul(s[i:i+W,:],m)
		# if i < 3:
			# print ("DNA seq",s[i:i+W,:])
			# print ("motif",m)
			# print ("mapping score: ",np.trace(sub))
		score_list.append(np.trace(sub))
	return score_list

def DNA_motif_scan(DNA_array,m1,m2):
	score_list = []
	# print (m1)
	# print (m2)
	for i in range(DNA_array.shape[0]):
		score_list_1 = motif_scan(DNA_array[i,:,:],m1)
		# print ("score_list_1",score_list_1)
		score_list_2 = motif_scan(DNA_array[i,:,:],m2)
		# print ("score_list_2",score_list_2)
		for j in range(len(score_list_1)):
			if score_list_2[j] > score_list_1[j]:
				score_list_1[j] = score_list_2[j]
		score_list.append(score_list_1)
	# print (score_list)
	out = np.array(score_list)
	print ("DNA scanning out shape",out.shape)
	return out
	

def get_roi(myList):
	## roi is region of interest, term used by janggu
	# chr19:13180899-13180900+
	# strand = [list(x)[-1] for x in myList]
	strand = [x[-1] for x in myList]
	# print (strand)
	chr = [x[:-1].split(":")[0] for x in myList]
	start = [int(x[:-1].split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x[:-1].split(":")[-1].split("-")[1]) for x in myList]
	roi_A = []
	roi = []
	for i in range(len(chr)):
		roi_A.append([chr[i],start[i],end[i],myList[i],".",strand[i]])
		roi.append([chr[i],start[i],end[i]])
	return roi_A,roi
	
def get_high_low_data(input,pos_cutoff,neg_cutoff):
	df = pd.read_csv(input,index_col=0)
	# pos = df[df['HbFBase']>=pos_cutoff].index.tolist()
	pos = df[df['HbFBase']>pos_cutoff].index.tolist()
	neg = df[df['HbFBase']<=neg_cutoff].index.tolist()
	print ("Pos size %s. Neg size %s"%(len(pos),len(neg)))
	return df.loc[pos+neg],pos,neg   

def roi2fasta(roi,genome_fa,flank):
	df = pd.DataFrame(roi)
	df[1] = df[1]-flank
	df[2] = df[2]+flank
	df.to_csv("tmp.bed",sep="\t",header=False,index=False)
	os.system("bedtools getfasta -fi %s -fo tmp.fa -bed tmp.bed -s -name"%(genome_fa))
	seq = read_fasta("tmp.fa")
	os.system("rm tmp.fa tmp.bed")
	return seq

## Define parameters

# high_hbf = 50
high_hbf = 0
low_hbf = 0
input = "Editable_A_scores.combined.scores.csv"
flank = 100
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
bw_file="/home/yli11/Projects/Li_gRNA/footprint/H1_H2_GM12878_Tn5_bw/Hudep2.bw"
meme_file = "selected_motifs.meme"
top_n=5 # number of features for each motif 

## read data

data,high,low = get_high_low_data(input,high_hbf,low_hbf)
roi_A,roi = get_roi(high+low)
seq = roi2fasta(roi_A,refgenome,flank)
test = pd.DataFrame.from_dict(seq,orient='index')
data['seq'] = test[0]

# 1. using janggu get DNA one-hot
## get one-hot data and ATAC feature matrix
dna_A = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_A,flank=flank)
Tn5 = Cover.create_from_bigwig('bigwig_coverage',bigwigfiles=bw_file,roi=roi,binsize=1,stepsize=1,flank=flank)  

## ReShape
dna_A=np.reshape(dna_A,(len(high+low),flank*2+1,4))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1))

## get motif PWM, 3. read meme get motif PWMs in both strands   
motifs = read_motif(meme_file)

# 4. scan motifs get score_list, max(pos_strand,neg_strand)
score_list_A = Parallel(n_jobs=-1)(delayed(DNA_motif_scan)(dna_A,motifs[m][0],motifs[m][1]) for m in motifs)


def get_footprint_score(s,l,footprint_score):
	flanking=2
	# print (s,l)
	left_start = s-flanking
	# print ("left_start:",left_start)
	if left_start >= 0:
		left = list(footprint_score[left_start:s])
	else:
		left = [np.nan]
	right_end = s+l+flanking
	# print ("right_end:",right_end)
	# print ("len(footprint_score):",len(footprint_score))
	if right_end <= len(footprint_score):
		right = list(footprint_score[s+l:right_end])
	else:
		right = [np.nan]		
	flanking = np.nanmean(left+right)
	# print ("left",left,"right",right)
	# print ("flanking",flanking,"left+right",left+right)
	occ = np.nanmean(footprint_score[s:s+l])
	# print ("all:",footprint_score[s:s+l],"occ:",occ)
	return flanking - occ

def get_top_n_motif_scores(score_list,top_n):
	"""score_list.shape = L * 1
	
	return
	------
	
	pos, value list
	"""
	return score_list.argsort()[-top_n:],score_list[score_list.argsort()[-top_n:]]
	
# 5. for each seq, get top N scores from (4) and their footprint score (given their positions), get adjusted score
def get_adjusted_motif_score(motif_score,footprint_score,n):
	"""motif_score and footprint_score are same shape, N * L"""
	out = []
	# print ("motif_score",motif_score)
	motif_length = footprint_score.shape[1] - motif_score.shape[1]
	for i in range(motif_score.shape[0]):
		pos,value =  get_top_n_motif_scores(motif_score[i],n)
		# print ("pos,:",pos)
		# print ("value,:",value)
		FOS_list = [get_footprint_score(s,motif_length,footprint_score[i]) for s in pos]
		# print ("FOS_list:",FOS_list)
		value = [value[i]*FOS_list[i] for i in range(len(value))]
		out.append(value)
	return out

adjusted_scores = Parallel(n_jobs=-1)(delayed(get_adjusted_motif_score)(motif_score,bw_values,top_n) for motif_score in score_list_A)
	

def set_col_names(motifs,top_n,label):
	out = []
	for i in motifs:
		for j in range(top_n):
			out.append("%s_%s_%s"%(label,i,j))
	return out

## get feature table
adjusted_scores = np.array(adjusted_scores)
adjusted_scores = np.swapaxes(adjusted_scores,0,1)
adjusted_scores = adjusted_scores.reshape((len(high+low),top_n*len(motifs)))
adjusted_scores = pd.DataFrame(adjusted_scores)
adjusted_scores.columns = set_col_names(motifs,top_n,"motif_footprint_score")
adjusted_scores.index = high+low
df = pd.concat([adjusted_scores,data],axis=1)
# df.to_csv("ML_data.csv")
df.to_csv("all_A_features.csv")






