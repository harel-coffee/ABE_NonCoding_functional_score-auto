import pandas as pd
pd.__version__
import sklearn
print sklearn.__version__
sklearn.__version__
import numpy
import numpy.core.multiarray
import numpy as np
-np.log2(0.0001)
import pandas as pd
df = pd.read_csv("BIOGRID-ALL-3.5.169.tab2.txt",sep="\t",index_col=0)
df.shape
df[df['Organism Interactor A'] == 9606].shape
df[(df['Organism Interactor A'] == 9606) & (df['Organism Interactor B'] == 9606)].shape
len(list(set(df['Entrez Gene Interactor A'].tolist()+df['Entrez Gene Interactor B'].tolist())))
import re
re.search("TT","ASFSDDDTT")
re.search("TT","ASFSDDDTT").start()
a="ASFSDDDTT"
a[7]
a[7:9]
re.search("TTF","ASFSDDDTT").start()
re.search("TTF","ASFSDDDTT")
re.search("DT","ASFSDDDTT").start()
a[4:6]
a[:6]
a=open("test","wt")
print (>>a,"s")
print>>a,"s"
print ("a",a)
a.close()
a="asd"
a+"i"
7406968 / 4
import pandas as pd
df = pd.read_csv("out.tsv",sep="\t",header=None)
df.head()
df[(df[2]==100)&(df[3]==20)]
df2 = df[(df[2]==100)&(df[3]==20)]
df2.shape
df2
df2.head()
df2.nunique()
df2[0].nunique()
df2[0].count_values()
df2[0].value_counts()
df2[df2[0]=='chr11:5705000-5705150%CTCTGTAAAATGGACCAATC']
df2[df2[0]=='chr11:5705000-5705150%CTCTGTAAAATGGACCAATC'][df2[1] == 'chr11']
df2.shape
df2.to_csv('gRNA_locations.csv')
import pandas as pd
df = pd.read_csv("out.tsv",index_col=0,header=None)
df = pd.read_csv("out.tsv",index_col=0,header=None,sep="\t")
df.head()
df = pd.read_csv("out.tsv",header=None,sep="\t")
df.head()
df[(df[2]==100)&(df[3]==23)]
df[(df[2]>95)&(df[3]==23)]
df = df[(df[2]>95)&(df[3]==23)]
df.head()
df.sort_values(2,ascending=False)
df.sort_values(2,ascending=False).head()
df.sort_values(2).head()
df.describe()
df.to_csv("gRNA_locations.csv")
import pandas as pd
df = pd.read_csv("out.tsv",header=None,sep="\t")
df = df[(df[2]>95)&(df[3]==23)]
df
df.head()
df[0].nunique()
import pandas as pd
df = pd.read_csv("out.tsv",sep="\t")
df.head()
import joblib
import markov_clustering as mc
import networkx as nx
import random
# number of nodes to use
numnodes = 200
# generate random positions as a dictionary where the key is the node id and the value
# is a tuple containing 2D coordinates
positions = {i:(random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
# use networkx to generate the graph
network = nx.random_geometric_graph(numnodes, 0.3, pos=positions)
# then get the adjacency matrix (in sparse form)
matrix = nx.to_scipy_sparse_matrix(network)
result = mc.run_mcl(matrix)           # run MCL with default parameters
clusters = mc.get_clusters(result)    # get clusters
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
mc.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="silver")
plt.savefig("test2.png")
import markov_clustering as mc
import networkx as nx
import random
import matplotlib.pyplot as plt
plt.switch_backend('agg')
HSC = nx.read_gml("G_HSC_filter_self_loops.gml")
matrix = nx.to_scipy_sparse_matrix(HSC)
result = mc.run_mcl(matrix)           # run MCL with default parameters
clusters = mc.get_clusters(result)    # get clusters
clusters
mc.draw_graph(matrix, clusters, node_size=10, with_labels=False, edge_color="silver")
plt.savefig("mcl_default_HSC.png")
for inflation in [i / 10 for i in range(15, 26)]:
    result = mc.run_mcl(matrix, inflation=inflation)
    clusters = mc.get_clusters(result)
    Q = mc.modularity(matrix=result, clusters=clusters)
    print("inflation:", inflation, "modularity:", Q)
for inflation in [i / 10 for i in range(20, 26)]:
    result = mc.run_mcl(matrix, inflation=inflation)
    clusters = mc.get_clusters(result)
    Q = mc.modularity(matrix=result, clusters=clusters)
    print("inflation:", inflation, "modularity:", Q)
import readline
readline.write_history_file("history_mcl.txt")
import markov_clustering as mc
import networkx as nx
import random
# number of nodes to use
numnodes = 200
# generate random positions as a dictionary where the key is the node id and the value
# is a tuple containing 2D coordinates
positions = {i:(random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
# use networkx to generate the graph
network = nx.random_geometric_graph(numnodes, 0.3, pos=positions)
# then get the adjacency matrix (in sparse form)
matrix = nx.to_scipy_sparse_matrix(netw
)
import markov_clustering as mc
import networkx as nx
import random
# number of nodes to use
numnodes = 200
# generate random positions as a dictionary where the key is the node id and the value
# is a tuple containing 2D coordinates
positions = {i:(random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
# use networkx to generate the graph
network = nx.random_geometric_graph(numnodes, 0.3, pos=positions)
# then get the adjacency matrix (in sparse form)
matrix = nx.to_scipy_sparse_matrix(network)
result = mc.run_mcl(matrix)           # run MCL with default parameters
clusters = mc.get_clusters(result)    # get clusters
clusters
HSC = nx.read_gml("G_HSC_filter_self_loops.gml")
matrix = nx.to_scipy_sparse_matrix(HSC)
result = mc.run_mcl(matrix)           # run MCL with default parameters
clusters = mc.get_clusters(result)    # get clusters
clusters
HSC
HSC.nodes()
matrix
G.nodes()
HSC.nodes()
len(HSC.nodes())
HSC.nodes()[0]
HSC.nodes()[1]
HSC.nodes().keys()
HSC.nodes().keys()[0]
HSC.nodes()
dir(HSC.nodes())
HSC.nodes().keys()
dir(HSC.nodes())
HSC.nodes().values()
HSC.nodes().keys()
HSC.nodes().values()
HSC.nodes()
HSC.nodes()[1]
dir(HSC.nodes())
list(G.nodes())
list(HSC.nodes())
result = mc.run_mcl(matrix,inflation=1)           # run MCL with default parameters
result = mc.run_mcl(matrix,inflation=1.001)           # run MCL with default parameters
cluster
clusters
HSC.nodes[6789]
HSC.nodes['6789']
HSC.nodes()
clusters
for i in clusters:
	print (i)
	print (list(i))
for i in clusters:
	print (list(i))
import markov_clustering as mc
from joblib import Parallel, delayed
import joblib'
import joblib
my_length = []
with open("SRR5997581.fastq") as f:
	for line in f:
		if line[0] == "@":
			line = line.strip()
			my_length.append(int(line.split("=")[-1]))
with open("SRR5997581.fastq") as f:
	for line in f:
		if line[0] == "+":
			line = line.strip()
			my_length.append(int(line.split("=")[-1]))
import numpy as np
np.mean(my_length)
f = open("test","wb")
print >>f,"asd"
f.close()
f = open("test","w")
print >>f,"asd"
import yagmail
yag = SMTP("liyc.stjude@gmail.com", oauth2_file="liyc.stjude.json")
yag = yagmail.SMTP("liyc.stjude@gmail.com", oauth2_file="liyc.stjude.json")
import yagmail
yag = yagmail.SMTP("liyc.stjude@gmail.com", oauth2_file="liyc.stjude.json")
import smtplib
dir(smtplib)
dir(smtplib.SMTP)
import os
os.system("ls")
yagmail.__version__
import mageck
print mageck.__path__
print mageck.__file__
mageck.__file__
import pygsheets
client = pygsheets.authorize(service_file='/home/yli11/Documents/Stjude-pipeline-report-345fc6038f45.json')
sh = client.open('spreadsheet-title')
wks = sh.sheet1
sh = client.open('spreadsheet-title')
sh = client.create('liyc-test')
wks = sh.sheet1
wks.update_value('A1', "Numbers on Stuff")
wks
sh = client.open('liyc-test')
wks = sh.sheet1
wks.update_value('A1', "Numbers on Stuff")
import pygsheets
pygsheets.__version__
import pygsheets
pygsheets.__version__
import pandas
import pandas as pd
df = pd.read_csv("HUDEP2_KLF1_combine_vs_input_peaks.revcomp.fa.fimo",sep="\t")
df['mapping_name'] = df['#pattern name'] + df['sequence name']
df=df.sort_values('p-value').drop_duplicates('mapping_name').sort_index()
print (df[df['#pattern name'] == '0_24.985'].shape)
print (df[df['#pattern name'] == '72_4.378'].shape)
print (df[df['#pattern name'] == '77_4.021'].shape)
print (df[df['#pattern name'] == '45_7.769'].shape)
import pandas as pd
df = pd.read_csv("HUDEP2_KLF1_combine_vs_input_peaks.revcomp.fa.fimo",sep="\t")
df['mapping_name'] = df['#pattern name'] + df['sequence name']
df=df.sort_values('p-value').drop_duplicates('mapping_name').sort_index()
print (df[df['#pattern name'] == '0_24.985'].shape)
print (df[df['#pattern name'] == '72_4.378'].shape)
print (df[df['#pattern name'] == '77_4.021'].shape)
print (df[df['#pattern name'] == '45_7.769'].shape)
import pandas as pd
df = pd.read_csv("HUDEP2_KLF1_combine_vs_input_peaks.revcomp.fa.fimo",sep="\t")
df['mapping_name'] = df['#pattern name'] + df['sequence name']
df=df.sort_values('p-value').drop_duplicates('mapping_name').sort_index()
print (df[df['#pattern name'] == '0_24.985'].shape)
print (df[df['#pattern name'] == '72_4.378'].shape)
print (df[df['#pattern name'] == '77_4.021'].shape)
print (df[df['#pattern name'] == '45_7.769'].shape)
import pandas as pd
import Dropbox
import dropbox
c = dropbox.client.DropboxClient(auth_token)
auth_token = 'vqpF8D3MTfAAAAAAAAAACxMLUVEGUMYqZJL5YIgb-sfiHgJC8mFWgxQlvihuB9TS'
c = dropbox.client.DropboxClient(auth_token)
dir(dropbox)
c = dropbox.dropbox(auth_token)
c = dropbox.Dropbox(auth_token)
c.share("/test.bed")
dropbox.__version__
c.sharing_create_shared_link_with_settings("/test.bed").url
import pyexamples
pyexamples.py_hello(b"world")
pyexamples.py_hello("world")
import idr
idr.__file__
import matplotlib
import pandas as pd
df = pd.read_csv("h1_vs_h2.transcript.final.combined.tpm.csv")
df.head()
df['xx'] = (df['treatment_mean']+1)/(df['control_mean']+1)
import numpy as np
df['xx'] = np.log2(df['xx'])
df.head()
module load python/2.7.13
from selenium import webdriver
webpage = r"http://www.crisprindelphi.design/single" # edit me
driver = webdriver.Chrome()
driver.get(webpage)
driver = webdriver.Chrome("/home/yli11/chromedriver")
import pandas as pd
df = pd.read_csv("test.csv",index_col=0)
import seaborn as sns
sns.clustermap(df)
import seaborn as sns
sns.__version__
import numpy as np
np.__version__
import matplotlib
matplotlib.__version__
import pandas as pd
import glob
files = glob.glob("*.tsv")
def to_csv(x):
	df = pd.read_csv(x,sep="\t",index_col=0)
	df.to_csv("%s.csv"%(x))
for f in files:
	to_csv(f)
import readlines
import readline
readline.write_history_file("to_csv.py")
import pandas as pd
label = pd.read_csv("labels.tsv",sep="\t")
label.head()
label.index = label['gRNA sequence']
df = pd.read_csv("gRNA.loci306.center.bed_homer.tsv",sep="\t")
df.head()
df.index = df['PeakID']
df.head()
sub = df.loc[label.index.tolist()]
sub.shape
label.shape
df2 = pd.concat([label,sub],axis=1)
df2.shape
df2.head()
from feature_selector import FeatureSelector
df2['ratio'] = df['HbF%']/df['control']
df2['ratio'] = df2['HbF%']/df2['control']
df2['class'] = df2['ratio'] >= 1.5
df2['class']
df2['class'] = df2['class'].astype(int)
df2['class']
train = df2.copy()
train_labels = train['class']
fs = FeatureSelector(data = train, labels = train_labels)
fs.identify_missing(missing_threshold=0.6)
fs.identify_single_unique()
fs.ops['single_unique']
fs.unique_stats.sample(5)
fs.identify_collinear(correlation_threshold=0.9)
correlated_features = fs.ops['collinear']
fs.plot_collinear()
plt.savefig("correlation.png")
import matplotlib.pylab as plt
plt.savefig("correlation.png")
fs.identify_zero_importance(task = 'classification', eval_metric = 'auc', n_iterations = 10, early_stopping = True)
fs.data_all.head(10)
train.columns
fs.plot_feature_importances(threshold = 0.99, plot_n = 12)
plt.savefig("importance.png")
fs.feature_importances.head(10)
fs.feature_importances.head(100)
train.columns.tolist()
rm_cols = ['HbF%', 'control', 'gRNA sequence', 'PeakID', 'Chr', 'Start', 'End', 'Strand', 'Peak Score', 'Focus Ratio/Region Size', 'Annotation', 'Detailed Annotation', 'Distance to TSS', 'Nearest PromoterID', 'Entrez ID', 'Nearest Unigene', 'Nearest Refseq', 'Nearest Ensembl', 'Gene Name', 'Gene Alias', 'Gene Description', 'Gene Type','ratio', 'class']
train = train.drop(rm_cols,axis=1)
fs = FeatureSelector(data = train, labels = train_labels)
fs.identify_all(selection_params = {'missing_threshold': 0.6, 'correlation_threshold': 0.7, 'task': 'classification', 'eval_metric': 'auc', 'cumulative_importance': 0.5})
train_removed_all_once = fs.remove(methods = 'all', keep_one_hot = True)
fs.feature_importances
fs.feature_importances.to_csv("importance.csv")
train['ChIP_seq_TF_HUDEP2_GATA1_S2_vs_input_peaks.narrowPeak']
df2.to_csv("df2.csv")
import readline
readline.write_history_file("history.py")
import pandas as pd
import matplotlib.pylab as plt
df = pd.read_csv("gRNA.loci306_upstream_400_600_with_label.csv")
df.head()
df.boxplot('ChIP_seq_TF_HUDEP2_sh276_KLF1_vs_input_peaks.narrowPeak',by='class')
plt.savefig("boxplot.png")
test = pd.read_csv("gRNA.loci306_upstream_400_600_importance.csv")
test.head()
test = pd.read_csv("gRNA.loci306_upstream_400_600_importance.csv",index_col=0)
test
import umap
from sklearn.datasets import load_digits
digits = load_digits()
digits.data
df = pd.DataFrame(digits.data)
import pandas as pd
df = pd.DataFrame(digits.data)
df.head()
embedding = umap.UMAP().fit_transform(df)
embedding = umap.UMAP().fit_transform(digits.data)
df.values
embedding = umap.UMAP().fit_transform(df.values)
import plotly.graph_objs as go
import plotly.io as plio
import plotly.graph_objs as go
go.Scatter.__doc__
go.Scatter().__doc__
print (go.Scatter.__doc__)
dir(go.Scatter)
dir(go)
from xvfbwrapper import Xvfb
vdisplay = Xvfb(width=1280, height=740)
vdisplay.start()
import tensorflow as tf
sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
import tensorflow
import tensorflow as tf
sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
sess = tf.compat.v1.Session(config=tf.ConfigProto(log_device_placement=True))
tf.test.is_gpu_available()
import tensorflow as tf
tf.test.is_gpu_available
tf.test.is_gpu_available()
import tensorflow as tf
tf.test.is_gpu_available()
import tensorflow as tf
tf.test.is_gpu_available()
import os
import numpy as np
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import Cover
from janggu.data import ReduceDim
from janggu.data import SqueezeDim
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
os.environ['JANGGU_OUTPUT'] = ''
order = 5
REFGENOME ='/home/yli11/Data/Human/hg19/fasta/hg19.fa'
pos = 'pos.bed'
neg = 'neg.bed'
Bioseq.create_from_refgenome(name='dna',refgenome=REFGENOME,roi=pos)
Bioseq("dna")
dna = Bioseq.create_from_refgenome(name='dna',refgenome=REFGENOME,roi=pos)
dna
dna.shape
dna[0]
pow(3,4)
pow(2,3)
import readline
readline.write_history_file("test.py")
import os
import numpy as np
from keras import Model
from keras import backend as K
from keras.layers import Conv2D
from keras.layers import Dense
from keras.layers import GlobalAveragePooling2D
from keras.layers import Input
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import Cover
from janggu.data import ReduceDim
from janggu.layers import DnaConv2D
from sklearn.metrics import roc_auc_score
REFGENOME = resource_filename('janggu', 'resources/pseudo_genome.fa')
# ROI contains regions spanning positive and negative examples
ROI_TRAIN_FILE = resource_filename('janggu', 'resources/roi_train.bed')
ROI_TEST_FILE = resource_filename('janggu', 'resources/roi_test.bed')
# PEAK_FILE only contains positive examples
PEAK_FILE = resource_filename('janggu', 'resources/scores.bed')
DNA = Bioseq.create_from_refgenome('dna', refgenome=REFGENOME,
                                   roi=ROI_TRAIN_FILE,
                                   binsize=200,
                                   order=order,
                                   cache=True)
order = 3
DNA = Bioseq.create_from_refgenome('dna', refgenome=REFGENOME,
                                   roi=ROI_TRAIN_FILE,
                                   binsize=200,
                                   order=order,
                                   cache=True)
DNA = Bioseq.create_from_refgenome('dna', refgenome=REFGENOME,
                                   roi=ROI_TRAIN_FILE,
                                   binsize=199,
                                   order=order,
                                   cache=True)
DNA = Bioseq.create_from_refgenome('dna', refgenome=REFGENOME,
                                   roi=ROI_TRAIN_FILE,
                                   binsize=199,
                                   order=order,
                                   cache=False)
DNA
DNA.shape
order = 5
REFGENOME ='/home/yli11/Data/Human/hg19/fasta/hg19.fa'
pos = 'pos.bed'
neg = 'neg.bed'
Bioseq.create_from_refgenome(name='dna',refgenome=REFGENOME,roi=pos)
Bioseq("dna")
dna = Bioseq.create_from_refgenome(name='dna',refgenome=REFGENOME,roi=pos)
dna.shape
DNA = Bioseq.create_from_refgenome('dna', refgenome=REFGENOME,
                                   roi=pos,
                                   binsize=199,
                                   stepsize=4,
                                   order=order,
                                   cache=False)
DNA.shape
LABELS = ReduceDim(Cover.create_from_bed('peaks', roi=ROI_TRAIN_FILE,
                               bedfiles=PEAK_FILE,
                               binsize=200,
                               resolution=200,
                               cache=True,
                               storage='sparse'))
os.system("module load bedtools")
LABELS = ReduceDim(Cover.create_from_bed('peaks', roi=ROI_TRAIN_FILE,
                               bedfiles=PEAK_FILE,
                               binsize=200,
                               resolution=200,
                               cache=True,
                               storage='sparse'))
import os
import numpy as np
from keras import Model
from keras import backend as K
from keras.layers import Conv2D
from keras.layers import Dense
from keras.layers import GlobalAveragePooling2D
from keras.layers import Input
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import Cover
from janggu.data import ReduceDim
from janggu.layers import DnaConv2D
from sklearn.metrics import roc_auc_score
REFGENOME = resource_filename('janggu', 'resources/pseudo_genome.fa')
# ROI contains regions spanning positive and negative examples
ROI_TRAIN_FILE = resource_filename('janggu', 'resources/roi_train.bed')
ROI_TEST_FILE = resource_filename('janggu', 'resources/roi_test.bed')
# PEAK_FILE only contains positive examples
PEAK_FILE = resource_filename('janggu', 'resources/scores.bed')
LABELS = ReduceDim(Cover.create_from_bed('peaks', roi=ROI_TRAIN_FILE,
                               bedfiles=PEAK_FILE,
                               binsize=200,
                               resolution=200,
                               cache=True,
                               storage='sparse'))
LABELS
LABELS.shape
LABELS[0]
import pickle
f = open("AT847CE","rb")
pickle.load(f)
from inspect import signature
"{}".format("asd")
"%s".format("asd")
import scipy.stats as ss
"Percentage of {0:.1%} in {1} overlapped with {2}".format(0.2,1,1)
import pandas as pd
df = pd.read_csv("~/dirs/hudep_footprint/H2.bed",sep="\t",header=None)
df.head()
df.describe()
df[df[4]>=100].to_csv("H2_FT_filtered.bed",sep="\t",header=False,index=False)
range(2)
"asd%s"%(2)
import Bio
import pandas as pd
import plotly
import plotly.express as px
px.colors.qualitative.swatches()
dir(px.colors)
px.colors.validate_colors
px.colors.validate_colors()
px.colors.validate_colors(12)
px.colors.n_colors
?px.colors.n_colors
?px.colors.n_colors()
px.colors.n_colors()
px.colors.n_colors(2)
dir(px.colors.n_colors)
dir(px.colors.n_colors())
dir(px.colors.n_colors)
px.colors.colorbrewer
dir(px.colors.colorbrewer)
px.colors.colorbrewer.Blues
import networkx
import pandas as pd
df = pd.read_csv("rna_seq_yli11_2019-07-02.transcript.tpm.csv")
df.head()
df.columns
df.columns = ['Transcript ID', 'Gene ID', 'Gene Name','Hudep1_S5','Hudep1_S3','Hudep1_S1','Hudep2_S4','Hudep2_S6','Hudep2_S2']
df = df.drop_duplicates("Gene ID")
import pandas as pd
df = pd.read_csv("/rgs01/project_space/tsaigrp/Genomics/common/projects/ATAC-seq/merged/joined_sites_CHANGEseq_GUIDEseq_casoffinder_20190618.tsv"
,sep="\t")
df.head()
df.columns
for c in df.columns:
	print (df[c].nunique())
for c in df.columns:
	print (c,df[c].nunique())
df.sample.head()
df['sample']head()
df['sample'].head()
df['target'].head()
df[df['target']==df['offtarget_sequence_GUIDEseq']]
df[df['distance_GUIDEseq']==0]
df[df['distance_GUIDEseq']==0][['target','offtarget_sequence_GUIDEseq']]
guide = df[df.distance_GUIDEseq>0]
ids = guide['genome_coordinate']
guide[ids.isin(ids[ids.duplicated()])]
guide[ids.isin(ids[ids.duplicated()])]['GUIDEseq_reads']
guide[ids.isin(ids[ids.duplicated()])][['site','offtarget_sequence_GUIDEseq','GUIDEseq_reads']]
guide[ids.isin(ids[ids.duplicated()])][['sample','offtarget_sequence_GUIDEseq','GUIDEseq_reads']]
df.columns
selected_cols_change = ['chr', 'start', 'end','sample','CHANGEseq_reads','strand']
selected_cols_guide = ['chr', 'start', 'end','sample','GUIDEseq_reads','strand']
CHANGEseq_index_list = df[df.distance_CHANGEseq>0].index.tolist()
GUIDEseq_index_list = df[df.distance_GUIDEseq>0].index.tolist()
change_not_in_guide = list(set(CHANGEseq_index_list)-set(GUIDEseq_index_list))
change_not_in_guide_df = df.loc[change_not_in_guide]
ids = change_not_in_guide_df['genome_coordinate']
change_not_in_guide_df[ids.isin(ids[ids.duplicated()])]
change_not_in_guide_df.groupby('genome_coordinate')['CHANGEseq_reads'].mean()
mean_readCounts = change_not_in_guide_df.groupby('genome_coordinate')['CHANGEseq_reads'].mean()
change_not_in_guide_df.index = change_not_in_guide_df['genome_coordinate']
change_not_in_guide_df = change_not_in_guide_df.drop_duplicates("genome_coordinate")
change_not_in_guide_df['mean_readCounts'] = mean_readCounts
selected_cols = ['chr', 'start', 'end','sample','mean_readCounts','strand']
change_not_in_guide_df[selected_cols].to_csv("CHANGEseq_not_GUIDEseq_mean_dis_above0_readCounts.bed",sep="\t",index=False,header=False)
guide_not_in_guide_df = df.loc[GUIDEseq_index_list]
mean_readCounts = guide_not_in_guide_df.groupby('genome_coordinate')['guideseq_reads'].mean()
guide_not_in_guide_df.index = guide_not_in_guide_df['genome_coordinate']
guide_not_in_guide_df = guide_not_in_guide_df.drop_duplicates("genome_coordinate")
guide_not_in_guide_df['mean_readCounts'] = mean_readCounts
guide_not_in_guide_df[selected_cols].to_csv("guideseq_mean_dis_above0_readCounts.bed",sep="\t",index=False,header=False)
guide_not_in_guide_df
guide_not_in_guide_df[selected_cols]
guide_not_in_guide_df = df.loc[GUIDEseq_index_list]
mean_readCounts = guide_not_in_guide_df.groupby('genome_coordinate')['guideseq_reads'].mean()
guide_not_in_guide_df.index = guide_not_in_guide_df['genome_coordinate']
guide_not_in_guide_df = df.loc[GUIDEseq_index_list]
guide_not_in_guide_df.shape
guide_not_in_guide_df.groupby('genome_coordinate')['guideseq_reads'].mean()
guide_not_in_guide_df = df.loc[GUIDEseq_index_list]
mean_readCounts = guide_not_in_guide_df.groupby('genome_coordinate')['GUIDEseq_reads'].mean()
guide_not_in_guide_df.index = guide_not_in_guide_df['genome_coordinate']
guide_not_in_guide_df = guide_not_in_guide_df.drop_duplicates("genome_coordinate")
guide_not_in_guide_df['mean_readCounts'] = mean_readCounts
guide_not_in_guide_df[selected_cols].to_csv("guideseq_mean_dis_above0_mean_readCounts.bed",sep="\t",index=False,header=False)
import pandas as pd
df = pd.read_csv("/rgs01/project_space/tsaigrp/Genomics/common/projects/ATAC-seq/merged/joined_sites_CHANGEseq_GUIDEseq_casoffinder_20190618.tsv"
,sep="\t")
df.columns
df[df['sample']=='CXCR4_site_1']
df[(df['sample']=='CXCR4_site_1')&(df['GUIDEseq_reads']>=0)]
df[(df['sample']=='CXCR4_site_1')&(df['GUIDEseq_reads']>=0)]['GUIDEseq_reads']
df[(df['sample']=='CXCR4_site_1')&(df['GUIDEseq_reads']>=0)][['start','GUIDEseq_reads']]
import umap
"1234567"[0:2]
"1234567"[1:3]
import pandas as pd
df = pd.read_csv("gRNA.loci306.bed",sep="\t",header=None)
start=0
length=13
def row_apply(x):
	chr = x[0]
	gRNA_start = x[1]
	gRNA_end = x[2]
	gRNA = x[3]
	value = x[4]
	strand = x[5]
	edit_start = gRNA_start+start
	edit_end = gRNA_start+start+length
	edit_gRNA_seq = gRNA[start:(start+length)]
	return [chr,edit_start,edit_end,edit_gRNA_seq,value,strand]
	
x=df.head().apply(row_apply,axis=1)
x
x.to_list()
pd.DataFrame(df.head().apply(row_apply,axis=1).tolist())
import pysam
import tabix
tb = tabix.open("Eigen_hg19_noncoding_annot_chr10.tab.bgz")
tb.query("chr10", 1000000, 1250000)
tb.query("chr10", 10000000, 12500000)
tb.query("chr10", 0,1)
tb.querys("chr10")
min(0,1)
import pyBigWig
file="~/dirs/hudep_ATAC_merged/atac_seq_yli11_2019-08-02/bam_files/atac_seq_footprint_yli11_2019-08-05/H2.bw"
bw = pyBigWig.open(file)
import os
os.path.isfile(file)
import pyBigWig
file="/home/yli11/dirs/hudep_ATAC_merged/atac_seq_yli11_2019-08-02/bam_files/atac_seq_footprint_yli11_2019-08-05/H2.bw"
os.path.isfile(file)
import os
os.path.isfile(file)
bw = pyBigWig.open(file)
import pandas as pd
def row_apply(x,bw,extend_length):
	return bw.values(x[0], min(0,x[1]-extend_length), x[2]+extend_length)
df = pd.read_csv("test.bed",sep="")
df = pd.read_csv("test.bed",sep="\t",header=None)
df
bw.values('chr1',1167431,1167441)
import matplotlib
import joblib
from matplotlib_venn import venn3
53.0/206
53.0/306
import pandas as pd
df = pd.read_csv("HUDEP2_G_WT_ZBTB7A_merged.vs.HUDEP2__G__WT_ZBTB7A_ChIPminusseq_merged_Input_peaks.rmblck.narrowPeak.fa.fimo",sep="\t")
df.head()
df[df['#pattern name'].str.contains("CCGGGAG")]
df[df['#pattern name'].isin(["CCGGGAG"])]
tmp=df[df['#pattern name'].isin(["CCGGGAG"])]
tmp[tmp['matched sequence'].isin(['CCGGAAG'])]
3644.0/79926
import pandas as pd
import umap
import pandas as pd
import umap
import numpy as np
df = pd.DataFrame(np.random.randint(0,100,size=(100, 4)), columns=list('ABCD'))
df
sub_df = df.loc[0:50]
sub_df
umap_obj = umap.UMAP(n_components=n,random_state=0)
umap_obj = umap.UMAP(n_components=2,random_state=0)
umap_obj.fit(sub_df)
s1 = umap.transform(sub_df)
s1 = umap_obj.transform(sub_df)
s1
pd.DataFrame(umap_obj.transform(sub_df)).head()
pd.DataFrame(umap_obj.transform(df)).head()
pd.DataFrame(umap_obj.transform(sub_df)).head()
sub_umap = pd.DataFrame(umap_obj.transform(sub_df))
def stand(dff):
    df = dff.copy()
    for c in df.columns:
        df[c] = (df[c] - df[c].mean())/df[c].std(ddof=0)
    print (df.head())
stand(sub_umap)
all_umap = pd.DataFrame(umap_obj.transform(df))
stand(all_umap)
umap.UMAP.fit_transform(sub_df)
umap.UMAP().fit_transform(sub_df)
pd.DataFrame(umap.UMAP().fit_transform(sub_df)).head()
umap_obj.fit_transform(sub_df)
pd.DataFrame(umap_obj.fit_transform(sub_df)).head()
pd.DataFrame(umap_obj.transform(sub_df)).head()
pd.DataFrame(umap_obj.transform(df)).head()
pd.DataFrame(umap_obj.transform(df)).describe()
pd.DataFrame(umap_obj.transform(sub_df)).describe()
import sklearn
man = sklearn.manifold.Isomap()
man.fit(sub_df)
pd.DataFrame(man.transform(sub_df)).head()
pd.DataFrame(man.transform(df)).head()
man = sklearn.manifold.locally_linear_embedding()
umap.__version__
import pandas as pd
df = pd.read_csv()
files = glob.glob("*.prob.csv")
import glob
files = glob.glob("*.prob.csv")
def parse_df(f):
	df = pd.read_csv(f.replace("../","").replace("/","_")+".prob.csv",index_col=0)
	return df
df_list = [parse_df(f) for f in files]
import pandas as pd
import glob
files = glob.glob("*.prob.csv")
files
def parse_df2(f):
	df = pd.read_csv(f,index_col=0)
	return df
df_list = [parse_df2(f) for f in files]
df_list[0]
df1 = pd.DataFrame(dict(x=np.random.randn(100), y=np.random.randint(0, 5, 100)))
df2 = pd.DataFrame(dict(x=np.random.randn(100), y=np.random.randint(0, 5, 100)))
import numpy as np
df1 = pd.DataFrame(dict(x=np.random.randn(100), y=np.random.randint(0, 5, 100)))
df2 = pd.DataFrame(dict(x=np.random.randn(100), y=np.random.randint(0, 5, 100)))
df1
df2
df_concat = pd.concat((df1, df2))
df_concat
from functools import reduce
d = reduce(lambda x, y: x.add(y, fill_value=0), df_list)
d
d = reduce(lambda x, y: x.add(y), df_list)
d
import pandas as pd
df = pd.read_csv("plot_df.csv")
df.head()
df = pd.read_csv("plot_df.csv",index_col=0)
df
df.corr()
import pandas as pd
df
import pandas as pd
df = pd.DataFrame()
df[1]=[1,2,3]
df[2]=[3,4,5]
df
df.at[1,1]
import pandas as pd
df = pd.DataFrame()
df[1]=[1,2,3]
df[2]=[3,4,5]
df.sum()
df[df[1]<2].sum()
df[df[1]<2].sum(axis=1)
df[df[1]<=2].sum(axis=1)
df[df[1]<=2].sum(axis=0)
df[df[1]<=2][2].sum()
df[df[1]<=-1][2].sum()
import numpy as np
pd.DataFrame(np.full((20,20), 0.0)).shape
pd.DataFrame(np.full((20,40), 0.0)).shape
pd.DataFrame(np.full((40,20), 0.0)).shape
from itertools import chain, combinations
def all_subsets(ss):
    return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))
for subset in all_subsets([0,1,2]):
	print (subset)
for subset in all_subsets([0,1,2]):
	print (list(subset))
pd.DataFrame(columns=list(range20))
import pandas as pd
df = pd.DataFrame(columns=list(range20))
df = pd.DataFrame(columns=list(range(20)))
df
df = pd.DataFrame(columns=list(range(20))+['freq'])
df
df.loc[0]=[0]*21
df
0.9**(1.0/3)
0.9**(1.0/20)
261889/
261889/13094
import pandas as pd
df = pd.read_csv("gRNA_peak.distance.bed",sep="\t")
df.head()
df = pd.read_csv("gRNA_peak.distance.bed",sep="\t",header=None)
df.head()
import pandas as pd
# bedtools intersect -a loci.bed -b H2_ideas.bed -wao > loci_intersect.bed
df = pd.read_csv("loci_intersect.bed",sep="\t",header=None)
df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
df = df.sort_values(9,ascending=False)
df = df.drop_duplicates('name')
all = pd.read_csv("sum.out",sep="\t")
all['name'] = all['chr']+":"+all['start'].astype(str)+"-"+all['end'].astype(str)
df.index = df.name
all.index = all.name
df[6].unique
df[6].unique()
df[6].unique().tolist()
for c in df[6].unique().tolist():
	all[c]=0
	for i in all.index:
		if df.at[i,6]==c:
			all.at[i,c]=1
all.head()
all.to_csv("combined.data.csv")
from sklearn.ensemble import RandomForestRegressor,GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics.scorer import make_scorer
from sklearn.base import TransformerMixin
from sklearn.datasets import make_regression
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import LinearRegression, Ridge
import re
from glasbey import Glasbey
import plotly.graph_objs as go
import plotly.io as plio
import plotly.express as px
import plotly
plotly.io.orca.config.executable = "/home/yli11/.conda/envs/py2/bin/orca"
import plotly.io as pio
pio.orca.config.use_xvfb = True
import datetime
import uuid
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
from joblib import Parallel, delayed
from os.path import isfile,isdir
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.manifold import TSNE
import matplotlib
import numpy as np
import scipy
import glob
import sys
import matplotlib.pyplot as plt
import os
import numpy as np
import getpass
import argparse
from matplotlib_venn import venn3,venn2
sys.setrecursionlimit(99999)
import pandas as pd
df = pd.read_csv("combined.data.csv")
df.sum()
import pandas as pd
df = pd.read_csv("combined.data.csv")
df.head()
df = pd.read_csv("combined.data.csv",index_col=0)
df[df[0]==1]
df[df["0"]==1]
df['atac'].describe()
import pandas as pd
df = pd.read_csv("combined.data.csv",index_col=0)
df.corr()
df.corr().to_csv("corr.csv")
import scipy
dir(scipy)
from scipy.signal import find_peaks
import pandas as pd
a=[['chr11', 5264348, 5264369], ['chr11', 5264371, 5264383], ['chr11', 5264385, 5264424]]
pd.DataFrame(a)
a=pd.DataFrame(a)
import numpy as np
import pandas as pd
def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df
import sys
import glob
def read(x):
	print (x)
	df = pd.read_csv(x,sep="\t",header=None)
	df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
	# df = df[df[4]!=0]
	df['%s'%(x)] = df[4]
	print (df[2]-df[1]).describe()
	df = df.drop([0,1,2,3,4],axis=1)
	df = df.set_index('name')
	# print df.shape
	print (df.head())
	return df
def merge(substring):
	df_list = [read(x) for x in glob.glob("*%s"%(substring))]
	df = pd.concat(df_list,axis=1)
	print (df.head())
	return df
df = merge(".bdg.sorted.bed12.smooth.bed")
import numpy as np
import pandas as pd
a=[1,2,3,"a"]
a
a.remove("a")
a
import pandas as pd
df = pd.read_csv("U2OS_CHANGE_seq_EMX1_CRL538_truncated_identified_matched.txt",sep="\t",header=None)
df.head()
df.corr()
df = pd.read_csv("U2OS_CHANGE_seq_EMX1_CRL538_truncated_identified_matched.txt",sep="\t",header=None)
df.head()
df[0]
df[1]
df.head()
df.corr()
df.corr().to_csv()
df.corr().to_csv("correlation.csv")
lines = open("intersect.bed").readlines()
len(lines)
lines[-1]
import os
os.cwd()
os.getcwd()
import Bio
import pandas as pd
df = pd.read_csv("yli11_2019-10-30_f6a8e/AAGAGAGAAATTACATCTGT.fimo",sep="\t",comment="#",header=None)
df.columns = ["motif_name","seq_name","start","end","strand","score","p-value","q-value","matched_sequence"]
df['start'] -= 1
df['end'] -= 1
df = df.drop(['score','q-value'],axis=1)
window_start = 15+1
window_end = 15+3+6+2
	df = df[df['end']>=window_start]
	df = df[df['start']<=window_end]
df = df[df['end']>=window_start]
df = df[df['start']<=window_end]
df
df = df.sort_values('p-value')
df['pair'] = df['motif_name']+df["seq_name"]
df = df.drop_duplicates('pair',keep='first')
df = df.drop(['pair'],axis=1)
df.head()
df['seq_name'].unique().tolist()
df.columns
df['info'] = df['motif_name'].astype(str)+","+df['seq_name'].astype(str)+","+df['start'].astype(str)+","+df['end'].astype(str)+","+df['strand']+","+df['matched_sequence']+","+df['p-value'].astype(str) ## this is used for generating figures
df.index = df['seq_name']
seq_list = df['seq_name'].unique().tolist()
# print (seq_list)
## p-value table and match info table
p_value_table = pd.DataFrame()
mapping_info_table = pd.DataFrame()
for m in motif_list:
	p_value_table[m],mapping_info_table[m] = get_motif_column(df,m,seq_list,fimo_cutoff)
	
def get_motif_column(df,m,index_list,default_p_value):
	p_value_list = df[df['motif_name']==m]['p-value'].get(index_list,default_p_value)
	info_list = df[df['motif_name']==m]['info'].get(index_list,"NA")
	return p_value_list,info_list
for m in motif_list:
	p_value_table[m],mapping_info_table[m] = get_motif_column(df,m,seq_list,fimo_cutoff)
motif_list = df['motif_name'].unique().tolist()
for m in motif_list:
	p_value_table[m],mapping_info_table[m] = get_motif_column(df,m,seq_list,fimo_cutoff)
for m in motif_list:
	p_value_table[m],mapping_info_table[m] = get_motif_column(df,m,seq_list,0.05)
p_value_table.index = seq_list
p_value_table
p_value_table.shape
p_value_table.columns
df
df.columns
mapping_info_table
p_value_table.append(pandas.Series(name='ref'))
p_value_table.append(pd.Series(name='ref'))
p_value_table
import myvariant
mv = myvariant.MyVariantInfo()
mv.query('dbsnp.rsid:rs58991260', fields='dbsnp')
import pandas as pd
df = pd.read_csv("~/Data/Human/hg19/SNP151/snp151_clean.bed",sep="\t",header=None)
q = pd.read_csv("gwas_trait_rsID.list",header=None)
q.head()
df.head()
df[df[3].isin(q)]
df.index = df[3]
df.loc[q[0].tolist()].shape
a=df.loc[q[0].tolist()]
a.head()
a.dropna().shape
a.dropna().to_csv("rs_SNP.tsv")
df.to_feather("SNP151.feather")
a.dropna().to_csv("rs_SNP.tsv",sep="\t",header=False,index=False)
import pandas as pd
snp_trait = pd.read_csv("gwas_trait_rsID.bed",sep="\t",header=None)
snp_bed = pd.read_csv("gwas_trait_rsID_hg19.bed",sep="\t",header=None)
snp_bed
snp_trait.head()
snp_trait.index = snp_trait[0]
snp_bed.index = snp_bed[3]
snp_trait['chr'] = snp_bed.loc[snp_trait.index.tolist()][0]
snp_trait['start'] = snp_bed.loc[snp_trait.index.tolist()][1]
snp_trait['end'] = snp_bed.loc[snp_trait.index.tolist()][2]
snp_trait
snp_trait.dropna()
snp_bed
import readline
readline.write_history_file("get_SNP_trait_bed.py")
import pandas as pd
snp_trait = pd.read_csv("gwas_trait_rsID.bed",sep="\t",header=None)
snp_bed = pd.read_csv("gwas_trait_rsID_hg19.bed",sep="\t",header=None)
snp_bed
snp_trait.head()
snp_trait.index = snp_trait[0]
snp_bed.index = snp_bed[3]
snp_trait['chr'] = snp_bed.loc[snp_trait.index.tolist()][0]
snp_trait['start'] = snp_bed.loc[snp_trait.index.tolist()][1]
snp_trait['end'] = snp_bed.loc[snp_trait.index.tolist()][2]
snp_trait
snp_trait.dropna()
snp_bed
import pandas as pd
snp_trait = pd.read_csv("gwas_trait_rsID.bed",sep="\t",header=None)
snp_bed = pd.read_csv("gwas_trait_rsID_hg19.bed",sep="\t",header=None)
snp_bed
snp_trait.head()
snp_trait.index = snp_trait[0]
snp_bed.index = snp_bed[3]
snp_trait['chr'] = snp_bed.loc[snp_trait.index.tolist()][0]
snp_trait['start'] = snp_bed.loc[snp_trait.index.tolist()][1]
snp_trait['end'] = snp_bed.loc[snp_trait.index.tolist()][2]
snp_trait
snp_trait.dropna()
snp_bed
snp_trait.dropna()[['chr','start','end',1]].to_csv("gwas_rsID_trait.bed",sep="\t",header=False,index=False)
import pandas as pd
snp_trait = pd.read_csv("gwas_trait_rsID.bed",sep="\t",header=None)
snp_bed = pd.read_csv("gwas_trait_rsID_hg19.bed",sep="\t",header=None)
snp_bed
snp_trait.head()
snp_trait.index = snp_trait[0]
snp_bed.index = snp_bed[3]
snp_trait['chr'] = snp_bed.loc[snp_trait.index.tolist()][0]
snp_trait['start'] = snp_bed.loc[snp_trait.index.tolist()][1].astype(int)
snp_trait['end'] = snp_bed.loc[snp_trait.index.tolist()][2].astype(int)
snp_trait
snp_trait.dropna()
snp_bed
snp_trait.dropna()[['chr','start','end',1]].to_csv("gwas_rsID_trait.bed",sep="\t",header=False,index=False)
import pandas as pd
snp_trait = pd.read_csv("gwas_trait_rsID.bed",sep="\t",header=None)
snp_bed = pd.read_csv("gwas_trait_rsID_hg19.bed",sep="\t",header=None)
snp_bed
snp_trait.head()
snp_trait.index = snp_trait[0]
snp_bed.index = snp_bed[3]
snp_trait['chr'] = snp_bed.loc[snp_trait.index.tolist()][0]
snp_trait['start'] = snp_bed.loc[snp_trait.index.tolist()][1]
snp_trait['end'] = snp_bed.loc[snp_trait.index.tolist()][2]
snp_trait
snp_trait.dropna()
snp_bed
df = snp_trait.dropna()
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)
df[['chr','start','end',1]].to_csv("gwas_rsID_trait.bed",sep="\t",header=False,index=False)
import pandas as pd
df = pd.read_csv("combined.gwas.hg19.bed",sep="\t",header=None)
df.head()
df.shape
df = df.dropna()
df.shape
a="123_"
a[-1:]
a[:-1]
def remove_last_underline(x):
	if x[-1] == "_":
		return x[:-1]
	return x
	
df[3] = [remove_last_underline(x) for x in df[3]]
df['pos'] = df[0]+":"+df[2].astype(str)
for s,d in df.groupby(3):
	d['pos'].to_csv("%s.list"%(s),index=False,header=False)
import pandas as pd
df = pd.read_csv("combined.gwas.hg19.bed",sep="\t",header=None)
df.head()
def remove_last_underline(x):
	if x[-1] == "_":
		return x[:-1]
	return x
	
df[3] = [remove_last_underline(x) for x in df[3]]
df['pos'] = df[0]+":"+df[2].astype(str)
df.head()
a="asd"
b="123asd3"
a in b
def check_substring(x,myList):
	tmp = dp(myList)
	tmp.remove(x)
	for j in tmp:
		if x in j:
			return 1
	return 0
myList = df[3].unique().tolist()
len(myList)
values = [check_substring(x,myList) for x in myList]
from copy import deepcopy as dp
values = [check_substring(x,myList) for x in myList]
map_substring = pd.DataFrame()
map_substring[0]=values
map_substring.index = myList
map_substring.head()
df[4] = map_substring.loc[df[3].tolist()][0]
map_substring.loc[df[3].tolist()].head()
df[4] = map_substring.loc[df[3].tolist()][0].tolist()
df.head()
map_substring[map_substring[0]==0].shape
map_substring.shape
for s,d in df.groupby(3):
	d['pos'].to_csv("%s.list"%(s),index=False,header=False)
import glob
glob.glob("*combine*")
import pandas as pd
df = pd.read_csv("combined.gwas.hg19.bed",sep="\t",header=None)
df.head()
def remove_last_underline(x):
	if x[-1] == "_":
		return x[:-1]
	return x
	
df[3] = [remove_last_underline(x) for x in df[3]]
df['pos'] = df[0]+":"+df[2].astype(str)
df.head()
for s,d in df.groupby(3):
	if df.shape[0] < 10:
		continue
	d['pos'].to_csv("%s.list"%(s),index=False,header=False)
"N"*20
a="TGCAGCACAGCAGCGAGGAAGGG"
a[:-3]
len(a[:-3])
import pandas as pd
df = pd.read_csv("test.bed",sep="\t",header=None)
df
df[0]!=df[1]
df[0] = [x.upper() for x in df[0]]
df[1] = [x.upper() for x in df[1]]
df[0]!=df[1]
sum(df[0]!=df[1])
import pandas as pd
df = pd.read_csv("test.bed",sep="\t",header=None)
df2 = pd.read_csv("candidate_gRNA.bed",sep="\t",header=None)
df.head()
df2.head()
df2.name = df2[0]+":"+df2[1].astype(str)+"-"+df2[2].astype(str)
df2['name'] = df2[0]+":"+df2[1].astype(str)+"-"+df2[2].astype(str)
df2.groupby(3)['name'].agg(', '.join)
df3 = pd.DataFrame(df2.groupby(3)['name'].agg(', '.join))
df3.shape
df3.head()
df['name'] = df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)
df.index = df.name
df['other matches'] = df3['name']
df.head()
df3df3['name']
df3['name']
df.index = df[3]
df
df['other matches'] = df3['name']
df
def remove_self_match(r):
	myList = r['other matches'].split(", ")
	myList.remove(r['name'])
	return ",".join(myList)
df['other matches'] = df.apply(remove_self_match,axis=1)	
df.head()
def get_numOfftarget(r):
	myList = r['other matches'].split(",")
	return len(myList)
df['numOffTargets'] = df.apply(get_numOfftarget,axis=1)
df
df.head()
def get_numOfftarget(r):
	if r['other matches'] == ""
		return 0
	myList = r['other matches'].split(",")
	return len(myList)
def get_numOfftarget(r):
	if r['other matches'] == ""
		return 0
	myList = r['other matches'].split(",")
	return len(myList)
def get_numOfftarget(r):
	if r['other matches'] == "":
		return 0
	myList = r['other matches'].split(",")
	return len(myList)
df['numOffTargets'] = df.apply(get_numOfftarget,axis=1)
df
df.head()
import pandas as pd
df = pd.read_csv("idr_peaks.fimo",sep="\t")
names = df['sequence name'].unique().tolist()
print (len(names))
df = pd.DataFrame()
df[0] = [x.split(":")[0] for x in names]
df[1] = [x.split(":")[1].split("-")[0] for x in names]
df[2] = [x.split(":")[1].split("-")[1] for x in names]
df.to_csv("idr_peaks_with_NFIX_motifs.bed",sep="\t",header=False,index=False)
2690.0/2707
import pandas as pd
df = pd.read_csv("NGG_gRNAs.off_targets.info.csv")
df.head()
df.numOffTargets.unique()
df = df[df['numOffTargets']==0]
df.shape
df.head()
df.to_csv("unique_matches.csv",index=False)
import pandas as pd
df = pd.read_csv("unique_matches.csv")
df.head()
df['G'] = [x[0] for x in df[3]]
df['G'] = [x[0] for x in df["3"]]
df.head()
df = df[df['G']=="G"]
df.shape
df = df[(df[4]>=0.3)&(df[4]<=0.7)]
df = df[(df["4"]>=0.3)&(df["4"]<=0.7)]
df.shape
df[[0,1,2,3,4,5]].to_csv("candidate_gRNA_strict.bed",sep="\t",header=False,index=False)
df[[str(x) for x in range(6)]].to_csv("candidate_gRNA_strict.bed",sep="\t",header=False,index=False)
import readline
readline.write_history_file("filter.py")
import pandas as pd
df = pd.read_csv("head fimo.cuts.freq.txt.raw_fimo.bed",sep="\t",header=None)
import pandas as pd
df = pd.read_csv("fimo.cuts.freq.txt.raw_fimo.bed",sep="\t",header=None)
df.head()
df.mean()
df.median()
df.mean()
import pandas as pd
df = pd.read_csv("fimo.txt",sep="\t")
df.head()
df[df['q-value']<=1e-6]
df[df['q-value']<=1e-5]
df[df['q-value']<=1e-4]
df[df['q-value']<=1e-4].shape
df[df['q-value']<=1e-3].shape
df[df['q-value']<=1e-3].head()
df[['sequence name','start','end']].to_csv("fimo_q_value_cutoff1e_3.bed",sep="\t",header=False,index=False)
df[['sequence name','start','stop']].to_csv("fimo_q_value_cutoff1e_3.bed",sep="\t",header=False,index=False)
import pandas as pd
df = pd.read_csv("/rgs01/project_space/tsaigrp/Genomics/common/projects/ATAC-seq/merged/joined_sites_CHANGEseq_GUIDEseq_casoffinder_20190919.tsv",sep="\t")
df.describe()
import pandas as pd
df  = pd.DataFrame()
df.head()
a=[0,1]
a.index(1)
a.index(0)
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
from mlxtend.classifier import StackingCVClassifier
import umap
import warnings
from sklearn.metrics import roc_curve,roc_auc_score
from sklearn.datasets import load_iris
from mlxtend.classifier import StackingCVClassifier
from mlxtend.feature_selection import ColumnSelector
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
from xgboost import XGBClassifier
from sklearn.feature_selection import RFECV
def sklearn_RF(par=False):
	est = RandomForestClassifier(n_estimators=800,random_state=0)
	if par:
		est = RandomForestClassifier(**par)
	myDict = {}
	myDict['max_depth']=[30,None]
	# myDict['criterion']=['gini','entropy'] 
	# myDict['max_features']=['sqrt',None] 
	myDict['min_samples_leaf']=[1,5,9,13] 
	# myDict['min_samples_split']=[2,8,14,20] 
	# myDict['min_weight_fraction_leaf']=[0,0.001,0.01,0.1] 
	# myDict['min_impurity_decrease']=[0,0.05,0.1,0.2]    
	myDict['warm_start']=[True,False]    
	return est, myDict
df = pd.read_csv("one_hot_chr_complete_feature_table.csv")
df.head()
df = pd.read_csv("one_hot_chr_complete_feature_table.csv",index_col=0)
df.head()
X = df.drop(["class"],axis=1)
Y = df['class']
model,params=sklearn_RF()
selector = RFECV(dp(model), step=1, cv=StratifiedKFold(3),scoring="roc_auc")
selector.fit(X,Y)
selector.fit(X,Y,n_jobs=-1)
selector = RFECV(dp(model), step=1, cv=StratifiedKFold(3),scoring="roc_auc",n_jobs=1)
selector = RFECV(dp(model), step=1, cv=StratifiedKFold(3),scoring="roc_auc",n_jobs=-1)
selector.fit(X,Y,n_jobs=-1)
selector.fit(X,Y)
selector.get_support()
X.columns[selector.get_support()]
X.columns[selector.get_support()].tolist()
selector.support_ 
selector.n_features
selector.n_features_
import autosklearn.classification
cls = autosklearn.classification.AutoSklearnClassifier()
import pandas as pd
df = pd.read_csv("one_hot_chr_complete_feature_table.csv",index_col=0)
df.head()
target = 'class'
X = df.drop([target],axis=1)
Y = df[target]
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=0)
X_train.head()
cls = autosklearn.classification.AutoSklearnClassifier(n_jobs=-1)
cls.fit(X_train, y_train)
cls.fit(X_train.values, y_train)
cls.fit(X_train.values, y_train.values)
import sklearn.model_selection
import sklearn.datasets
import sklearn.metrics
X, y = sklearn.datasets.load_digits(return_X_y=True)
X_train, X_test, y_train, y_test = \
        sklearn.model_selection.train_test_split(X, y, random_state=1)
automl = autosklearn.classification.AutoSklearnClassifier()
automl.fit(X_train, y_train)
X = df.drop([target],axis=1)
Y = df[target]
X_train
y_train
X.values
Y.values
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=0)
X_train
X_train.values
y_train.values
automl.fit(X_train, y_train)
cls = autosklearn.classification.AutoSklearnClassifier()
cls.fit(X_train, y_train)
cls = autosklearn.classification.AutoSklearnClassifier(n_jobs=10)
cls.fit(X_train, y_train)
predictions = cls.predict(X_test)
predictions
predictions = cls.predict_proba(X_test)
predictions
from sklearn.metrics import roc_curve,roc_auc_score
roc_auc_score(y_test.tolist(),[x[1] for x in predictions])
y_test
roc_auc_score(y_test.tolist(),[x[1] for x in predictions])
import pandas as pd
bed = pd.read_csv("ensembl_hg19_genes.bed.gz",sep="\t",header=None)
gene = pd.read_csv("ensembl2geneName.bed.gz",sep="\t",header=None)
gene.head()
gene = pd.read_csv("ensembl2geneName.bed.gz",sep="\t")
gene.head()
bed.head()
gene = pd.read_csv("ensembl2geneName.bed.gz",sep="\t",index_col=0)
gene.head()
gene['value'].to_dict()
myDict = gene['value'].to_dict()
bed[3]=bed[3].replace(myDict)
bed[3]=bed[3].map(myDict)
bed.head()
bed.to_csv("hg19_gene.bed",sep="\t",header=False,index=False)
bed = bed.drop_duplicates(3)
bed.shape
bed.to_csv("hg19_gene_rmdup.bed",sep="\t",header=False,index=False)
32606595-502
from pkg_resources import resource_filename
from janggu.data import Bioseq
?resource_filename
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
roi="gRNA_bin_100.bed"
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi)
dna.shape
dna[0]
dna[0][0, :10, 0, :]
bw_file="/home/yli11/dirs/hudep_ATAC_merged/atac_seq_yli11_2019-08-02/bam_files/atac_seq_footprint_yli11_2019-08-05/H2.bw"
cover = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1)
from janggu.data import Cover
cover = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1)
cover[0]
cover.shape
dna.shape
13147*201
from janggu.data import ReduceDim
ReduceDim
new_cover = ReduceDim(cover)
new_cover.shape
cover.gindexer
cover.gindexer()
dir(cover.gindexer)
cover.gindexer.collapse
cover.gindexer.export_to_bed
cover.gindexer.export_to_bed()
cover.gindexer.export_to_bed("test.bed")
type(cover)
cover.shape
cover.head()
dir(cover)
cover.ndim
cover.name
cover.reshape
new_cover = ReduceDim(cover,(200,200))
cover[0]
dna[0][0, :10, 0, :]
dna[0][0, :10, 0, 1]
dna[0][0, :10, 0, :]
dna.shape
dna[0].shape
dna[0][0, :10, 0, 4]
dna[0][0]
dna[0][0].shape
dna[0][0][0].shape
dna[0][0][0]
dna[0][0, :10, 0, 4]
dna[0][0, :10, 0, :]
x=dna[0][0][0]
x
x.append(2)
dna[0][0, :10, 0, :]
np.append(dna[0][0][0],[[3]])
import numpy as np
np.append(dna[0][0][0],[[3]])
dna[0][0][0]
np.append(dna[0][0][0],[3])
np.append(dna[0][0][0],3)
cover.shape
dna.shape
x=np.reshape(cover,(13147,201))
x.shape
x
x[0]
x[0].shape
cover[0:10]
x[0][0:10]
dna.shape
x.shape
y = np.append(dna,x)
y.shape
y = np.concat(dna,x,axis=3)
y = np.concatenate(dna,x,axis=3)
x.shape
x=np.reshape(cover,(13147,201,1,1))
x.shape
x[0]
x[0][0]
x[0][0:10]
x[0][0:10,1,1]
x[0][0:10,0,0]
cover[0:10]
cover[0:10,0,0,0]
cover[0:10,0,0]
cover[0:10,0]
cover[0:10]
x[0][0:10,0,0]
x.shape
dna.shape
y = np.concatenate(dna,x,axis=3)
y = np.concatenate(dna,x,axis=1)
y = np.concatenate(dna,x,axis=2)
type(dna)
dna.shape
d = np.array(dna)
d.shape
y = np.stack((dna,x), axis=4)
dna.shape
x.shape
y = np.stack((dna,x), axis=3)
dir(dna)
dd = np.reshape(d,(13147,201,1,4))
dd.sha[e
;
dd.shape
y = np.concatenate(dd,x,axis=3)
y = np.stack((dd,x), axis=3)
dd.shape
x.shape
y = np.stack((dd,x), axis=-1)
dna.shape
x.shape
y = np.append(dna,x,axis=3)
y = np.append(dd,x,axis=3)
y.shape
y[0]
from keras.layers import Conv2D
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
(x_train, y_train), (x_test, y_test) = mnist.load_data()
y_train
y_train.shape
y_train = np.random.randint(1,size=(100))
y_train
y_train = np.random.randint(2,size=(100))
y_train
y_train.shape
x_train = y[0:100,:,:,:]
x_test = y[100:200,:,:,:]
y_train = np.random.randint(2,size=(100))
y_test = np.random.randint(2,size=(100))
array([0, 0, 0, 0, 0, 0, 0
num_classes=2
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
y_test
input_shape = (201, 1, 5)
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))
y.shape
yy = np.reshape(y,(13147,201,5))
yy.shape
yy[0][0]
y[0][0][0]
x_train = yy[0:100,:,:]
x_test = yy[100:200,:,:]
y_train = np.random.randint(2,size=(100))
y_test = np.random.randint(2,size=(100))
num_classes=2
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
# (13147, 201, 1, 5)
input_shape = (201, 5)
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))
yy = np.reshape(y,(13147,201,5,1))
yy[0][0]
x_train = yy[0:100,:,:,:]
x_test = yy[100:200,:,:,:]
y_train = np.random.randint(2,size=(100))
y_test = np.random.randint(2,size=(100))
num_classes=2
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
# (13147, 201, 1, 5)
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(MaxPooling2D(pool_size=(2, 2), dim_ordering="th"))
model.add(MaxPooling2D(pool_size=(2, 2), dim_ordering="tf"))
from __future__ import print_function
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
batch_size = 128
num_classes = 10
epochs = 12
# input image dimensions
img_rows, img_cols = 28, 28
# the data, split between train and test sets
(x_train, y_train), (x_test, y_test) = mnist.load_data()
if K.image_data_format() == 'channels_first':
    x_train = x_train.reshape(x_train.shape[0], 1, img_rows, img_cols)
    x_test = x_test.reshape(x_test.shape[0], 1, img_rows, img_cols)
    input_shape = (1, img_rows, img_cols)
else:
    x_train = x_train.reshape(x_train.shape[0], img_rows, img_cols, 1)
    x_test = x_test.reshape(x_test.shape[0], img_rows, img_cols, 1)
    input_shape = (img_rows, img_cols, 1)
x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
x_train /= 255
x_test /= 255
print('x_train shape:', x_train.shape)
print(x_train.shape[0], 'train samples')
print(x_test.shape[0], 'test samples')
# convert class vectors to binary class matrices
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(x_test, y_test))
score = model.evaluate(x_test, y_test, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2), dim_ordering="tf"))
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
# model.add(Dense(num_classes, activation='softmax'))
model.add(Lambda(lambda x: K.tf.nn.softmax(x)))
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='sigmoid'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
x_train = yy[0:100,:,:,:]
x_test = yy[100:200,:,:,:]
y_train = np.random.randint(2,size=(100))
y_test = np.random.randint(2,size=(100))
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
num_classes=2
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='sigmoid'))
# model.add(Lambda(lambda x: K.tf.nn.softmax(x)))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
x_train = yy[0:100,:,:,:]
x_test = yy[100:200,:,:,:]
y_train = np.random.randint(2,size=(100))
y_test = np.random.randint(2,size=(100))
num_classes=2
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
x_train = yy[0:100,:,:,:]
x_test = yy[100:200,:,:,:]
y_train = np.random.randint(2,size=(100))
y_test = np.random.randint(2,size=(100))
num_classes=2
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(2, activation='sigmoid'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
import readline 
readline.write_history_file("janggu_test.py")
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
roi="gRNA_bin_100.bed"
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi)
# dna.shape
# dna[0]
# dna[0][0, :10, 0, :]
bw_file="/home/yli11/dirs/hudep_ATAC_merged/atac_seq_yli11_2019-08-02/bam_files/atac_seq_footprint_yli11_2019-08-05/H2.bw"
from janggu.data import Cover
cover = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1)
x=np.reshape(cover,(13147,201,1,1))
d=np.reshape(dna,(13147,201,1,4))
d.shape
x.shape
cover[0]
cover[0:10]
x[0][0:10]
dna[0][0:10,0,:]
dna[0][0,0,:]
dna[0][1,0,:]
dna.shape
dna[0][0,0,:]
dna[0][0,0,0]
dna[0][0, :10, 0, :]
x.shape
x[0][0, :10, 0, :]
d.shape
d[0][0:10]
d[0][0,0:10]
d[0][0:10]
x[0][0:10]
y = np.append(d,x,axis=3)
       [[ 0.49785003]]
y.shape
from keras.layers import Conv2D
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
yy = np.reshape(y,(13147,201,5,1))
yy[0][0:10]
x_train = yy[0:100,:,:,:]
x_test = yy[100:200,:,:,:]
y_train = np.random.randint(2,size=(100))
y_test = np.random.randint(2,size=(100))
num_classes=2
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
input_shape = (201, 5,1)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(2, activation='sigmoid'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
score = model.evaluate(x_test, y_test, verbose=0,callbacks=['cor', 'mae', 'mse', 'var_explained'])
score = model.evaluate(x_test, y_test, verbose=0)
score
y_test.shape
from keras.callbacks import Callback
y_pred = model.predict(x_test)
y_pred
y_pred = model.predict_prob(x_test)
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(2, activation='softmax'))
import tensorflow as tf
tf.__version__
keras.__version__
from keras.layers.core import Dense,Dropout,Activation,Flatten,Lambda
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
# model.add(Dense(2, activation='sigmoid'))
model.add(Lambda(lambda x: K.tf.nn.softmax(x)))
_version__
'1.3.0'
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
model = Sequential()
model.add(Conv2D(32, kernel_size=(30, 3),
                 activation='relu',
                 input_shape=input_shape))
model.add(Conv2D(64, (30, 2), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
# model.add(Dense(2, activation='sigmoid'))
# model.add(Lambda(lambda x: K.tf.nn.softmax(x)))
model.add(Dense(2, activation = (tf.nn.softmax)))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train,
          batch_size=10,
          epochs=30,
          verbose=1,
          validation_data=(x_test, y_test))
y_pred = model.predict(x_test)
y_pred
import argparse
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
from keras.layers import Conv2D
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
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
import umap
import warnings
from sklearn.metrics import roc_curve,roc_auc_score
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
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
roi="gRNA_bin_100.bed"
roi = ['chr6:135505212-135505213', 'chr19:13140918-13140919']
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi)
roi
roi = ['chr6135505212135505213', 'chr191314091813140919']
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi)
roi = [['chr6',135505212,135505213]]
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi)
dna
dna[0]
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi)
roi = [['chr6',135505212,135505213,"x",'.',"-"]]
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi)
dna[0]
dna = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi,flank=100)
dna[0]
from keras.utils import plot_model
dir(plot_model)
plot_model()
plot_model.__name__
plot_model.__defaults__
plot_model.__doc__
import keras
from keras.layers import Input, LSTM, Dense
from keras.models import Model
tweet_a = Input(shape=(280, 256))
tweet_b = Input(shape=(280, 256))
shared_lstm = LSTM(64)
# When we reuse the same layer instance
# multiple times, the weights of the layer
# are also being reused
# (it is effectively *the same* layer)
encoded_a = shared_lstm(tweet_a)
encoded_b = shared_lstm(tweet_b)
# We can then concatenate the two vectors:
merged_vector = keras.layers.concatenate([encoded_a, encoded_b], axis=-1)
# And add a logistic regression on top
predictions = Dense(1, activation='sigmoid')(merged_vector)
# We define a trainable model linking the
# tweet inputs to the predictions
model = Model(inputs=[tweet_a, tweet_b], outputs=predictions)
model.compile(optimizer='rmsprop',
              loss='binary_crossentropy',
              metrics=['accuracy'])
shared_lstm = LSTM(64)
# When we reuse the same layer instance
# multiple times, the weights of the layer
# are also being reused
# (it is effectively *the same* layer)
encoded_a = shared_lstm(tweet_a)
encoded_b = shared_lstm(tweet_b)
shared_lstm = LSTM(64)
encoded_a = shared_lstm(tweet_a)
keras.__version__
import keras
from keras.layers import Input, LSTM, Dense
from keras.models import Model
tweet_a = Input(shape=(280, 256))
tweet_b = Input(shape=(280, 256))
shared_lstm = LSTM(64)
encoded_a = shared_lstm(tweet_a)
encoded_b = shared_lstm(tweet_b)
merged_vector = keras.layers.concatenate([encoded_a, encoded_b], axis=-1)
# And add a logistic regression on top
predictions = Dense(1, activation='sigmoid')(merged_vector)
# We define a trainable model linking the
# tweet inputs to the predictions
model = Model(inputs=[tweet_a, tweet_b], outputs=predictions)
model.compile(optimizer='rmsprop',
              loss='binary_crossentropy',
              metrics=['accuracy'])
import argparse
import matplotlib
matplotlib.use('agg')
import tensorflow as tf
import seaborn as sns
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
from keras.layers import Conv2D,LSTM
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
from janggu.data import Cover
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
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
from sklearn.metrics import roc_curve,roc_auc_score
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
from keras.utils import plot_model
from copy import deepcopy as dp
def get_neg_data(df,pos_list):
	neg_df = df[df['HbFBase']<=10]
	neg_index = []
	pos_position_list = pd.DataFrame(df.loc[pos_list].groupby('main-pos').size())
	for i in pos_position_list.index:
		neg_index += neg_df[neg_df['main-pos']==i].sample(n=pos_position_list.at[i,0]).index.tolist()
	return neg_index
def get_pos_neg_data(df):
	# get pos neg data
	pos_list = df[df['HbFBase']>=100].index.tolist()
	neg_list = get_neg_data(df,pos_list) 
	return pos_list,neg_list
## read data
freq = pd.read_csv("positional_model.csv",header=None)
freq = freq.fillna(1e-4)
freq = freq[1].tolist()
freq[13:20]=[1e-100]*7
freq = pd.DataFrame(freq)
df = pd.read_csv("feature_table.csv",index_col=0)
def get_main_pos(x):
	myList = []
	for i in x.split(", "):
		myList.append(int(i))  
	return freq.loc[myList][0].idxmax()
df['main-pos'] = [get_main_pos(y) for y in df.pos]
def get_neg_data(df,pos_list):
	neg_df = df[df['HbFBase']<=10]
	neg_index = []
	pos_position_list = pd.DataFrame(df.loc[pos_list].groupby('main-pos').size())
	for i in pos_position_list.index:
		neg_index += neg_df[neg_df['main-pos']==i].sample(n=pos_position_list.at[i,0]).index.tolist()
	return neg_index
def get_pos_neg_data(df):
	# get pos neg data
	pos_list = df[df['HbFBase']>=100].index.tolist()
	neg_list = get_neg_data(df,pos_list) 
	return pos_list,neg_list
def get_pos_neg_data(df):
	pos_list = df[df['HbFBase']>=100].index.tolist()
	neg_list = get_neg_data(df,pos_list) 
	return pos_list,neg_list
freq = pd.read_csv("positional_model.csv",header=None)
freq = freq.fillna(1e-4)
freq = freq[1].tolist()
freq[13:20]=[1e-100]*7
freq = pd.DataFrame(freq)
df = pd.read_csv("feature_table.csv",index_col=0)
def get_main_pos(x):
	myList = []
	for i in x.split(", "):
		myList.append(int(i))  
	return freq.loc[myList][0].idxmax()
df['main-pos'] = [get_main_pos(y) for y in df.pos]
pos,neg = get_pos_neg_data(df)
def get_roi(myList):
	# chr19:13180899-13180900
	chr = [x.split(":")[0] for x in myList]
	start = [int(x.split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x.split(":")[-1].split("-")[1]) for x in myList]
	roi_F = []
	roi_R = []
	for i in range(len(chr)):
		roi_F.append([chr[i],start[i],end[i],".",".","+"])
		roi_R.append([chr[i],start[i],end[i],".",".","-"])
	return roi_F,roi_R
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
bw_F_file="/home/yli11/Projects/Li_gRNA/footprint/H1_H2_GM12878_Tn5_bw/Hudep2_forward.wig.bw"
bw_R_file="/home/yli11/Projects/Li_gRNA/footprint/H1_H2_GM12878_Tn5_bw/Hudep2_reverse.wig.bw"
flank = 100
roi_F,roi_R = get_roi(pos+neg)
dna_F = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_F,flank=flank)
dna_R = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_R,flank=flank)
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=roi_F,
                                 binsize=1,
                                 stepsize=1,flank=flank)
def get_roi(myList):
	# chr19:13180899-13180900
	chr = [x.split(":")[0] for x in myList]
	start = [int(x.split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x.split(":")[-1].split("-")[1]) for x in myList]
	roi_F = []
	roi_R = []
	roi = []
	for i in range(len(chr)):
		roi_F.append([chr[i],start[i],end[i],".",".","+"])
		roi_R.append([chr[i],start[i],end[i],".",".","-"])
		roi.append([chr[i],start[i],end[i],".",".","-"])
	return roi_F,roi_R,roi
roi_F,roi_R,roi = get_roi(pos+neg)
dna_F[0].shape
dna_F[0][0,0:10,:,:]
dna_R[0][0,0:10,:,:]
dna_R[0][0,-10:,:,:]
dna_F[0][0,0:10,:,:]
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1,flank=flank)
roi
def get_roi(myList):
	# chr19:13180899-13180900
	chr = [x.split(":")[0] for x in myList]
	start = [int(x.split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x.split(":")[-1].split("-")[1]) for x in myList]
	roi_F = []
	roi_R = []
	roi = []
	for i in range(len(chr)):
		roi_F.append([chr[i],start[i],end[i],".",".","+"])
		roi_R.append([chr[i],start[i],end[i],".",".","-"])
		roi.append([chr[i],start[i],end[i]])
	return roi_F,roi_R,roi
roi_F,roi_R,roi = get_roi(pos+neg)
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1,flank=flank)
roi
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1,flank=1)
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1,flank=flank)
bw_file="/home/yli11/dirs/hudep_ATAC_merged/atac_seq_yli11_2019-08-02/bam_files/atac_seq_footprint_yli11_2019-08-05/H2.bw"
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1,flank=flank)
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=[['chr11', 4629102, 4629103]],
                                 binsize=1,
                                 stepsize=1,flank=flank)
roi
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1,flank=flank)
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=[['chr11', 4629102, 4629103], ['chr11', 4629103, 4629104], ['chr11', 4629107, 4629108], ['chr11', 4629112, 4629113], ['chr11', 5173162, 5173163], ['chr11', 5173167, 5173168], ['chr11', 5173168, 5173169], ['chr11', 5173171, 5173172], ['chr11', 5173173, 5173174], ['chr11', 5173174, 5173175], ['chr11', 5246880, 5246881], ['chr11', 5246884, 5246885], ['chr11', 5248322, 5248323], ['chr11', 5248326, 5248327], ['chr11', 5248327, 5248328], ['chr11', 5248328, 5248329], ['chr11', 5248329, 5248330], ['chr11', 5248331, 5248332], ['chr11', 5248341, 5248342], ['chr11', 5248369, 5248370], ['chr11', 5248373, 5248374], ['chr11', 5248374, 5248375], ['chr11', 5248450, 5248451], ['chr11', 5255910, 5255911], ['chr11', 5255914, 5255915], ['chr11', 5255916, 5255917], ['chr11', 5255918, 5255919], ['chr11', 5264357, 5264358], ['chr11', 5264363, 5264364], ['chr11', 5264364, 5264365], ['chr11', 5264366, 5264367], ['chr11', 5264367, 5264368], ['chr11', 5264372, 5264373], ['chr11', 5268537, 5268538], ['chr11', 5268543, 5268544], ['chr11', 5268546, 5268547], ['chr11', 5268547, 5268548], ['chr11', 5302014, 5302015], ['chr11', 5302018, 5302019], ['chr11', 5302020, 5302021], ['chr11', 5302021, 5302022], ['chr11', 5305925, 5305926], ['chr11', 5305929, 5305930], ['chr11', 5305936, 5305937], ['chr11', 5646429, 5646430], ['chr11', 5646430, 5646431], ['chr11', 5646431, 5646432], ['chr11', 5646435, 5646436], ['chr11', 5646442, 5646443], ['chr11', 5646445, 5646446], ['chr19', 12952049, 12952050]],
                                 binsize=1,
                                 stepsize=1,flank=flank)
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_F_file,
                                 roi=[['chr6', 135608445, 135608446], ['chr2', 57987101, 57987102], ['chr19', 13076623, 13076624], ['chr11', 5181541, 5181542], ['chr6', 135395042, 135395043]],
                                 binsize=1,
                                 stepsize=1,flank=flank)
cover_F = Cover.create_from_bigwig('bigwig_coverage',bigwigfiles=bw_F_file,roi=[['chr6', 135608445, 135608446]],binsize=1,stepsize=1,flank=flank)					
cover_F = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_file,
                                 roi=roi,
                                 binsize=1,
                                 stepsize=1,flank=flank)				
cover_R = Cover.create_from_bigwig('bigwig_coverage',
                                 bigwigfiles=bw_file,
                                 roi=roi_F,
                                 binsize=1,
                                 stepsize=1,flank=flank)
len(roi)
roi
len(roi_F)
len(pos)
bw_values_F=np.reshape(cover_F,(len(pos+neg),flank*2+1,1,1))
bw_values_R=np.reshape(cover_R,(len(pos+neg),flank*2+1,1,1))
dna_F=np.reshape(dna_F,(len(pos+neg),flank*2+1,4,1))
dna_R=np.reshape(dna_R,(len(pos+neg),flank*2+1,4,1))
import seya
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5 Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5 Reverse")
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
	shared_drop_out = Dropout(0.7)
	DNA_model_F = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_f)
	DNA_model_R = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_r)
	combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
	DNA_LSTM = LSTM(32)(combined_DNA)
	
	combined_LSTM = DNA_LSTM
	
	combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
	
	output = Dense(2, activation = "softmax")(combined_dense)
	return output
	
myModel = ATAC_model(201)
import argparse
import matplotlib
matplotlib.use('agg')
import tensorflow as tf
import seaborn as sns
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
from keras.layers import Conv2D,LSTM
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
from janggu.data import Cover
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D,Input,GlobalAveragePooling2D
from keras import backend as K
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
from sklearn.metrics import roc_curve,roc_auc_score
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
from keras.utils import plot_model
from copy import deepcopy as dp
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5 Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5 Reverse")
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
	shared_drop_out = Dropout(0.7)
	DNA_model_F = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_f)
	DNA_model_R = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_r)
	combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
	DNA_LSTM = LSTM(32)(combined_DNA)
	
	combined_LSTM = DNA_LSTM
	
	combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
	
	output = Dense(2, activation = "softmax")(combined_dense)
	return output
	
myModel = ATAC_model(201)
dna_length=100
dna_f = Input(shape=(dna_length,4,1),name="DNAForward")
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5 Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5 Reverse")
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
	shared_drop_out = Dropout(0.7)
	DNA_model_F = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_f)
	DNA_model_R = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_r)
	combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
	DNA_LSTM = LSTM(32)(combined_DNA)
	
	combined_LSTM = DNA_LSTM
	
	combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
	
	output = Dense(2, activation = "softmax")(combined_dense)
	return output
	
myModel = ATAC_model(201)
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
	shared_drop_out = Dropout(0.7)
	DNA_model_F = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_f)
	DNA_model_R = shared_drop_out(DNA_shared_MaxPool_1)(DNA_shared_Conv2D_1)(dna_r)
	combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
	DNA_LSTM = LSTM(32)(combined_DNA)
	
	combined_LSTM = DNA_LSTM
	
	combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
	
	output = Dense(2, activation = "softmax")(combined_dense)
	return output
	
myModel = ATAC_model(201)
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
	DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1)(dna_f)
	DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1)(dna_r)
	combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
	DNA_LSTM = LSTM(32)(combined_DNA)
	
	combined_LSTM = DNA_LSTM
	
	combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
	
	output = Dense(2, activation = "softmax")(combined_dense)
	return output
myModel = ATAC_model(201)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1)(dna_f)
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1)(dna_r)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1)(dna_f)
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1)(dna_r)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
DNA_shared_MaxPool_1(DNA_shared_Conv2D_1)(dna_f)
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
	DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
	DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
	combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
	DNA_LSTM = LSTM(32)(combined_DNA)
	
	combined_LSTM = DNA_LSTM
	
	combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
	
	output = Dense(2, activation = "softmax")(combined_dense)
	return output
myModel = ATAC_model(201)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
combined_DNA.shap-e
combined_DNA.shape
LSTM.__doc__
DNA_LSTM = LSTM(32)(combined_DNA)
DNA_LSTM = LSTM()(combined_DNA)
DNA_LSTM = LSTM(32,return_sequences=True)(combined_DNA)
combined_DNA.shape
model = Sequential()
model.add(Convolution1D(input_dim=4,
                        input_length=1000,
                        nb_filter=320,
                        filter_length=26,
                        border_mode="valid",
                        activation="relu",
                        subsample_length=1))
model.add(MaxPooling1D(pool_length=13, stride=13))
model.add(Dropout(0.2))
model = Sequential()
model.add(Convolution1D(input_dim=4,
                        input_length=1000,
                        nb_filter=320,
                        filter_length=26,
                        border_mode="valid",
                        activation="relu",
                        subsample_length=1))
from keras.layers.convolutional import Convolution1D, MaxPooling1D
model = Sequential()
model.add(Convolution1D(input_dim=4,
                        input_length=1000,
                        nb_filter=320,
                        filter_length=26,
                        border_mode="valid",
                        activation="relu",
                        subsample_length=1))
model.add(MaxPooling1D(pool_length=13, stride=13))
model.add(Dropout(0.2))
model.summary()
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
combined_DNA.summary()
combined_DNA.shape
DNA_model_F.shape
Lambda
from keras.layers import Layer, Input, Lambda
Flatten(combined_DNA).shape
test = Flatten(combined_DNA)
test = Flatten(DNA_model_R)
test = Flatten()(combined_DNA)
test.shape
DNA_model_F.shape
Flatten()(DNA_model_F).shape
DNA_shared_Conv2D_1(dna_f).shape
from keras.layers import Layer, Input, Lambda,TimeDistributed
TimeDistributed(combined_DNA)
TimeDistributed(combined_DNA).shape
TimeDistributed(combined_DNA)
TimeDistributed()(combined_DNA)
TimeDistributed(DNA_model_F)
TimeDistributed(DNA_model_F).shape
Flatten()(TimeDistributed(DNA_model_F)).shape
x = Conv2D(8, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1')(inp_layer)
x = Conv2D(16, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1_1')(x)
x = Conv2D(32, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1_2')(x)
input_shape = (300,300,1)
conv_filters = 16
kernel_size = (3, 3)
x = Conv2D(8, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1')(inp_layer)
x = Conv2D(16, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1_1')(x)
x = Conv2D(32, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1_2')(x)
x = Conv2D(8, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1')(inp_layer)
inp_layer = Input(name='the_input', shape=input_shape, dtype='float32',batch_shape=(1,300,300,1))
inp_layer.shape
dna_f.shape
x = Conv2D(8, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1')(inp_layer)
x = Conv2D(16, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1_1')(x)
x = Conv2D(32, kernel_size, padding='same',
                   activation='tanh', kernel_initializer='he_normal',
                    name='conv1_2')(x)
x.shape
from keras.layers import Reshape, Lambda
conv_to_LSTM_dims = (1,300,300,32)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
combined_DNA = keras.layers.concatenate([DNA_model_F, DNA_model_R], axis=-1)
combined_DNA.shape
DNA_model_F.shape
DNA_shared_MaxPool_1.shape
DNA_shared_Conv2D_1(dna_f).shape
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
DNA_model_F.shape
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(15, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(20, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
DNA_model_R.shape
dna_length - 15 - 20
dna_length - 15 - 20 + 1 + 1
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
conv_length=15
pool_length = 20
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
DNA_model_F.shape
DNA_LSTM_F = LSTM(20)(DNA_model_F)
DNA_LSTM_R = LSTM(20)(DNA_model_R)
combined_DNA = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
conv_length=15
pool_length = 20
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
DNA_LSTM_F = LSTM(20)(DNA_model_F)
DNA_LSTM_R = LSTM(20)(DNA_model_R)
combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
combined_LSTM = combined_DNA_LSTM
combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
output = Dense(2, activation = "softmax")(combined_dense)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
conv_length=15
pool_length = 20
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
DNA_LSTM_F = LSTM(20)(DNA_model_F)
DNA_LSTM_R = LSTM(20)(DNA_model_R)
combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
combined_LSTM = combined_DNA_LSTM
combined_dense = Dense(400, activation='relu')(Flatten(combined_LSTM))
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
conv_length=15
pool_length = 20
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
DNA_LSTM_F = LSTM(20)(DNA_model_F)
DNA_LSTM_R = LSTM(20)(DNA_model_R)
combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
combined_LSTM = combined_DNA_LSTM
combined_LSTM.shape
DNA_LSTM_F.shape
DNA_model_F.shape
combined_DNA_LSTM.shape
combined_dense = Dense(400, activation='relu')(combined_LSTM)
dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
conv_length=15
pool_length = 20
DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
DNA_LSTM_F = LSTM(20)(DNA_model_F)
DNA_LSTM_R = LSTM(20)(DNA_model_R)
combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
combined_LSTM = combined_DNA_LSTM
combined_dense = Dense(400, activation='relu')(combined_LSTM)
output = Dense(2, activation = "softmax")(combined_dense)
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
	conv_length=15
	pool_length = 20
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
	DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
	DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
	conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
	DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
	DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
	DNA_LSTM_F = LSTM(20)(DNA_model_F)
	DNA_LSTM_R = LSTM(20)(DNA_model_R)
	combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
	combined_LSTM = combined_DNA_LSTM
	combined_dense = Dense(400, activation='relu')(combined_LSTM)
	output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
	return output
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
	conv_length=15
	pool_length = 20
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
	DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
	DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
	conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
	DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
	DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
	DNA_LSTM_F = LSTM(20)(DNA_model_F)
	DNA_LSTM_R = LSTM(20)(DNA_model_R)
	combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
	combined_LSTM = combined_DNA_LSTM
	combined_dense = Dense(400, activation='relu')(combined_LSTM)
	output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
	return output
myModel = ATAC_model(201)
myModel.summary()
myModel.shape
seq = Model(outputs=myModel)
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
	conv_length=15
	pool_length = 20
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
	DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
	DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
	conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
	DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_F)
	DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='reshapeconvtolstm')(DNA_model_R)
	DNA_LSTM_F = LSTM(20)(DNA_model_F)
	DNA_LSTM_R = LSTM(20)(DNA_model_R)
	combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
	combined_LSTM = combined_DNA_LSTM
	combined_dense = Dense(400, activation='relu')(combined_LSTM)
	output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
	return [dna_f,dna_r],output
input,output = ATAC_model(201)
seq = Model(inputs=input, outputs=output)
def ATAC_model(dna_length):
	dna_f = Input(shape=(dna_length,4,1),name="DNA_Forward")
	dna_r = Input(shape=(dna_length,4,1),name="DNA_Reverse")
	Tn5_f = Input(shape=(dna_length,1,1),name="Tn5_Forward")
	Tn5_r = Input(shape=(dna_length,1,1),name="Tn5_Reverse")
	conv_length=15
	pool_length = 20
	DNA_shared_Conv2D_1 = Conv2D(20, kernel_size=(conv_length, 4),activation='relu',name="DNA_Conv2d_shared_1")
	DNA_shared_MaxPool_1 = MaxPooling2D(pool_size=(pool_length, 1),strides=(1,1),name="DNA_MaxPool_shared_1")
	DNA_model_F = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_f))
	DNA_model_R = DNA_shared_MaxPool_1(DNA_shared_Conv2D_1(dna_r))
	conv_to_LSTM_dims=(dna_length - conv_length - pool_length + 2,pool_length)
	DNA_model_F = Reshape(target_shape=conv_to_LSTM_dims, name='DNA_CNN_model_F')(DNA_model_F)
	DNA_model_R = Reshape(target_shape=conv_to_LSTM_dims, name='DNA_CNN_model_R')(DNA_model_R)
	DNA_LSTM_F = LSTM(20)(DNA_model_F)
	DNA_LSTM_R = LSTM(20)(DNA_model_R)
	combined_DNA_LSTM = keras.layers.concatenate([DNA_LSTM_F, DNA_LSTM_R], axis=-1)
	combined_LSTM = combined_DNA_LSTM
	combined_dense = Dense(400, activation='relu')(combined_LSTM)
	output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
	return [dna_f,dna_r],output
input,output = ATAC_model(201)
seq = Model(inputs=input, outputs=output)
seq.summary()
import statsmodels
df = pd.DataFrame()
import pandas as pd
df = pd.DataFrame()
df[0]=[1,2,3]
df[1]=[3,4,5]
df[2]=[3,4,5]
df[[0,1,2]].tolist()
df[[0,1,2]]
from pybedtools import BedTool
import pybedtools
import pandas as pd
pybedtools.BedTool("GUIDEseq_matched_subset_20191212.csv")
df = pd.read_csv("GUIDEseq_matched_subset_20191212.csv")
df.head()
df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].to_list()
df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].tolist()
df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].values.tolist()
test1 = df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].values.tolist()
bed = pybedtools.create_interval_from_list(test1)
test1 = df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].drop_na().values.tolist()
test1 = df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].dropna().values.tolist()
bed = pybedtools.create_interval_from_list(test1)
test1
df['Site_SubstitutionsOnly.Start'] = df['Site_SubstitutionsOnly.Start'].astype(int)
test1 = df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].dropna().astype(int,errors='ignore').values.tolist()
test1
bed = pybedtools.create_interval_from_list(test1)
bed = pybedtools.create_interval_from_list(['chr9', 133951317, 133951340])
test
bed
bed = pybedtools.create_interval_from_list(['chr9', 133951317, 133951340,'chr9', 133951317, 133951340])
bed
bed = pybedtools.create_interval_from_list([['chr19', 55116295, 55116318], ['chr1', 212161544, 212161567], ['chr19', 45600527, 45600550]])
bed_file1 = "tmp1"
df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End']].dropna().astype(int,errors='ignore').to_csv(bed_file1,sep="\t",header=False,index=False)
bed1 = pybedtools.Bedtool(bed_file1)
bed1 = pybedtools.BedTool(bed_file1)
bed1
bed1.head()
df[['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End','Site_SubstitutionsOnly.Sequence','sample','Site_SubstitutionsOnly.Strand']].dropna().astype(int,errors='ignore').to_csv(bed_file1,sep="\t",header=False,index=False)
bed1 = pybedtools.BedTool(bed_file1)
bed1.head()
df2 = pd.read_csv("CHANGEseq_GUIDEseq_classes_validation_set_20191212.csv")
df2.head()
bed_file2="tmp2"
df2.columns
df[['chr','start','end','sample']].dropna().astype(int,errors='ignore').to_csv(bed_file2,sep="\t",header=False,index=False)
df2[['chr','start','end','sample']].dropna().astype(int,errors='ignore').to_csv(bed_file2,sep="\t",header=False,index=False)
bed2 = pybedtools.BedTool(bed_file2)
bed2.head()
bed2.intersect(bed1)
bed3 = bed2.intersect(bed1)
bed3.head()
bed3 = bed2.intersect(bed1,wab=True)
bed3 = bed2.intersect(bed1,wa=True,wb=True)
bed3.head()
bed3 = bed2.intersect(bed1,wa=True,wb=True,d=True)
bed3 = bed2.intersect(bed1,wa=True,wao=True)
bed3 = bed2.intersect(bed1,wa=True,wo=True,wa=True)
bed3 = bed2.intersect(bed1,wa=True,wo=True)
bed3 = bed2.intersect(bed1,wao=True)
bed3.head()
bed3.columns()
bed3.columns
pybedtools.contrib.venn_maker.venn_maker(bed1+bed2)
pybedtools.contrib.venn_maker.venn_maker(bed1,bed2)
bed1 + bed2
(bed1 + bed2).head()
bed1.head()
bed3.head()
bed3.to_csv()
dir(bed3)
bed3.to_dataframe()
a=bed3.to_dataframe()
a.columns
bed3.head()
import readline
readline.write_history_file("pybedtools_example.py")
import pandas as pd
df = pd.DataFrame()
df[1]=['a','b']
df
df.index = df.index.astype(str)+","+df[1]
df
df.index = df.index.astype(str)+","+df[1]
df
df.index = df.index.astype(str)+","+df[1]
import pandas as pd
df = pd.read_csv("joined_CHANGEseq_GUIDEseq_classes_validation_set_updated_20191212.csv")
df2 = pd.read_csv("output.csv")
df.head()
df2.head()
df2 = pd.read_csv("output.csv",index_col=0)
df2.head()
df2['genome_coordinates']
df2.columns
df2['genome_coordinate']
df2.index = df2['genome_coordinate']
df2.index.unique()
df.index = df['genome_coordinate']
df
df2.head()
df['test_sub'] = df2['Site_SubstitutionsOnly.Sequence']
df.head()
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']]
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']]['Site_SubstitutionsOnly.Sequence','test_sub']
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']].columns
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']][['Site_SubstitutionsOnly.Sequence','test_sub']]
NaN != Nan
import readline
readline.write_history_file("xx.py")
GCCAAGGCACTGGATAGCTGAGG
import pandas as pd
df = pd.read_csv("joined_CHANGEseq_GUIDEseq_classes_validation_set_updated_20191212.csv")
df2 = pd.read_csv("output.csv")
df.head()
df2.head()
df2 = pd.read_csv("output.csv",index_col=0)
df2.head()
df2['genome_coordinates']
df2.columns
df2['genome_coordinate']
df2.index = df2['genome_coordinate']
df2.index.unique()
df.index = df['genome_coordinate']
df
df2.head()
df['test_sub'] = df2['Site_SubstitutionsOnly.Sequence']
df.head()
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']]
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']]['Site_SubstitutionsOnly.Sequence','test_sub']
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']].columns
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']][['Site_SubstitutionsOnly.Sequence','test_sub']]
df[df['Site_SubstitutionsOnly.Sequence']!=df['test_sub']][['Site_SubstitutionsOnly.Sequence','test_sub']].to_csv("check_sub.csv")
df.columns
df['test_gap'] = df2['Site_GapsAllowed.Sequence']
df.head()
df[df['Site_GapsAllowed.Sequence']!=df['test_gap']]
df[df['Site_GapsAllowed.Sequence']!=df['test_gap']]['Site_GapsAllowed.Sequence','test_gap']
df[df['Site_GapsAllowed.Sequence']!=df['test_gap']].columns
df[df['Site_GapsAllowed.Sequence']!=df['test_gap']][['Site_GapsAllowed.Sequence','test_gap']]
df[df['Site_GapsAllowed.Sequence']!=df['test_gap']][['Site_GapsAllowed.Sequence','test_gap']].to_csv("check_gap.csv")
import string
tab = string.maketrans("ACTG", "TGAC")
dir(string)
import str
import bytes
bytes
a = "Learning Tranlate() Methods"
print (a.translate(bytes.maketrans(b"aeiou", b"12345")))
import pandas as pd
df = pd.DataFrame()
df[1]=[5,6,7]
df.index[df[1].isin([5,6])]
df.index[df[1].isin([5,6])].tolist()
import string
tab = bytes.maketrans(b"ACTG", b"TGAC")
def revcomp(seq):
	tab = bytes.maketrans(b"ACTG", b"TGAC")
	return seq.translate(tab)[::-1]
revcomp(""ACGT)
revcomp("ACGT")
from keras.models import Model
import argparse
import matplotlib
matplotlib.use('agg')
import tensorflow as tf
import seaborn as sns
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
from keras.layers import Conv2D,LSTM,Conv1D
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
from janggu.data import Cover
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D,Input,GlobalAveragePooling2D,MaxPooling1D
from keras import backend as K
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
from sklearn.metrics import roc_curve,roc_auc_score
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
from keras.utils import plot_model
from copy import deepcopy as dp
from keras.layers import Layer, Input, Lambda,TimeDistributed
from keras.layers import Reshape, Lambda,Bidirectional
## Define parameters
high_hbf = 50
low_hbf = 0
input = "9638_A_include_old_score.csv"
flank = 200
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
bw_file="/home/yli11/Projects/Li_gRNA/footprint/H1_H2_GM12878_Tn5_bw/Hudep2.bw"
## read data
high,low = get_high_low_data(input,high_hbf,low_hbf)
roi_F,roi_R,roi = get_roi(high+low)
## get one-hot data and ATAC feature matrix
dna_F = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_F,flank=flank)
dna_R = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_R,flank=flank)
Tn5 = Cover.create_from_bigwig('bigwig_coverage',bigwigfiles=bw_file,roi=roi,binsize=1,stepsize=1,flank=flank)				
## ReShape
dna_F=np.reshape(dna_F,(len(high+low),flank*2+1,4))
dna_R=np.reshape(dna_R,(len(high+low),flank*2+1,4))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1))
high_hbf = 50
low_hbf = 0
input = "9638_A_include_old_score.csv"
flank = 200
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
bw_file="/home/yli11/Projects/Li_gRNA/footprint/H1_H2_GM12878_Tn5_bw/Hudep2.bw"
high,low = get_high_low_data(input,high_hbf,low_hbf)
roi_F,roi_R,roi = get_roi(high+low)
def get_high_low_data(input,pos_cutoff,neg_cutoff):
	df = pd.read_csv(input,index_col=0)
	pos = df[df['HbFBase']>=pos_cutoff].index.tolist()
	neg = df[df['HbFBase']<=neg_cutoff].index.tolist()
	return pos,neg
def get_roi(myList):
	# chr19:13180899-13180900
	chr = [x.split(":")[0] for x in myList]
	start = [int(x.split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x.split(":")[-1].split("-")[1]) for x in myList]
	roi_F = []
	roi_R = []
	roi = []
	for i in range(len(chr)):
		roi_F.append([chr[i],start[i],end[i],".",".","+"])
		roi_R.append([chr[i],start[i],end[i],".",".","-"])
		roi.append([chr[i],start[i],end[i]])
	return roi_F,roi_R,roi
high,low = get_high_low_data(input,high_hbf,low_hbf)
roi_F,roi_R,roi = get_roi(high+low)
dna_F = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_F,flank=flank)
dna_R = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_R,flank=flank)
Tn5 = Cover.create_from_bigwig('bigwig_coverage',bigwigfiles=bw_file,roi=roi,binsize=1,stepsize=1,flank=flank)				
dna_F=np.reshape(dna_F,(len(high+low),flank*2+1,4))
dna_R=np.reshape(dna_R,(len(high+low),flank*2+1,4))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1))
dna_F.shape
dna_F[0,200,:]
dna_F[:,200,:]
pd.DataFrame(dna_F[:,200,:])
pd.DataFrame(dna_F[:,200,:]).to_csv("test.csv")
dna_F[:,200,3]=0
dna_F[:,200,0]=1
pd.DataFrame(dna_F[:,200,:]).to_csv("test.csv")
flank=200
dna_F[[1,2],200,3]
dna_F[[1,2,3,4],200,3]
list(range(0,10,2))
from keras.models import Model
import argparse
import matplotlib
matplotlib.use('agg')
import tensorflow as tf
import seaborn as sns
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
from keras.layers import Conv2D,LSTM,Conv1D
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
from janggu.data import Cover
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D,Input,GlobalAveragePooling2D,MaxPooling1D
from keras import backend as K
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
from sklearn.metrics import roc_curve,roc_auc_score
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
from keras.utils import plot_model
from copy import deepcopy as dp
from keras.layers import Layer, Input, Lambda,TimeDistributed
from keras.layers import Reshape, Lambda,Bidirectional
import sklearn
sklearn.__version__
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
import imblearn; imblearn.show_versions()
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
from keras.models import Model
import argparse
import matplotlib
matplotlib.use('agg')
import tensorflow as tf
import seaborn as sns
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from janggu.data import Bioseq
from janggu.data import ReduceDim
import numpy as np
from keras.layers import Conv2D,LSTM,Conv1D
from keras.layers import AveragePooling2D
from janggu import inputlayer
from janggu import outputconv
from janggu import DnaConv2D
from janggu.data import ReduceDim
from janggu.data import Cover
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D,Input,GlobalAveragePooling2D,MaxPooling1D
from keras import backend as K
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
from sklearn.metrics import roc_curve,roc_auc_score
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
from keras.utils import plot_model
from copy import deepcopy as dp
from keras.layers import Layer, Input, Lambda,TimeDistributed
from keras.layers import Reshape, Lambda,Bidirectional
def auc(y_true, y_pred):
    auc = tf.metrics.auc(y_true, y_pred)[1]
    K.get_session().run(tf.local_variables_initializer())
    return auc
auc([0,0,1],[1,1,0])
def auc(y_true, y_pred):
    auc = tf.metrics.auc(y_true, y_pred)[1]
    # K.get_session().run(tf.local_variables_initializer())
    return auc
auc([0,0,1],[1,1,0])
auc([0,0,1],[1,1,0])[0]
auc([0,0,1],[1,1,0])[1]
auc([0,0,1],[1,1,0])
def get_roi(myList):
	# chr19:13180899-13180900+
	strand = [x[:-1] for x in myList]
	chr = [x[:-1].split(":")[0] for x in myList]
	start = [int(x[:-1].split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x[:-1].split(":")[-1].split("-")[1]) for x in myList]
	roi_A = []
	roi_T = []
	roi = []
	for i in range(len(chr)):
		roi_A.append([chr[i],start[i],end[i],".",".",strand[i]])
		roi_T.append([chr[i],start[i],end[i],".",".",get_opposite_strand(strand[i])])
		roi.append([chr[i],start[i],end[i]])
	return roi_A,roi_T,roi
	
def get_opposite_strand(x):
	if x == "+":
		return "-"
	if x == "-":
		return "+"
def get_high_low_data(input,pos_cutoff,neg_cutoff):
	df = pd.read_csv(input,index_col=0,sep="\t")
	df.index = df['coord']
	pos = df[df['HbFBase']>=pos_cutoff].index.tolist()
	neg = df[df['HbFBase']<=neg_cutoff].index.tolist()
	return pos,neg
## Define parameters
high_hbf = 50
low_hbf = 0
input = "9638_As_EBM_FDR.tsv"
flank = 200
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
bw_file="/home/yli11/Projects/Li_gRNA/footprint/H1_H2_GM12878_Tn5_bw/Hudep2.bw"
## read data
high,low = get_high_low_data(input,high_hbf,low_hbf)
roi_A,roi_T,roi = get_roi(high+low)
## get one-hot data and ATAC feature matrix
dna_A = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_A,flank=flank)
dna_T = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_T,flank=flank)
Tn5 = Cover.create_from_bigwig('bigwig_coverage',bigwigfiles=bw_file,roi=roi,binsize=1,stepsize=1,flank=flank)				
## ReShape
dna_A=np.reshape(dna_A,(len(high+low),flank*2+1,4))
dna_T=np.reshape(dna_T,(len(high+low),flank*2+1,4))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1))
from keras.callbacks import EarlyStopping
getattr
locals()
from bayes_opt import BayesianOptimization
keras.utils.to_categorical([0,0,1], 2)
keras.utils.to_categorical([1,1,1,0], 2)
5//2
4//2
6/2
6//2
def mutate_DNA(DNA_array,ref,alt,pos):
	out = dp(DNA_array)
	for i in out.shape[0]:
		if out[i,pos,ref]== 1:
			out[i,pos,ref] = 0
			out[i,pos,alt] = 1
		else:
			print ("Something wrong, existing...")
			exit()
	return outputconv(			
)
def mutate_DNA(DNA_array,ref,alt,pos):
	out = dp(DNA_array)
	for i in out.shape[0]:
		if out[i,pos,ref]== 1:
			out[i,pos,ref] = 0
			out[i,pos,alt] = 1
		else:
			print ("Something wrong, existing...")
			exit()
	return outputconv(
)
def mutate_DNA(DNA_array,ref,alt,pos):
	out = dp(DNA_array)
	for i in out.shape[0]:
		if out[i,pos,ref]== 1:
			out[i,pos,ref] = 0
			out[i,pos,alt] = 1
		else:
			print ("Something wrong, existing...")
			exit()
	return out
def get_pos_neg_data(DNA_array,Tn5_array,ref,alt,pos,high_low_labels):
	## only the mutated G in high HbF is positive
	X = []
	y = []
	Z = []
	## input data size here is double-ed, given ref alt
	for i in range(len(high_low_labels)):
		current_DNA = DNA_array[i,:,:]
		current_Tn5 = Tn5_array[i,:]
		if high_low_labels[i] == 1:
			## in high group, mutated G is positive and reference A is negative
			mutated_DNA = mutate_DNA(current_DNA,ref,alt,pos)
			X.append(mutated_DNA)
			X.append(current_DNA)
			y.append(1)
			y.append(0)
			Z.append(current_Tn5)
		if high_low_labels[i] == 0:
			## in low group, mutated G is negative and reference A is negative
			mutated_DNA = mutate_DNA(current_DNA,ref,alt,pos)
			X.append(mutated_DNA)
			X.append(current_DNA)
			y.append(0)
			Z.append(current_Tn5)
			Z.append(current_Tn5)			
	return np.concatenate(X),np.concatenate(Z),y
		
def cv_evaluation(dna_A,dna_T,Tn5,high,low,flank):
	## this is not a generic function
	## A G is hard coded
	outer = StratifiedKFold(n_splits=3)
	my_pred=[]
	my_true=[]
	y = np.array([1]*len(high)+[0]*len(low)) ## label for A, the true pos and neg label is different
	for train_index, test_index in outer.split(np.zeros(len(y)),y):
		dna_A_train, dna_A_test = dna_A[train_index,:,:], dna_A[test_index,:,:]
		dna_T_train, dna_T_test = dna_T[train_index,:,:], dna_T[test_index,:,:]
		Tn5_train, Tn5_test = Tn5[train_index,:], Tn5[test_index,:]		
		## get training testing data
		## one hot encoding: A=0, C=1, G=2, T=3
		X_A_train,y_A_train,Z_A_train = get_pos_neg_data(dna_A_train,Tn5_train,0,2,flank,y)
		X_A_test,y_A_test,Z_A_test = get_pos_neg_data(dna_A_test,Tn5_test,0,2,flank,y)
		X_T_train,y_T_train,Z_T_train = get_pos_neg_data(dna_T_train,Tn5_train,3,1,flank,y)
		X_T_test,y_T_test,Z_T_test = get_pos_neg_data(dna_T_test,Tn5_test,3,1,flank,y)	
		X_train = np.concatenate([X_A_train,X_T_train])
		Z_train = np.concatenate([Z_A_train,Z_T_train])
		y_train = y_A_train+y_T_train
		X_test = np.concatenate([X_A_test,X_T_test])
		Z_test = np.concatenate([Z_A_test,Z_T_test])
		## CNN model is a class, should have fit and predict method
		myModel = ABE_HbF_model()
		myModel.fit(X_train,Z_train,y_train)
		pred_y = myModel.predict(X_test,Z_test)
		my_pred += pred_y
		my_true += y[test_index]
	df = pd.DataFrame()
	df['true']=my_true
	df['pred']=my_pred
	return df
class ABE_HbF_model:
	"""Deep learning model for A base editor
	
	Attributes:
	
	"""
	def __init__(self,mode="design_1"):
		self.model = globals()[mode]()
		
	def fit(self,dna_train,Tn5_train,y_train):
		y_train = keras.utils.to_categorical(y_train, 2)
		self.model.compile(loss=keras.losses.binary_crossentropy,
						optimizer=keras.optimizers.Adadelta(),
						metrics=['accuracy',keras_auROC])
		es = EarlyStopping(monitor='keras_auROC', mode='max', verbose=1)
		self.model.fit([dna_train,Tn5_train],y_train,batch_size=32,epochs=50,verbose=1,validation_split=0.2,callbacks=[es])
		
	def predict(self,X_test,Z_test):
		my_pred = myModel.predict(X_test,Z_test)
		my_pred = [x[1] for x in my_pred]
		
		my_pred_A = my_pred[:len(my_pred)//2]
		my_pred_T = my_pred[len(my_pred)//2:]
		if len(my_pred_A) != len(my_pred_T):
			print ("length are not equal, existing")
			exit()
		mean_diff_list = []
		for i in range(0,len(my_pred_A),2):
			A_score = my_pred_A[i]
			T_score = my_pred_T[i]
			## mutated score
			G_score = my_pred_A[i+1]
			C_score = my_pred_T[i+1]
			mean_diff = np.mean([G_score-A_score,C_score-T_score])
			mean_diff_list.append(mean_diff)
		return mean_diff_list
		
def auc(y_true, y_pred):
    auc = tf.metrics.auc(y_true, y_pred)[1]
    K.get_session().run(tf.local_variables_initializer())
    return auc
	
def keras_auROC(y_true, y_pred):
	"""suppose each input is N*2, where [1] is the positive label"""
	return roc_auc_score([x[1] for x in y_true],[x[1] for x in y_pred])
def auc(y_true, y_pred):
    auc = tf.metrics.auc(y_true, y_pred)[1]
    K.get_session().run(tf.local_variables_initializer())
    return auc
def keras_auROC(y_true, y_pred):
	"""suppose each input is N*2, where [1] is the positive label"""
	return roc_auc_score([x[1] for x in y_true],[x[1] for x in y_pred])
dna_input = Input(shape=(dna_length,4),name="DNA_input")
conv_length_1=10
conv_length_2=15
conv_length_3=20
pool_length_DNA = 12
pool_length_Tn5 = 10
conv_length_Tn5 = 8
p = 0.8
DNA_shared_Conv1D_1 = Conv1D(20, kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_shared_Conv1D_2 = Conv1D(20, kernel_size=conv_length_2,activation='relu',name="DNA_Conv1d_2")(dna_input)
DNA_shared_Conv1D_3 = Conv1D(20, kernel_size=conv_length_3,activation='relu',name="DNA_Conv1d_3")(dna_input)
DNA_shared_MaxPool_1 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_1")(DNA_shared_Conv1D_1)
DNA_shared_MaxPool_2 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_2")(DNA_shared_Conv1D_2)
DNA_shared_MaxPool_3 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_3")(DNA_shared_Conv1D_3)
dna_length=201
dna_input = Input(shape=(dna_length,4),name="DNA_input")
conv_length_1=10
conv_length_2=15
conv_length_3=20
pool_length_DNA = 12
pool_length_Tn5 = 10
conv_length_Tn5 = 8
p = 0.8
DNA_shared_Conv1D_1 = Conv1D(20, kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_shared_Conv1D_2 = Conv1D(20, kernel_size=conv_length_2,activation='relu',name="DNA_Conv1d_2")(dna_input)
DNA_shared_Conv1D_3 = Conv1D(20, kernel_size=conv_length_3,activation='relu',name="DNA_Conv1d_3")(dna_input)
DNA_shared_MaxPool_1 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_1")(DNA_shared_Conv1D_1)
DNA_shared_MaxPool_2 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_2")(DNA_shared_Conv1D_2)
DNA_shared_MaxPool_3 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_3")(DNA_shared_Conv1D_3)
combined_DNA_pool = keras.layers.concatenate([Dropout(p)(DNA_shared_MaxPool_1),
											Dropout(p)(DNA_shared_MaxPool_2),
											Dropout(p)(DNA_shared_MaxPool_3)], axis=-1)
DNA_shared_MaxPool_1.shape
DNA_shared_Conv1D_1.shape
DNA_Conv1D_1 = Conv1D(20, kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_Conv1D_2 = Conv1D(20, kernel_size=conv_length_2,activation='relu',name="DNA_Conv1d_2")(dna_input)
DNA_Conv1D_3 = Conv1D(20, kernel_size=conv_length_3,activation='relu',name="DNA_Conv1d_3")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_3 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_3")(DNA_Conv1D_3)
DNA_MaxPool_1.shape
Flatten(DNA_MaxPool_1)
Flatten()(DNA_MaxPool_1)
Flatten()(DNA_MaxPool_1).shape
x = Input(batch_shape=(16, 10, 10))
print(x)
x = Input(batch_shape=(16, 10, 10))
x = Flatten()(x)
print(x)
x.shape
x = Input(shape = (12,100,10))
x = Dense(32)(x)
f = Flatten()(x)
f.shape
input = Input((35, 50))
output = Flatten()(input)
output.shape
DNA_Conv1D_1.shape
model = Sequential()
model.add(Conv2D(64, (3, 3),
                 input_shape=(3, 32, 32), padding='same',))
# now: model.output_shape == (None, 64, 32, 32)
model.add(Flatten())
model.shape
model
model.summary()
model = Sequential()
model.add(dna_input)
dna_input = Input(shape=(dna_length,4),name="DNA_input")
conv_length_1=10
conv_length_2=15
conv_length_3=20
pool_length_DNA = 12
pool_length_Tn5 = 10
conv_length_Tn5 = 8
p = 0.8
DNA_Conv1D_1 = Conv1D(20, kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_Conv1D_2 = Conv1D(20, kernel_size=conv_length_2,activation='relu',name="DNA_Conv1d_2")(dna_input)
DNA_Conv1D_3 = Conv1D(20, kernel_size=conv_length_3,activation='relu',name="DNA_Conv1d_3")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_3 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_3")(DNA_Conv1D_3)
combined_DNA_pool = keras.layers.concatenate([Flatten(DNA_MaxPool_1),
											Flatten(DNA_MaxPool_2),
											Flatten(DNA_MaxPool_3)], axis=-1)
dna_input = Input(shape=(dna_length,4),name="DNA_input")
conv_length_1=10
conv_length_2=15
conv_length_3=20
pool_length_DNA = 12
pool_length_Tn5 = 10
conv_length_Tn5 = 8
p = 0.8
DNA_Conv1D_1 = Conv1D(20, kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_Conv1D_2 = Conv1D(20, kernel_size=conv_length_2,activation='relu',name="DNA_Conv1d_2")(dna_input)
DNA_Conv1D_3 = Conv1D(20, kernel_size=conv_length_3,activation='relu',name="DNA_Conv1d_3")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_3 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_3")(DNA_Conv1D_3)
combined_DNA_pool = keras.layers.concatenate([Flatten()(DNA_MaxPool_1),
											Flatten()(DNA_MaxPool_2),
											Flatten()(DNA_MaxPool_3)], axis=-1)
combined_DNA_pool.shape
TN5_input = Input(shape=(dna_length),name="Tn5_input")
TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")
TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_model = TN5_MaxPool_1
TN5_input = Input(shape=(dna_length),name="Tn5_input")
TN5_input = Input(shape=(dna_length,1),name="Tn5_input")
TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")
TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_model = TN5_MaxPool_1
TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")
TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length_Tn5,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_input = Input(shape=(dna_length,1),name="Tn5_input")
TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")
TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length_Tn5,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
pool_length_Tn5
TN5_Conv1D_1.shape
TN5_input = Input(shape=(dna_length,1),name="Tn5_input")
TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")(TN5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length_Tn5,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_model = TN5_MaxPool_1
combined_DNA_Tn5 = keras.layers.concatenate([combined_DNA_pool,Flatten()(TN5_model)], axis=-1)	
combined_dense = Dense(100, activation='relu')(Dropout(p)(combined_dense))
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
combined_dense = Dense(1000, activation='relu')(Dropout(p)(combined_DNA_Tn5))
combined_dense = Dense(100, activation='relu')(Dropout(p)(combined_dense))
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")(Tn5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length_Tn5,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_model = TN5_MaxPool_1
# TN5_LSTM = Dropout(0.5)(Bidirectional(LSTM(32))(TN5_model))
combined_DNA_Tn5 = keras.layers.concatenate([combined_DNA_pool,Flatten()(TN5_model)], axis=-1)	
combined_dense = Dense(1000, activation='relu')(Dropout(p)(combined_DNA_Tn5))
combined_dense = Dense(100, activation='relu')(Dropout(p)(combined_dense))
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
model.summary()
def design_1(dna_length):
	dna_input = Input(shape=(dna_length,4),name="DNA_input")
	conv_length_1=10
	conv_length_2=15
	conv_length_3=20
	pool_length_DNA = 12
	pool_length_Tn5 = 10
	conv_length_Tn5 = 8
	p = 0.8
	DNA_Conv1D_1 = Conv1D(20, kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
	DNA_Conv1D_2 = Conv1D(20, kernel_size=conv_length_2,activation='relu',name="DNA_Conv1d_2")(dna_input)
	DNA_Conv1D_3 = Conv1D(20, kernel_size=conv_length_3,activation='relu',name="DNA_Conv1d_3")(dna_input)
	DNA_MaxPool_1 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_1")(DNA_Conv1D_1)
	DNA_MaxPool_2 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_2")(DNA_Conv1D_2)
	DNA_MaxPool_3 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_3")(DNA_Conv1D_3)
	combined_DNA_pool = keras.layers.concatenate([Flatten()(DNA_MaxPool_1),
												Flatten()(DNA_MaxPool_2),
												Flatten()(DNA_MaxPool_3)], axis=-1)
	Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
	TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")(Tn5_input)
	TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length_Tn5,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
	TN5_model = TN5_MaxPool_1
	combined_DNA_Tn5 = keras.layers.concatenate([combined_DNA_pool,Flatten()(TN5_model)], axis=-1)	
	combined_dense = Dense(1000, activation='relu')(Dropout(p)(combined_DNA_Tn5))
	combined_dense = Dense(100, activation='relu')(Dropout(p)(combined_dense))
	output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
	model = Model(inputs=[dna_input,Tn5_input], outputs=output)
	return model
def design_1(dna_length):
	dna_input = Input(shape=(dna_length,4),name="DNA_input")
	conv_length_1=10
	conv_length_2=15
	conv_length_3=20
	pool_length_DNA = 12
	pool_length_Tn5 = 10
	conv_length_Tn5 = 8
	p = 0.8
	DNA_Conv1D_1 = Conv1D(20, kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
	DNA_Conv1D_2 = Conv1D(20, kernel_size=conv_length_2,activation='relu',name="DNA_Conv1d_2")(dna_input)
	DNA_Conv1D_3 = Conv1D(20, kernel_size=conv_length_3,activation='relu',name="DNA_Conv1d_3")(dna_input)
	DNA_MaxPool_1 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_1")(DNA_Conv1D_1)
	DNA_MaxPool_2 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_2")(DNA_Conv1D_2)
	DNA_MaxPool_3 = MaxPooling1D(pool_size=pool_length_DNA,name="DNA_MaxPool_3")(DNA_Conv1D_3)
	combined_DNA_pool = keras.layers.concatenate([Flatten()(DNA_MaxPool_1),Flatten()(DNA_MaxPool_2),Flatten()(DNA_MaxPool_3)], axis=-1)
	Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
	TN5_Conv1D_1 = Conv1D(60, kernel_size=conv_length_Tn5,activation='relu',name="TN5_Conv1d")(Tn5_input)
	TN5_MaxPool_1 = MaxPooling1D(pool_size=pool_length_Tn5,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
	TN5_model = TN5_MaxPool_1
	combined_DNA_Tn5 = keras.layers.concatenate([combined_DNA_pool,Flatten()(TN5_model)], axis=-1)	
	combined_dense = Dense(1000, activation='relu')(Dropout(p)(combined_DNA_Tn5))
	combined_dense = Dense(100, activation='relu')(Dropout(p)(combined_dense))
	output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
	model = Model(inputs=[dna_input,Tn5_input], outputs=output)
	return model
true_predict_df = cv_evaluation(dna_A,dna_T,bw_values,high,low,flank)
dna_A=np.reshape(dna_A,(len(high+low),flank*2+1,4))
dna_T=np.reshape(dna_T,(len(high+low),flank*2+1,4))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1))
high_hbf = 50
low_hbf = 0
input = "9638_As_EBM_FDR.tsv"
flank = 200
refgenome="/home/yli11/Data/Human/hg19/fasta/hg19.fa"
bw_file="/home/yli11/Projects/Li_gRNA/footprint/H1_H2_GM12878_Tn5_bw/Hudep2.bw"
## read data
high,low = get_high_low_data(input,high_hbf,low_hbf)
roi_A,roi_T,roi = get_roi(high+low)
## get one-hot data and ATAC feature matrix
dna_A = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_A,flank=flank)
dna_T = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_T,flank=flank)
Tn5 = Cover.create_from_bigwig('bigwig_coverage',bigwigfiles=bw_file,roi=roi,binsize=1,stepsize=1,flank=flank)				
## ReShape
dna_A=np.reshape(dna_A,(len(high+low),flank*2+1,4))
dna_T=np.reshape(dna_T,(len(high+low),flank*2+1,4))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1))
roi_T
def get_roi(myList):
	# chr19:13180899-13180900+
	strand = [x[:-1] for x in myList]
	print (strand)
	chr = [x[:-1].split(":")[0] for x in myList]
	start = [int(x[:-1].split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x[:-1].split(":")[-1].split("-")[1]) for x in myList]
	roi_A = []
	roi_T = []
	roi = []
	for i in range(len(chr)):
		roi_A.append([chr[i],start[i],end[i],".",".",strand[i]])
		roi_T.append([chr[i],start[i],end[i],".",".",get_opposite_strand(strand[i])])
		roi.append([chr[i],start[i],end[i]])
	return roi_A,roi_T,roi
	
get_roi(["chr2:57948196-57948197-"])
def get_roi(myList):
	# chr19:13180899-13180900+
	strand = [list(x)[-1] for x in myList]
	print (strand)
	chr = [x[:-1].split(":")[0] for x in myList]
	start = [int(x[:-1].split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x[:-1].split(":")[-1].split("-")[1]) for x in myList]
	roi_A = []
	roi_T = []
	roi = []
	for i in range(len(chr)):
		roi_A.append([chr[i],start[i],end[i],".",".",strand[i]])
		roi_T.append([chr[i],start[i],end[i],".",".",get_opposite_strand(strand[i])])
		roi.append([chr[i],start[i],end[i]])
	return roi_A,roi_T,roi
	
get_roi(["chr2:57948196-57948197-"])
"chr2:57948196-57948197-"[:-1]
"chr2:57948196-57948197-"[-1]
myList=["chr2:57948196-57948197-"]
[x[-1] for x in myList]
def get_roi(myList):
	# chr19:13180899-13180900+
	# strand = [list(x)[-1] for x in myList]
	strand = [x[-1] for x in myList]
	print (strand)
	chr = [x[:-1].split(":")[0] for x in myList]
	start = [int(x[:-1].split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x[:-1].split(":")[-1].split("-")[1]) for x in myList]
	roi_A = []
	roi_T = []
	roi = []
	for i in range(len(chr)):
		roi_A.append([chr[i],start[i],end[i],".",".",strand[i]])
		roi_T.append([chr[i],start[i],end[i],".",".",get_opposite_strand(strand[i])])
		roi.append([chr[i],start[i],end[i]])
	return roi_A,roi_T,roi
	
get_roi(["chr2:57948196-57948197-"])
def get_roi(myList):
	# chr19:13180899-13180900+
	# strand = [list(x)[-1] for x in myList]
	strand = [x[-1] for x in myList]
	# print (strand)
	chr = [x[:-1].split(":")[0] for x in myList]
	start = [int(x[:-1].split(":")[-1].split("-")[0]) for x in myList]
	end = [int(x[:-1].split(":")[-1].split("-")[1]) for x in myList]
	roi_A = []
	roi_T = []
	roi = []
	for i in range(len(chr)):
		roi_A.append([chr[i],start[i],end[i],".",".",strand[i]])
		roi_T.append([chr[i],start[i],end[i],".",".",get_opposite_strand(strand[i])])
		roi.append([chr[i],start[i],end[i]])
	return roi_A,roi_T,roi
	
roi_A,roi_T,roi = get_roi(high+low)
roi_A
roi_T
## get one-hot data and ATAC feature matrix
dna_A = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_A,flank=flank)
dna_T = Bioseq.create_from_refgenome(name='dna',refgenome=refgenome,roi=roi_T,flank=flank)
Tn5 = Cover.create_from_bigwig('bigwig_coverage',bigwigfiles=bw_file,roi=roi,binsize=1,stepsize=1,flank=flank)				
## ReShape
dna_A=np.reshape(dna_A,(len(high+low),flank*2+1,4))
dna_T=np.reshape(dna_T,(len(high+low),flank*2+1,4))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1))
bw_values=np.reshape(Tn5,(len(high+low),flank*2+1,1))
true_predict_df = cv_evaluation(dna_A,dna_T,bw_values,high,low,flank)
def mutate_DNA(DNA_array,ref,alt,pos):
	out = dp(DNA_array)
	for i in range(out.shape[0]):
		if out[i,pos,ref]== 1:
			out[i,pos,ref] = 0
			out[i,pos,alt] = 1
		else:
			print ("Something wrong, existing...")
			exit()
	return out
true_predict_df = cv_evaluation(dna_A,dna_T,bw_values,high,low,flank)
def mutate_DNA(DNA_array,ref,alt,pos):
	out = dp(DNA_array)
	for i in range(out.shape[0]):
		print (i)
		if out[i,pos,ref]== 1:
			out[i,pos,ref] = 0
			out[i,pos,alt] = 1
		else:
			print ("Something wrong, existing...")
			exit()
	return out
true_predict_df = cv_evaluation(dna_A,dna_T,bw_values,high,low,flank)
def mutate_DNA(DNA_array,ref,alt,pos):
	out = dp(DNA_array)
	for i in range(out.shape[0]):
		print (i)
		print (out.shape)
		if out[i,pos,ref]== 1:
			out[i,pos,ref] = 0
			out[i,pos,alt] = 1
		else:
			print ("Something wrong, existing...")
			exit()
	return out
mutate_DNA(dna_A,0,2,200)
true_predict_df = cv_evaluation(dna_A,dna_T,bw_values,high,low,flank)
dna_A[1,:,:].shape
dna_A[0,:,:].shape
dna_A[0:0,:,:].shape
dna_A[0:1,:,:].shape
dna_A[0:1,1,1].shape
dna_A[0:1,1,1]
dna_A[0,1,1]
dna_A[1,1,1]
dna_A.shape
dna_A[1530:1531,1,1]
dna_A[1530:1532,1,1]
dna_A[1531:1532,1,1]
dna_A[1532:1532,1,1]
def get_pos_neg_data(DNA_array,Tn5_array,ref,alt,pos,high_low_labels):
	## only the mutated G in high HbF is positive
	X = []
	y = []
	Z = []
	## input data size here is double-ed, given ref alt
	for i in range(len(high_low_labels)):
		current_DNA = DNA_array[i:i+1,:,:]
		current_Tn5 = Tn5_array[i:i+1,:,:]
		if high_low_labels[i] == 1:
			## in high group, mutated G is positive and reference A is negative
			mutated_DNA = mutate_DNA(current_DNA,ref,alt,pos)
			X.append(mutated_DNA)
			X.append(current_DNA)
			y.append(1)
			y.append(0)
			Z.append(current_Tn5)
		if high_low_labels[i] == 0:
			## in low group, mutated G is negative and reference A is negative
			mutated_DNA = mutate_DNA(current_DNA,ref,alt,pos)
			X.append(mutated_DNA)
			X.append(current_DNA)
			y.append(0)
			Z.append(current_Tn5)
			Z.append(current_Tn5)			
	return np.concatenate(X),np.concatenate(Z),y
		
true_predict_df = cv_evaluation(dna_A,dna_T,bw_values,high,low,flank)
class ABE_HbF_model:
	"""Deep learning model for A base editor
	
	Attributes:
	
	"""
	def __init__(self,dna_length,mode="design_1"):
		self.model = globals()[mode](dna_length)
		print (self.model.summary())
	def fit(self,dna_train,Tn5_train,y_train):
		y_train = keras.utils.to_categorical(y_train, 2)
		self.model.compile(loss=keras.losses.binary_crossentropy,
						optimizer=keras.optimizers.Adadelta(),
						metrics=['accuracy',keras_auROC])
		es = EarlyStopping(monitor='keras_auROC', mode='max', verbose=1)
		self.model.fit([dna_train,Tn5_train],y_train,batch_size=32,epochs=50,verbose=1,validation_split=0.2,callbacks=[es])
		
	def predict(self,X_test,Z_test):
		my_pred = myModel.predict(X_test,Z_test)
		my_pred = [x[1] for x in my_pred]
		
		my_pred_A = my_pred[:len(my_pred)//2]
		my_pred_T = my_pred[len(my_pred)//2:]
		if len(my_pred_A) != len(my_pred_T):
			print ("length are not equal, existing")
			exit()
		mean_diff_list = []
		for i in range(0,len(my_pred_A),2):
			A_score = my_pred_A[i]
			T_score = my_pred_T[i]
			## mutated score
			G_score = my_pred_A[i+1]
			C_score = my_pred_T[i+1]
			mean_diff = np.mean([G_score-A_score,C_score-T_score])
			mean_diff_list.append(mean_diff)
		return mean_diff_list
		
def cv_evaluation(dna_A,dna_T,Tn5,high,low,flank):
	## this is not a generic function
	## A G is hard coded
	outer = StratifiedKFold(n_splits=3)
	my_pred=[]
	my_true=[]
	y = np.array([1]*len(high)+[0]*len(low)) ## label for A, the true pos and neg label is different
	for train_index, test_index in outer.split(np.zeros(len(y)),y):
		dna_A_train, dna_A_test = dna_A[train_index,:,:], dna_A[test_index,:,:]
		dna_T_train, dna_T_test = dna_T[train_index,:,:], dna_T[test_index,:,:]
		Tn5_train, Tn5_test = Tn5[train_index,:], Tn5[test_index,:]		
		## get training testing data
		## one hot encoding: A=0, C=1, G=2, T=3
		X_A_train,y_A_train,Z_A_train = get_pos_neg_data(dna_A_train,Tn5_train,0,2,flank,y)
		X_A_test,y_A_test,Z_A_test = get_pos_neg_data(dna_A_test,Tn5_test,0,2,flank,y)
		X_T_train,y_T_train,Z_T_train = get_pos_neg_data(dna_T_train,Tn5_train,3,1,flank,y)
		X_T_test,y_T_test,Z_T_test = get_pos_neg_data(dna_T_test,Tn5_test,3,1,flank,y)	
		X_train = np.concatenate([X_A_train,X_T_train])
		Z_train = np.concatenate([Z_A_train,Z_T_train])
		y_train = y_A_train+y_T_train
		X_test = np.concatenate([X_A_test,X_T_test])
		Z_test = np.concatenate([Z_A_test,Z_T_test])
		## CNN model is a class, should have fit and predict method
		myModel = ABE_HbF_model(2*flank+1)
		myModel.fit(X_train,Z_train,y_train)
		pred_y = myModel.predict(X_test,Z_test)
		my_pred += pred_y
		my_true += y[test_index]
	df = pd.DataFrame()
	df['true']=my_true
	df['pred']=my_pred
	return df
true_predict_df = cv_evaluation(dna_A,dna_T,bw_values,high,low,flank)
source
execfile
execfile()
exec(open("model1.py").read())
keras.utils.to_categorical([0,1,1], 2)
exec(open("model1.py").read())
X_A_train
exec(open("model1.py").read())
a=numpy.array([1,2,3])
b=numpy.array([4,5,6])
c=numpy.concatenate((a,b),axis=0)
a=np.array([1,2,3])
b=np.array([4,5,6])
c=np.concatenate((a,b),axis=0)
c.shape
a=np.array([1,2,3])
b=np.array([1,2,3])
c=np.concatenate([a,b],axis=0)
c.shape
a.shape
a=np.array([1,[2],[3]])
.shape
a.shape
exec(open("model1.py").read())
np.ones((1,2,2))
a=np.ones((1,2,2))
b=np.zeros((1,2,2))
c=np.concatenate([a,b,a,b])
c.shape
exec(open("model1.py").read())
a.shape
a=np.ones((1,2,2))
b=np.zeros((1,2,2))
c=np.concatenate([a,b,a,b,a,a,a])
c.shape
exec(open("model1.py").read())
df
true_predict_df
1100/(1100+150)
exec(open("model4_JASPAR.py").read())
high_hbf
exec(open("model4_JASPAR.py").read())
dna_length=201
motif_dict={}
parse_meme("motifs.meme",motif_dict,"")
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
conv_length_1=30
nb_filter_1 = 90
nb_filter_2 = 90
conv_length_2=20
pool_length = 10
p = 0.1
DNA_Conv1D = Conv1D(nb_filter_1, strides=3,kernel_size=conv_length_1,activation='relu',name="DNA_Conv1d")(dna_input)
TN5_Conv1D = Conv1D(nb_filter_1, strides=1,kernel_size=conv_length_1,activation='relu',name="TN5_Conv1d")(Tn5_input)
DNA_Conv1D.shape
dna_length-conv_length_1+1
(dna_length-conv_length_1+1)/3
ceil(dna_length-conv_length_1+1)/3
np.ceil(dna_length-conv_length_1+1)/3
np.ceil((dna_length-conv_length_1+1)/3)
exec(open("model4_JASPAR.py").read())
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
DNA_conv_length_1=26
DNA_conv_length_2=10
DNA_pool_length_1 = 13
DNA_pool_length_2 = 5
DNA_filter_1 = 100
DNA_filter_2 = 150
Tn5_filter_1 = 50
Tn5_conv_length = 8
Tn5_pool_length = 4
Dense_unit_1 = 1000
Dense_unit_2 = 100
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_Conv1D_2 = Conv1D(DNA_filter_2, kernel_size=DNA_conv_length_2,activation='relu',name="DNA_Conv1d_2")(DNA_MaxPool_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=DNA_pool_length_2,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_2 = Dropout(0.2)(DNA_MaxPool_2)
DNA_features =  Flatten()(DNA_MaxPool_2)
TN5_Conv1D_1 = Conv1D(Tn5_filter_1, kernel_size=Tn5_conv_length,activation='relu',name="TN5_Conv1d")(Tn5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=Tn5_pool_length,stride=1,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_MaxPool_1 = Dropout(0.2)(TN5_MaxPool_1)
Tn5_features =  Flatten()(TN5_MaxPool_1)
combined_DNA_Tn5 = keras.layers.concatenate([DNA_features,Tn5_features], axis=-1)	
combined_dense = Dense(Dense_unit_1, activation='relu')(combined_DNA_Tn5)
combined_dense = Dense(Dense_unit_2, activation='relu')(combined_dense)
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
dna_length=201
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
DNA_conv_length_1=26
DNA_conv_length_2=10
DNA_pool_length_1 = 13
DNA_pool_length_2 = 5
DNA_filter_1 = 100
DNA_filter_2 = 150
Tn5_filter_1 = 50
Tn5_conv_length = 8
Tn5_pool_length = 4
Dense_unit_1 = 1000
Dense_unit_2 = 100
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_Conv1D_2 = Conv1D(DNA_filter_2, kernel_size=DNA_conv_length_2,activation='relu',name="DNA_Conv1d_2")(DNA_MaxPool_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=DNA_pool_length_2,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_2 = Dropout(0.2)(DNA_MaxPool_2)
DNA_features =  Flatten()(DNA_MaxPool_2)
TN5_Conv1D_1 = Conv1D(Tn5_filter_1, kernel_size=Tn5_conv_length,activation='relu',name="TN5_Conv1d")(Tn5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=Tn5_pool_length,stride=1,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_MaxPool_1 = Dropout(0.2)(TN5_MaxPool_1)
Tn5_features =  Flatten()(TN5_MaxPool_1)
combined_DNA_Tn5 = keras.layers.concatenate([DNA_features,Tn5_features], axis=-1)	
combined_dense = Dense(Dense_unit_1, activation='relu')(combined_DNA_Tn5)
combined_dense = Dense(Dense_unit_2, activation='relu')(combined_dense)
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
DNA_conv_length_1=26
DNA_conv_length_2=10
DNA_pool_length_1 = 13
DNA_pool_length_2 = 5
DNA_filter_1 = 100
DNA_filter_2 = 150
Tn5_filter_1 = 50
Tn5_conv_length = 8
Tn5_pool_length = 4
Dense_unit_1 = 1000
Dense_unit_2 = 100
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_Conv1D_2 = Conv1D(DNA_filter_2, kernel_size=DNA_conv_length_2,activation='relu',name="DNA_Conv1d_2")(DNA_MaxPool_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=DNA_pool_length_2,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_2 = Dropout(0.2)(DNA_MaxPool_2)
DNA_features =  Flatten()(DNA_MaxPool_2)
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_Conv1D_2 = Conv1D(DNA_filter_2, kernel_size=DNA_conv_length_2,activation='relu',name="DNA_Conv1d_2")(DNA_MaxPool_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=DNA_pool_length_2,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_Conv1D_2.shape
DNA_MaxPool_1.shape
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
DNA_conv_length_1=26
DNA_conv_length_2=10
DNA_pool_length_1 = 13
DNA_pool_length_2 = 5
DNA_filter_1 = 100
DNA_filter_2 = 150
Tn5_filter_1 = 50
Tn5_conv_length = 8
Tn5_pool_length = 4
Dense_unit_1 = 1000
Dense_unit_2 = 100
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_features =  Flatten()(DNA_MaxPool_1)
TN5_Conv1D_1 = Conv1D(Tn5_filter_1, kernel_size=Tn5_conv_length,activation='relu',name="TN5_Conv1d")(Tn5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=Tn5_pool_length,stride=1,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_MaxPool_1 = Dropout(0.2)(TN5_MaxPool_1)
Tn5_features =  Flatten()(TN5_MaxPool_1)
combined_DNA_Tn5 = keras.layers.concatenate([DNA_features,Tn5_features], axis=-1)	
combined_dense = Dense(Dense_unit_1, activation='relu')(combined_DNA_Tn5)
combined_dense = Dense(Dense_unit_2, activation='relu')(combined_dense)
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
model.summary()
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
DNA_conv_length_1=26
DNA_conv_length_2=10
DNA_pool_length_1 = 13
DNA_pool_length_2 = 5
DNA_filter_1 = 100
DNA_filter_2 = 150
Tn5_filter_1 = 50
Tn5_conv_length = 8
Tn5_pool_length = 4
Dense_unit_1 = 1000
Dense_unit_2 = 100
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_features =  Flatten()(DNA_MaxPool_1)
TN5_Conv1D_1 = Conv1D(Tn5_filter_1, kernel_size=Tn5_conv_length,activation='relu',name="TN5_Conv1d")(Tn5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=Tn5_pool_length,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_MaxPool_1 = Dropout(0.2)(TN5_MaxPool_1)
Tn5_features =  Flatten()(TN5_MaxPool_1)
combined_DNA_Tn5 = keras.layers.concatenate([DNA_features,Tn5_features], axis=-1)	
combined_dense = Dense(Dense_unit_1, activation='relu')(combined_DNA_Tn5)
combined_dense = Dense(Dense_unit_2, activation='relu')(combined_dense)
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
model.summary()
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
DNA_conv_length_1=26
DNA_conv_length_2=10
DNA_pool_length_1 = 13
DNA_pool_length_2 = 5
DNA_filter_1 = 100
DNA_filter_2 = 150
Tn5_filter_1 = 50
Tn5_conv_length = 8
Tn5_pool_length = 4
Dense_unit_1 = 1000
Dense_unit_2 = 100
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,stride=3,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_Conv1D_2 = Conv1D(DNA_filter_2, kernel_size=DNA_conv_length_2,activation='relu',name="DNA_Conv1d_2")(DNA_MaxPool_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=DNA_pool_length_2,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_2 = Dropout(0.2)(DNA_MaxPool_2)
DNA_features =  Flatten()(DNA_MaxPool_2)
TN5_Conv1D_1 = Conv1D(Tn5_filter_1, kernel_size=Tn5_conv_length,activation='relu',name="TN5_Conv1d")(Tn5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=Tn5_pool_length,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_MaxPool_1 = Dropout(0.2)(TN5_MaxPool_1)
Tn5_features =  Flatten()(TN5_MaxPool_1)
combined_DNA_Tn5 = keras.layers.concatenate([DNA_features,Tn5_features], axis=-1)	
combined_dense = Dense(Dense_unit_1, activation='relu')(combined_DNA_Tn5)
combined_dense = Dense(Dense_unit_2, activation='relu')(combined_dense)
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
model.summary()
dna_input = Input(shape=(dna_length,4),name="DNA_input")
Tn5_input = Input(shape=(dna_length,1),name="Tn5_input")
DNA_conv_length_1=26
DNA_conv_length_2=10
DNA_pool_length_1 = 13
DNA_pool_length_2 = 5
DNA_filter_1 = 100
DNA_filter_2 = 150
Tn5_filter_1 = 50
Tn5_conv_length = 8
Tn5_pool_length = 4
Dense_unit_1 = 500
Dense_unit_2 = 100
DNA_Conv1D_1 = Conv1D(DNA_filter_1, kernel_size=DNA_conv_length_1,activation='relu',name="DNA_Conv1d_1")(dna_input)
DNA_MaxPool_1 = MaxPooling1D(pool_size=DNA_pool_length_1,stride=3,name="DNA_MaxPool_1")(DNA_Conv1D_1)
DNA_MaxPool_1 = Dropout(0.2)(DNA_MaxPool_1)
DNA_Conv1D_2 = Conv1D(DNA_filter_2, kernel_size=DNA_conv_length_2,activation='relu',name="DNA_Conv1d_2")(DNA_MaxPool_1)
DNA_MaxPool_2 = MaxPooling1D(pool_size=DNA_pool_length_2,name="DNA_MaxPool_2")(DNA_Conv1D_2)
DNA_MaxPool_2 = Dropout(0.2)(DNA_MaxPool_2)
DNA_features =  Flatten()(DNA_MaxPool_2)
TN5_Conv1D_1 = Conv1D(Tn5_filter_1, kernel_size=Tn5_conv_length,activation='relu',name="TN5_Conv1d")(Tn5_input)
TN5_MaxPool_1 = MaxPooling1D(pool_size=Tn5_pool_length,stride=3,name="TN5_MaxPool")(TN5_Conv1D_1)
TN5_MaxPool_1 = Dropout(0.2)(TN5_MaxPool_1)
Tn5_features =  Flatten()(TN5_MaxPool_1)
combined_DNA_Tn5 = keras.layers.concatenate([DNA_features,Tn5_features], axis=-1)	
combined_dense = Dense(Dense_unit_1, activation='relu')(combined_DNA_Tn5)
combined_dense = Dense(Dense_unit_2, activation='relu')(combined_dense)
output = Dense(2, activation = (tf.nn.softmax))(combined_dense)
model = Model(inputs=[dna_input,Tn5_input], outputs=output)
model.summary()
exec(open("model5.py").read())
high,low = get_high_low_data(input,high_hbf,low_hbf)
exec(open("model5.py").read())
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
X, y = load_iris(return_X_y=True)
clf = LogisticRegression(random_state=0).fit(X, y)
clf.get_params()
clf.density()
clf.densify()
clf.coef_
clf.coef_.shape
X.shape
Y.shape
y.shape
y
a=[1,2,3,4]
a[-2:]
exec(open("model5.py").read())
dna_A.shape
from Bio.Seq import Seq
from Bio import SeqIO
try:
	from StringIO import StringIO ## for Python 2
except ImportError:
	from io import StringIO ## for Python 3
import uuid
from joblib import Parallel, delayed
try:
	from StringIO import StringIO ## for Python 2
except ImportError:
	from io import StringIO ## for Python 3
def read_fasta(f):
	my_dict = {}
	for r in SeqIO.parse(f, "fasta"):
		my_dict[r.id] = str(r.seq).upper()
	return my_dict
	
def read_motif(meme_file):
	revcomp_file = uuid.uuid4()
	os.system("meme-get-motif -rc -all %s > /tmp/%s"%(meme_file,revcomp_file))
	original_motif_label = "++original++"
	revcomp_motif_label = "--revcomp--"
	
	dict1 = parse_meme(meme_file,label=original_motif_label)
	dict2 = parse_meme(revcomp_file,label=revcomp_motif_label)
	myDict = {}
	for k in dict1:
		motif_name = k.replace(original_motif_label,"")
		myDict[motif_name]=[dict1[k].T.values,dict2[k.replace(original_motif_label,revcomp_motif_label)].T.values]
	return myDict
	
exec(open("test_motif_scanning.py").read())
onehot_encoder
onehot_encoder.categories_
dna
dna=[['A',"G"],['A','C']]
onehot_encoder.fit_transform(dna)
dna
>>> enc = OneHotEncoder(handle_unknown='ignore')
>>> X = [['Male', 1], ['Female', 3], ['Female', 2]]
>>> enc.fit(X)
OneHotEncoder(handle_unknown='ignore')
>>> enc.categories_
[array(['Female', 'Male'], dtype=object), array([1, 2, 3], dtype=object)]
>>> enc.transform([['Female', 1], ['Male', 4]]).toarray()
array([[1., 0., 1., 0., 0.],
       [0., 1., 0., 0., 0.]])
>>> enc.inverse_transform([[0, 1, 1, 0, 0], [0, 0, 0, 1, 0]])
array([['Male', 1],
       [None, 2]], dtype=object)
enc = OneHotEncoder(handle_unknown='ignore')
enc.fit(dna)
enc.fit_transform(dna)
enc.fit_transform(dna).to_array()
enc.fit_transform(dna).toarray()
enc = OneHotEncoder(handle_unknown='ignore',categories=['A','C','G','T'])
enc.fit_transform(dna).toarray()
exec(open("test_motif_scanning.py").read())
dna_dict
myDNA_list.shape
np.stack(myDNA_list.shape)
np.stack(myDNA_list)
np.stack(myDNA_list).shape
exec(open("test_motif_scanning.py").read())
motifs
motifs['m1']
motifs['m1'][0]
motifs['m1'][0].shape
DNA_motif_scan(dna_A,motifs['m1'][0],motifs['m1'][1])
DNA_motif_scan(dna_A,motifs['m2'][0],motifs['m2'][1])
def DNA_motif_scan(DNA_array,m1,m2):
	score_list = []
	for i in range(DNA_array.shape[0]):
		score_list_1 = motif_scan(DNA_array[i,:,:],m1)
		score_list_2 = motif_scan(DNA_array[i,:,:],m2)
		# for j in range(len(score_list_1)):
		# 	if score_list_2[j] > score_list_1[j]:
		# 		score_list_1[j] = score_list_2[j]
		score_list.append(score_list_1)
	out = np.array(score_list)
	print ("DNA scanning out shape",out.shape)
	return out
	
DNA_motif_scan(dna_A,motifs['m2'][0],motifs['m2'][1])
dna_A
motifs['m2'][0]
>test1
ACGTGTGT
>test2
AAAAAAAA
>test4
CTCTACTC
>test1
ACGTGTGT
>test2
AAAAAAAA
>test4
CTCTACTC
>test1
ACGTGTGT
>test2
AAAAAAAA
>test4
CTCTACTC
ACGTGTGTmotifs['m2'][0]
motifs['m2'][0]
import numpy as np
a=np.nan
np.nanmean([a,2,5])
7/3.0
a = np.array([3,5,1,2,6])
a.argmax()
a.argmax(3)
a.argsort()
a.argsort(reverse=True)
a.argsort()[-2:]
a[a.argsort()[-2:]]
exec(open("feature_extraction.py").read())
from janggu.data import Bioseq
exec(open("feature_extraction.py").read())
1492/460
current_model
import sklearn
sklearn.__version__
import sklearn
sklearn .__version__
exec(open("main_classification.py").read())
plot_df
ddf
precision_recall_curve(ddf['true'],ddf['pred'])
a,b,c = precision_recall_curve(ddf['true'],ddf['pred'])
a
a[:4]
b[:4]
b
a
a[0]=0
a
a[:10]
b[:10]
sorted(a)
exec(open("main_classification.py").read())
ls
import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import string
import glob
import numpy as np
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
cadd = pd.read_csv("DeepSEA_CADD_GERP/gRNA_all_A.CADD.vcf",sep="\t")
deepsea = pd.read_csv("DeepSEA_CADD_GERP/DeepSEA.out.funsig",index_col=0)
gerp =pd.read_csv("DeepSEA_CADD_GERP/gRNA_all_GERP.tsv",sep="\t")
gerp.index = gerp['chrom']+":"+gerp['start'].astype(str)+"-"+gerp['end'].astype(str)
deepsea['name'] = deepsea['chr']+":"+(deepsea['pos']-1).astype(str)+"-"+deepsea['pos'].astype(str)
deepsea.index = deepsea['name']
cadd['name'] = "chr"+cadd['#Chrom'].astype(str)+":"+(cadd['Pos']-1).astype(str)+"-"+cadd['Pos'].astype(str)
cadd.index = cadd['name']
df = pd.read_csv("Editable_A_scores.tsv",sep="\t",index_col=0)
df['CADD'] = cadd['PHRED']
df['DeepSEA'] = deepsea['Functional significance score']
df['GERP'] = gerp['gerp_bp_score']
df['DeepSEA'] = df['DeepSEA'].apply(lambda x:-np.log10(x))
df.index = df.coord
df[['CADD','DeepSEA','GERP','HbFBase']].to_csv("Editable_A_scores.combined.scores.csv")
df = df.read_csv("Editable_A_scores.combined.scores.csv",index_col=0)
cadd = pd.read_csv("DeepSEA_CADD_GERP/gRNA_all_A.CADD.vcf",sep="\t")
deepsea = pd.read_csv("DeepSEA_CADD_GERP/DeepSEA.out.funsig",index_col=0)
gerp =pd.read_csv("DeepSEA_CADD_GERP/gRNA_all_GERP.tsv",sep="\t")
gerp.index = gerp['chrom']+":"+gerp['start'].astype(str)+"-"+gerp['end'].astype(str)
deepsea['name'] = deepsea['chr']+":"+(deepsea['pos']-1).astype(str)+"-"+deepsea['pos'].astype(str)
deepsea.index = deepsea['name']
cadd['name'] = "chr"+cadd['#Chrom'].astype(str)+":"+(cadd['Pos']-1).astype(str)+"-"+cadd['Pos'].astype(str)
cadd.index = cadd['name']
df = pd.read_csv("Editable_A_scores.tsv",sep="\t",index_col=0)
df['CADD'] = cadd['PHRED']
df['DeepSEA'] = deepsea['Functional significance score']
df['GERP'] = gerp['gerp_bp_score']
df['DeepSEA'] = df['DeepSEA'].apply(lambda x:-np.log10(x))
df.index = df.coord
df[['CADD','DeepSEA','GERP','HbFBase']].to_csv("Editable_A_scores.combined.scores.csv")
df = pd.read_csv("Editable_A_scores.combined.scores.csv",index_col=0)
df
df.shape
from decimal import Decimal
sns.set_style("whitegrid")
top_n = df[df['HbFBase']>=50]['DeepSEA'].tolist()
bot_n = df[df['HbFBase']==0]['DeepSEA'].tolist()
plot_df = pd.DataFrame([top_n,bot_n]).T
plot_df.columns = ['High',"Low"]
print (plot_df.describe())
plot_df = pd.melt(plot_df)
color_dict={}
color_dict['High'] = "#213fff"
color_dict['Low'] = "#6e899c"
sns.violinplot(x="variable",y='value',data=plot_df,palette =color_dict,linewidth=3,width=0.7,cut=3)
import matplotlib.pyplot as plt
y=5.2
h=0.3
print (scipy.stats.mannwhitneyu(top_n,bot_n).pvalue)
plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c="black")
plt.text(0.5, y+h+0.05, "Mann–Whitney U test: %.2E" % scipy.stats.mannwhitneyu(top_n,bot_n).pvalue, ha='center', va='bottom', color="black")
plt.ylim(-0.5,6)
plt.xticks([0,1],['High HbF score','Low HbF score'])
plt.xlabel("HbFBase scores")
plt.ylabel("DeepSEA scores (-log10)")
plt.savefig("DeepSEA-HbFBase-high-low.pdf", bbox_inches='tight')
import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import string
import glob
import numpy as np
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
df = pd.read_csv("Editable_A_scores.combined.scores.csv",index_col=0)
df['chr'] = [x[:-1].split(":")[0] for x in df.index]
df['start'] = [int(x[:-1].split(":")[-1].split("-")[0]) for x in df.index]
df['end'] = [int(x[:-1].split(":")[-1].split("-")[1]) for x in df.index]
print (df.head())
names = ['HbFBase','CADD',"DeepSEA"]
for n in names:  
    df['%s_rank'%(n)] = -df[n].rank(ascending=False)
print (df.head())
df = df.sort_values(['chr','start'])
print (df.head())
def to_matrix(df,n,top=500):
    c=n
    size=100
    tmp = df.copy()
    tmp = tmp.sort_values(n,ascending=False).head(n=top)
    my_index_list = tmp.index.tolist()
    out = []
    for i in my_index_list:
        line = []
        current_chr = df.at[i,"chr"]
        for j in range(-size,size+1):
            try:
                chr = df.at[i+j,"chr"]
            except:
                chr = "None"
            if chr == current_chr:
                value = df.at[i+j,"%s_rank"%(c)]
            else:
                value = np.nan
            line.append(value)
        out.append(line)
    out_df = pd.DataFrame(out)
    sel_cols = out_df.columns.tolist()
    out_df.index = my_index_list
    out_df['chr'] = tmp['chr']
    out_df['start'] = tmp['start']
    out_df['end'] = tmp['end']
    out_df['name'] = tmp['name']
    out_df['value'] = "."
    out_df['strand'] = "."
    print (df["%s_rank"%(c)].mean())
    out_df = out_df.fillna(df["%s_rank"%(c)].mean())
    out_df[['chr','start','end','name','value','strand']+sel_cols].to_csv("%s.computeMatrix.bed"%(c),header=False,index=False,sep="\t")
    return out_df[sel_cols] 
color_dict ={}
color_dict['HbFBase']='red'
color_dict['CADD']='green'
color_dict['DeepSEA']='blue'
fig, ax = plt.subplots()
for n in names:
    print (n)
    result_df = to_matrix(df,n)
    mean_line = pd.DataFrame(result_df.mean())
    test = pd.melt(result_df)
    sns.lineplot(x="variable", y="value", data=test,c=color_dict[n],ax=ax,label=n)
    ax.set_xticklabels(['']+list(range(-100,101,25)))
    ax.set_yticklabels(['']+list(range(5000,999,-1000))+[1])
    plt.ylabel("Average rank")
    plt.xlabel("Downstream / Upstream neighbors")   
color_dict ={}
color_dict['HbFBase']='red'
color_dict['CADD']='green'
color_dict['DeepSEA']='blue'
fig, ax = plt.subplots()
for n in names:
    print (n)
    result_df = to_matrix(df,n)
    mean_line = pd.DataFrame(result_df.mean())
    test = pd.melt(result_df)
    sns.lineplot(x="variable", y="value", data=test,c=color_dict[n],ax=ax,label=n)
    ax.set_xticklabels(['']+list(range(-100,101,25)))
    ax.set_yticklabels(['']+list(range(5000,999,-1000))+[1])
    plt.ylabel("Average rank")
    plt.xlabel("Downstream / Upstream neighbors")   
[x[:-1].split(":")[0] for x in df.index]
[x[:-1] for x in df.index]
import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import string
import glob
import numpy as np
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
df = pd.read_csv("Editable_A_scores.combined.scores.csv",index_col=0)
df['chr'] = [x[:-1].split(":")[0] for x in df.index]
df['start'] = [int(x[:-1].split(":")[-1].split("-")[0]) for x in df.index]
df['end'] = [int(x[:-1].split(":")[-1].split("-")[1]) for x in df.index]
df['name'] = [x[:-1] for x in df.index]
print (df.head())
df
pd.set_option('display.max_columns', None)
df
exec(open("peak_vs_neighbors.py").read())
df.head()
df = df.reset_index()
df = df.reset_index(drop=True)
df
exec(open("peak_vs_neighbors.py").read())
exec(open("feature_extraction.py").read())
data
exec(open("feature_extraction.py").read())
seq
test = pd.DataFrame.from_dict(seq)
test = pd.DataFrame.from_dict(seq,orient='index')
test
data
data['seq'] = test[0]
exec(open("feature_extraction.py").read())
data
test
data
data.shape
exec(open("feature_extraction.py").read())
data
data.shape
data.columns
data.describe()
exec(open("feature_extraction.py").read())
adjusted_scores
adjusted_scores.shape
adjusted_scores.index = high+low
adjusted_scores.shape
adjusted_scores
adjusted_scores.shape
data.shape
df = pd.concat([adjusted_scores,data],axis=1)
df.shape
df.isnull().any().any()
df.to_csv("ML_data.csv")
df
exec(open("main_classification.py").read())
df
exec(open("main_classification.py").read())
df
df.label
df[["label",'HbFBase']]
df[["label",'CADD',"DeepSEA"]].groupby('label').describe()
df[["label",'CADD',"DeepSEA"]].groupby('label').describe().T
Y
exec(open("main_classification.py").read())
import sklearn
sklearn.__version__
exec(open("main_classification.py").read())
exec(open("test.py").read())
df
exec(open("test.py").read())
exec(open("main_classification.py").read())
np.random([1,2,3])
exec(open("main_classification.py").read())
X_gkm
X_gkm.to_dict()
exec(open("main_classification.py").read())
Y
Y.index
(Y==1)
Y[Y==1]
exec(open("main_classification.py").read())
Y[Y==1]
Y[Y==1].index
Y[Y==1].index.tolist()
X_gkm.loc[Y[Y==1].index.tolist()].to_dict()
exec(open("main_classification.py").read())
a=pd.DataFrame()
a[0]=[1,2]
a[1]=[1,2]
a['a']=[1,2]
a
pd.melt(a)
a={}
a["a"]=[123]
a["a"]=[1,2,3]
a["b"]=[4,5,6]
a["c"]=[1,2,3]
pd.DataFrame(a)
exec(open("main_classification.py").read())
Y
seed = Y.sample(frac=0.1).index.tolist()
Y[Y.isin(seed)]
seed
Y[Y.index.isin(seed)]
Y[~Y.index.isin(seed)]
exec(open("main_classification.py").read())
np.random.randint(0,99999)
exec(open("main_classification.py").read())
for item in out_list:
	for k in X_dict:
		auROC_dict[k]+=item[k][1]
		auPRC_dict[k]+=item[k][2]
		df_list.append(item[k][0])
plot_df = pd.concat(df_list)
plot_auROC_multi(plot_df,color_dict)
plot_auPRC_multi(plot_df,color_dict)
boxplot_paired_t_test(pd.DataFrame(auROC_dict),color_dict,"Area under ROC curve","auROC_boxplot")
boxplot_paired_t_test(pd.DataFrame(auPRC_dict),color_dict,"Area under Precision-Recall curve","auPRC_boxplot")
exec(open("main_classification.py").read())
import pandas as pd
df = pd.DataFrame()
df[1]=[1,2,3]
df[2]=[3,4,5]
df.sum()
df.sum()>7
df[df.sum()>7]
df.loc[df.sum()>7]
df.loc[df.sum(axis=1)>7]
df[df.sum(axis=1)>7]
df.sum()
df.sum(axis=1)
df
df.columns[df.sum()>7]
df[df.columns[df.sum()>7]]
df.columns[df.sum()>7].tolist()
import pandas as pd
df = pd.read_csv("merged_rwu_GSE116177.tsv",sep="\t",index_col=0)
df.head()
df = df.T
df.head()
tmp = df.copy()
M = tmp.shape[1]
M
tmp[tmp<cutoff]=0
cutoff=1
tmp[tmp<cutoff]=0
tmp[tmp>=cutoff]=1
tmp.head()
tmp.sum().head()
(M/2.0-0.5)
import plotly.express as px
px.scatter
import pandas as pd
df = pd.read_csv("merged_rwu_GSE116177.tsv",sep="\t",index_col=0)
from sklearn import preprocessing
X_scaled = preprocessing.scale(df)
X_scaled
X_scaled.shape
df.shape
exec(open("d_score.py").read())
dna
exec(open("d_score.py").read())
df = pd.read_csv('TERT.GATK.hg19_multianno.vcf.gz',sep="\t",header=None,comment="#")
df[0] = "chr"+df[0].astype(str)
df = df[df[0]=="chr5"]
df = df[df[1].between(900907, 1668916, inclusive=True)]
print (df.shape)
df.to_csv("Genotyping.csv")
import pandas as pd
df = pd.read_csv("WTHUDEP2_GATA1.bw",sep="\t",header=None)
df[0].unique()
exec(open("main_classification.py").read())
X_RF
tmp_X
tmp_X = tmp_X.transform(lambda x:np.log2(x+1))
tmp_X
tmp_df = pd.concat([tmp_df,tmp_X],axis=1)
tmp_df.index= tmp_df['tmp_index'].tolist()
tmp_df = tmp_df.drop(['tmp_index'],axis=1)
X_RF = tmp_df
exec(open("main_classification.py").read())
from pytfmpval import tfmp
mat = (" 3  7  9  3 11 11 11  3  4  3  8  8  9  9 11  2"
        " 5  0  1  6  0  0  0  3  1  4  5  1  0  5  0  7"
        " 4  3  1  4  3  2  2  2  8  6  1  4  2  0  3  0"
        " 2  4  3  1  0  1  1  6  1  1  0  1  3  0  0  5"
       )
m = tfmp.read_matrix(mat)
tfmp.score2pval(m, 8.7737)
import pandas as pd
df = pd.read_csv("results.KO_vs_WT.txt",sep="\t",index_col=0)
df.head()
df = df[df['logFC'].abs()>=2]
df.shape
df = df[df['adj.P.Val']<=0.01]
df.shape
df = pd.read_csv("results.KO_vs_WT.txt",sep="\t",index_col=0)
df = df[df['logFC'].abs()>=1]
df.shape
df = df[df['adj.P.Val']<=0.01]
df.shape
df.to_csv("DEGs_logFC_1_FDR_2.csv")
df['gene'] = df.index.tolist()
df['gene'].to_csv("DEG.list",index=False,header=False)
import pandas as pd
name = pd.read_csv("mm9.ensembl_v67.t2g",sep="\t",header=None)
df = pd.read_csv("mm9.ensembl_v67.bed",sep="\t",header=None)
df.head()
name.head()
name = pd.read_csv("mm9.ensembl_v67.t2g",sep="\t")
name.head()
name.index = name['target_id']
df[3] = df[4].map(name['ext_gene'].to_dict())
df.head()
df[4]="."
df.head()
df.to_csv("mm9.ensembl_v67.gene_name.bed",sep="\t",header=False,index=False)
df.head()
df.index = df[3]
df.index = df[3].upper()
df.index = df[3].str.upper()
df.loc['ERRFI1']
a="12345"
a
a[3:4]
a[2:3]
a[2:3]="g"
b=list(a)
b[2]
b[2:3]
b[2:3]="g"
b
string(b)
b[2:4]="gc"
b
b[2:3]
b[2:3] == "g"
b[2:3] == list("g")
b[2:4] == list("fg")
b[2:4] = list("fg")
b
a="12345"
a[2]
a[2:3]
import numpy as np
a=[1,2,3]
a
b=[]
b.append(a)
b
c=[]
c.append(b)
c
x = np.array(c)
x
x.shape
np.concatenate([x,x])
np.concatenate([x,x]).shape
np.concatenate([x,x,x,x]).shape
c.append(b)
c
x = np.array(c)
x.shape
np.concatenate([x,x,x,x]).shape
np.stack([x,x])
np.stack([x,x]).shape
np.stack([x,x,x,x]).shape
def find_boundary(scores,candidate_start,candidate_end):
	max=0
	max_pos = -1
	for i in range(candidate_start,candidate_end+1):
		if scores[i] > max:
			max = scores[i]
			max_pos = i
	return max, max_pos
a=[4,2,1,5,10,3,1,4]
find_boundary(a,2,7)
27**(1/3)
27**(1/3.0)
27**(1/0)
a={}
a[1]=2
a
del a[1]
a
conda install -c bioconda pybigwig
import pandas as pd
df = pd.read_csv("Editable_A_scores.combined.scores.csv")
df.head()
def to_bed(x):
	chr,tmp = x.split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	strand = tmp[-1]
	return [chr,end,strand]
df['coord'].apply(to_bed).apply(pd.Series)
df['coord'].apply(to_bed)
def to_bed(x):
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	return [chr,end,strand]
df['coord'].apply(to_bed).apply(pd.Series)
df[['chr','pos','strand']] = df['coord'].apply(to_bed).apply(pd.Series)
df.head()
def to_bed(x):
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	if strand == "+":
		ref = "A"
		alt = "G"
	if strand == "-":
		ref = "T"
		alt = "C"		
	return [chr,end,ref,alt]
df[['chr','pos','ref','alt']] = df['coord'].apply(to_bed).apply(pd.Series)
df.head()
pos = df[df['HbFBase']>=50]
pos['id']='pos'
pos[['chr','pos','id','ref','alt']].to_csv("pos.vcf",sep="\t",header=False,index=False)
neg = df[df['HbFBase']==0]
neg['id']='pos'
neg[['chr','pos','id','ref','alt']].to_csv("neg.vcf",sep="\t",header=False,index=False)
pos.shape
neg.shape
import readline
readline.write_history_file("log.py")
import pandas as pd
def to_bed(x):
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	if strand == "+":
		ref = "A"
		alt = "G"
	if strand == "-":
		ref = "T"
		alt = "C"		
	return [chr,end,ref,alt]
def to_bed(x):
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	if strand == "+":
		ref = "A"
		alt = "G"
	if strand == "-":
		ref = "T"
		alt = "C"		
	return [chr,int(end),ref,alt]
df = pd.read_csv("Editable_A_scores.combined.scores.csv")
df[['chr','pos','ref','alt']] = df['coord'].apply(to_bed).apply(pd.Series)
df.head()
pos = df[df['HbFBase']>=50]
pos['id']='pos'
pos[['chr','pos','id','ref','alt']].to_csv("pos.vcf",sep="\t",header=False,index=False)
neg = df[df['HbFBase']==0]
neg['id']='pos'
neg[['chr','pos','id','ref','alt']].to_csv("neg.vcf",sep="\t",header=False,index=False)
df
def to_bed(x):
	print (x)
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	if strand == "+":
		ref = "A"
		alt = "G"
	if strand == "-":
		ref = "T"
		alt = "C"		
	return [chr,int(end),ref,alt]
df[['chr','pos','ref','alt']] = df['coord'].apply(to_bed).apply(pd.Series)
def to_bed(x):
	
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	if strand == "+":
		ref = "A"
		alt = "G"
	if strand == "-":
		ref = "T"
		alt = "C"		
	print (x,chr,int(end),ref,alt)
	return [chr,int(end),ref,alt]
df[['chr','pos','ref','alt']] = df['coord'].apply(to_bed).apply(pd.Series)
def to_bed(x):
	
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	if strand == "+":
		ref = "A"
		alt = "G"
	if strand == "-":
		ref = "T"
		alt = "C"		
	print (x,tmp,chr,int(end),ref,alt)
	return [chr,int(end),ref,alt]
df[['chr','pos','ref','alt']] = df['coord'].apply(to_bed).apply(pd.Series)
def to_bed(x):
	
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp
	if strand == "+":
		ref = "A"
		alt = "G"
	if strand == "-":
		ref = "T"
		alt = "C"		
	print (x,tmp,chr,int(end),ref,alt)
	return [chr,int(end),ref,alt]
df[['chr','pos','ref','alt']] = df['coord'].apply(to_bed).apply(pd.Series)
df.head()
pos = df[df['HbFBase']>=50]
pos['id']='pos'
pos[['chr','pos','id','ref','alt']].to_csv("pos.vcf",sep="\t",header=False,index=False)
neg = df[df['HbFBase']==0]
neg['id']='pos'
neg[['chr','pos','id','ref','alt']].to_csv("neg.vcf",sep="\t",header=False,index=False)
0.000929 < 1e-3
import networkx as nx
import pandas as pd
df = pd.read_csv("complete",sep="\t")
df.head()
df.columns
import pandas as pd
df = pd.read_csv("all_A_features.csv",index_col=0)
df2 = pd.read_csv("editable_A_all_bw_plus_HiChIP.csv",index_col=0)
df.index = [x[:-1] for x in df.index]
df.head()
df = pd.concat([df,df2],axis=1)
df.head()
df.shape
df.columns
df.isnull().any().any()
import readline
readline.write_history_file("combine_data.py")
(1+0.01)*30
(1+0.01)**30
(1+0.01)**15
1.1*10
1.1**10
1.1**3
1.1**5
1.1**8
1.08**12
1.04**12
1.05**12
1.06**12
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
# import umap
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
color_dict = {}
color_dict['RF:all']="#fa2111"
color_dict['RF:TFBS']="#fa8911"
color_dict['RF:Epi']="#fc42a5"
color_dict['DeepSEA']="#213fff"
color_dict['CADD']="#00bd3c"
color_dict['LS-GKM']="#464d4f"
glob.glob("*.csv")
glob.glob("*auROC*.csv")
glob.glob("*auROC*iter*.csv")
files = glob.glob("*auROC*iter*.csv")
df = pd.read_csv(files[0])
df.head()
pd.melt(df)
df  =pd.melt(df)
df.head()
glob.glob("*auROC*iter*.csv")
def parse_csv(x):
	df = pd.read_csv(x)
	df = pd.melt(df)
	df['resolution'] = x.split("_")[2]
	return df
df_list = [parse_csv(x) for x in files]
df = pd.concat(df_list)
def parse_csv(x):
	df = pd.read_csv(x)
	df = pd.melt(df)
	df['resolution'] = x.split("_")[2]
	return df
df_list = [parse_csv(x) for x in files]
df = pd.concat(df_list)
df.head()
df.shape
plt.figure()
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=median)
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median)
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
plt.figure()
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df)
ax.get_legend().remove()
plt.ylim(0.5,1)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
plt.figure()
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df)
plt.legend(loc='upper left',title="")
plt.ylim(0.5,1)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
plt.figure()
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["0","5","50","100","500"],hue_order=['RF:all','RF:TFBS','RF:Epi','LS-GKM','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,1)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
plt.figure()
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["0","5","50","100","500"],hue_order=['RF:all','RF:TFBS','RF:Epi','LS-GKM','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,0.85)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
plt.figure()
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["0","5","50","100","500"],hue_order=['RF:all','RF:TFBS','RF:Epi','LS-GKM','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,0.8)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
plt.figure(figsize=(12,5))
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["0","5","50","100","500"],hue_order=['RF:all','RF:TFBS','RF:Epi','LS-GKM','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,0.8)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
plt.figure(figsize=(12,5))
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["0","5","50","100","500"],hue_order=['RF:all','RF:TFBS','RF:Epi','LS-GKM','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,0.75)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
def parse_csv(x):
	df = pd.read_csv(x)
	df = pd.melt(df)
	df['resolution'] = "±"+x.split("_")[2]+"bp"
	return df
df_list = [parse_csv(x) for x in files]
df = pd.concat(df_list)
plt.figure(figsize=(12,5))
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["0","5","50","100","500"],hue_order=['RF:all','RF:TFBS','RF:Epi','LS-GKM','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,0.75)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
df_list = [parse_csv(x) for x in files]
df = pd.concat(df_list)
plt.figure(figsize=(12,5))
ax=sns.barplot(x='resolution',y='value',hue='variable',palette=color_dict,estimator=np.median,data=df,order=["±0bp","±5bp","±50bp","±100bp","±500bp"],hue_order=['RF:all','RF:TFBS','RF:Epi','LS-GKM','CADD','DeepSEA'])
plt.legend(loc='upper left',title="")
plt.ylim(0.5,0.75)
plt.savefig("barplot_auROC.pdf", bbox_inches='tight')
import readline
readline.write_history_file("barplot.py")
