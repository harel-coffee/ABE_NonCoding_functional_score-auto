
import pandas as pd
import os
import sys
"""
for gRNA on the negative strand, need to make sure the editable bp position is in terms of positive strand.

strategy:

1. get editable window bed file

use getfasta to check if genome fasta is the same as the editable gRNA substring

2. merge bed file and get fasta

Force strandedness. That is, only merge features that are the same strand. By default, this is disabled.

3. count all As and Cs


workflow:

- 1. get an editable region bed file, and use `get_fasta` tp check if this editbale bed file is correct (the sequence should be the same as we substring gRNA)

- 2. get force strandness merged editable bed file, use `get_fasta` to get the sequence for the merged regions

- 3. check: all of our substring-ed gRNA should be a substring of this new sequence

- 4. count A and C in this sequence


"""


df = pd.read_csv("gRNA.all.bed",sep="\t",header=None)
start=int(sys.argv[1])
length=int(sys.argv[2])
print ("Editing Window starts at:%s, length is: %s"%(start,length))

def row_apply(x):
	chr = x[0]
	gRNA_start = x[1]
	gRNA_end = x[2]
	gRNA = x[3]
	value = x[4]
	strand = x[5]
	if strand == "-":
		edit_end = gRNA_end-start
		edit_start = gRNA_end-(start+length)
		edit_gRNA_seq = gRNA[start:(start+length)]
	if strand == "+":
		edit_start = gRNA_start+start
		edit_end = gRNA_start+start+length
		edit_gRNA_seq = gRNA[start:(start+length)]
	return [chr,edit_start,edit_end,edit_gRNA_seq,value,strand]
			
df = pd.DataFrame(df.apply(row_apply,axis=1).values.tolist())
df.to_csv("edit.bed",sep="\t",header=False,index=False)

get_fasta = "bedtools getfasta -fi /home/yli11/Data/Human/hg19/fasta/hg19.fa -bed edit.bed -fo edit.fa -tab -s -name"
os.system(get_fasta)
df = pd.read_csv("edit.fa",sep="\t",header=None)
df['flag'] = df[0].str.upper()==df[1].str.upper()
print ("The following output should be an empty dataframe")
print (df[df['flag']==False].head())

sort_command = "sort -k1,1 -k2,2n edit.bed > edit.bed.sorted"
merge = "bedtools merge -i edit.bed.sorted -s -c 4 -o collapse > edit.merged.bed"
os.system(sort_command)
os.system(merge)

command = """awk -F "\\t" '{print $1"\\t"$2"\\t"$3"\\t"$5"\\t0\\t"$4}' edit.merged.bed > edit.merged.bed6"""
os.system(command)
get_fasta = "bedtools getfasta -fi /home/yli11/Data/Human/hg19/fasta/hg19.fa -bed edit.merged.bed6 -fo edit.merged.fa -s -name -tab"
os.system(get_fasta)

## check
df = pd.read_csv("edit.merged.fa",sep="\t",header=None)
def row_apply(x):
	gRNAs = x[0].split(",")
	for g in gRNAs:
		if not g in x[1].upper():
			print ("some error")
			print (x[1])
	
df.apply(row_apply,axis=1)

editable_A = 0
editable_C = 0

for s in df[1]:
	for i in s.upper():
		if i == "A":
			editable_A += 1
		if i == "C":
			editable_C += 1

print ("editable_A:",editable_A,"editable_C:",editable_C,"Total:",editable_A+editable_C)

os.system("rm edit*")

