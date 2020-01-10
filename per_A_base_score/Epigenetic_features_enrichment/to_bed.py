def to_bed(x):
	strand = x[-1]
	chr,tmp = x[:-1].split(":")
	start,tmp = tmp.split("-")
	end = tmp[:-1]
	return [chr,end,strand]


df[['chr','pos','strand']] = df['coord'].apply(to_bed).apply(pd.Series)


