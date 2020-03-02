import os
import sys
import pandas

os.chdir("/data/joy/BBL/studies/alpraz/rawData/output/.heudiconv")
df = pandas.DataFrame()

for x in os.listdir("."):
	info_files = os.listdir(x+"/info/")
	info_files = [f for f in info_files if "dicominfo" in f]
	
	info_files = [os.getcwd() + '/' + x + "/info/" + f for f in info_files]
	
	for f in info_files:
		temp_df = pandas.read_csv(f, sep='\t')
		df = df.append(temp_df)

#a = list(df)
#print(a)
print(df.groupby(['image_type']).size())
