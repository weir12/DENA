import sklearn
from sklearn.feature_selection import RFE
from sklearn.feature_selection import VarianceThreshold
from sklearn.svm import SVR
import textwrap
import argparse
import multiprocessing
import pandas as pd
import numpy as np
import h5py
import os,tqdm,pdb
from itertools import product
#from feature_selector import FeatureSelector
print(sklearn.__version__)
features=['base_q','length','mean','med','sd']
def get_chunk_size(M,N):
	assert all(type(x) is int for x in [M,N]) and M>N>0,"Please enter two positive integers and M>N"
	return np.array([(M//N)+1]*(M%N)+[M//N]*(N-(M%N)))
def single_file(file_fn,folder_fn,progressbar):
	with h5py.File(os.path.join(folder_fn,file_fn[0]),'r') as h5:
		res=np.vstack([h5[key][10-7:11+7] for key in features]).flatten()
		base="".join(["ATCG"[np.argmax(x)] for x in h5['motif'][:].reshape(4,-1).T[10-7:11+7]])
		progressbar.update(1)
	result=pd.Series(np.append(res,int(file_fn[0].split('_')[0])),dtype=object).append(pd.Series(base))
	return result	
def worker(df,folder_fn,pid,q):
	progressbar = tqdm.tqdm(total=df.shape[0], desc='progress '+str(pid), position=pid)		
	result=df.apply(single_file,args=(folder_fn,progressbar),axis=1)
	q.put(result)
	return None
def Extract(args):
	print("start extract features")
	all_files=pd.read_csv(args.s,header=None)
	chunks=get_chunk_size(all_files.shape[0],args.p)		

	with multiprocessing.Pool(processes = args.p) as pool:
		index=0
		pid=1
		q = multiprocessing.Manager().Queue()
		for i in chunks:
			w = pool.apply_async(worker,(all_files.iloc[index:index+i,:],args.i,pid,q))
			index+=i
			pid+=1
		pool.close()
		pool.join()
	header_df = pd.DataFrame(columns=[str(f)+"_"+str(p) for f in features for p in list(range(-7,8))]+["label"]+['motif'])
	header_df.to_csv(args.o,sep="\t",index=0)
	for i in range(q.qsize()) :
		values=q.get()
		values.to_csv(args.o ,mode='a',sep="\t",index=0,header=0)
	
if __name__ == "__main__":
	description = textwrap.dedent("""
		This script is a program to select best-n features of dataset 
		We recommend that you execute the following sub-commands in order:
			Extract    Extract RNN features from dataset
			Filter     Remove features with low variance
			RFE        SVM-FRE algorithm is used to select the best features
			lightGBM   
		See 'python ./main.py {sub-command} -h' to read about a options details.
		""")
	parser = argparse.ArgumentParser(
			description=description,
			formatter_class=argparse.RawDescriptionHelpFormatter)
	subparsers = parser.add_subparsers(title="Sub-command",dest='command')
	parser_a = subparsers.add_parser('Extract',formatter_class=argparse.RawDescriptionHelpFormatter,help='Extract RNN features from dataset')
	parser_a.add_argument('-s',required=True, default=None,
						help=("txt file contains names of valid dataset files"))
	parser_a.add_argument('-i',required=True,default=None,
						help=("Abs path of valid dataset(*.hdf5)"))
	parser_a.add_argument('-p',  default=32,
						help=("number of core used"))
	parser_a.add_argument('-o', default='./RNN_features.txt',
						help=("output file"))
	parser_a.set_defaults(func=Extract)
	
	args = parser.parse_args()
	args.func(args)

		
#	parser.add_argument('-m', '--model', default='SVM-RFE',
#                        help=("Algorithms used to select features"))

#	parser.add_argument('-v', '--VarianceThreshold',default=0.8,
#                       help=("Remove features with low variance or not(set 0)"))
	