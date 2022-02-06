# -*- coding: utf-8 -*-
# @Author: liangou
# @Date:   2021-02-26 11:14:00
# @Last Modified by:   liangou
# @Last Modified time: 2021-04-06 10:52:31


'''
1.filter and combination all the motif site 
2.apply model to features table
'''
import itertools
import pickle,subprocess
from collections import defaultdict
from collections import Counter
import numpy as np
import glob,os,pdb
import difflib
import argparse
import time
from pathlib import Path
import torch
import torch.nn as nn
#define const
def get_version():
    return 'DENA_3.3'
MOTIF="RRACH"
class LSTMModel(nn.Module):
	def __init__(self, in_dim, hidden_dim, n_layer, n_classes,drop_out=False):
		super(LSTMModel, self).__init__()
		self.lstm = nn.LSTM(in_dim, hidden_dim, n_layer, batch_first=True,bidirectional=True)
		self.classifier=nn.Sequential(
			nn.Linear(hidden_dim*2,hidden_dim*4),
			nn.ReLU6(inplace=True),
			nn.Linear(hidden_dim*4,2))

	def forward(self, x_event):
		out, (h_n, c_n) = self.lstm(x_event)
		x_event = out[:, -1, :]
		x = self.classifier(x_event)
		return x
tmp_fn="./predict_temp_{}".format(time.time())

SINGLE_LETTER_CODE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':['C','G','T'],
    'D':['A','G','T'], 'H':['A','C','T'], 'K':['G','T'], 'M':['A','C'],
    'N':['A','C','G','T'], 'R':['A','G'], 'S':['C','G'], 'V':['A','C','G'],
    'W':['A','T'], 'Y':['C','T']}


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v',  action='version', version=get_version(),help='Display version')
	parser.add_argument("-i", required=True,help = "Directory that contains the feature table")
	parser.add_argument("-m",  required=True,help = "Directory that contains the model file")
	parser.add_argument("-o",  required=True,help = "Output directory")
	parser.add_argument("-p",  required=True,help = "prefix of predicted result")
	parser.add_argument("-d",action='store_true',default=False,help = "Output the modification probability with read level")
	args = parser.parse_args()
	path=Path(args.m)
	models=[file for file in path.rglob("*.pth")]+[file for file in path.rglob("*.dat")]
	print(models)
	assert len(set([x.suffix for x in models]))==1,"The model files found seem to have different suffixes"
	print("We found {} well-trained model file".format(len(models)))
	models=[str(file) for file in models]
	inputs=glob.glob(os.path.join(args.i,"[0-9]*_tmp"))
	if not os.path.exists(args.o):
		os.makedirs(args.o)
		print("Creating output folder")
	f3=open(os.path.join(args.o,args.p+'.tsv'),'w')
	if args.d :
		f4=open(os.path.join(args.o,args.p+'_details.tsv'),'w')
	for pattern in itertools.product(*[ SINGLE_LETTER_CODE[base] for base in MOTIF]):
		pattern="".join(pattern)
		cmd=["awk","'/%s/{flag=1;print;next}/^>/{flag=0}flag'"%pattern,*inputs,">",tmp_fn]
		subprocess.check_call(" ".join(cmd),shell=True)
		with open(difflib.get_close_matches(pattern,models,1,1e-10)[0],'rb') as f1:
			if models[0].split(".")[-1]=="dat":
				model = pickle.loads(f1.read())
			else:
				model=LSTMModel(5,256,3,2)
				model_check=torch.load(f1,map_location='cpu')['state_dict']
				model.load_state_dict(model_check)
				model.eval()
		X=defaultdict(list)
		X2=defaultdict(list)
		with open(tmp_fn,'r') as f:
			for line in f.readlines():
				if line.startswith(">"):
					key=line.strip().strip(">")
				else:
					X[key].append(np.array(line.split()[1].split(","),dtype=np.float64))
					X2[key].append(line.split()[0])
		for site,body in X.items():
			if models[0].split(".")[-1]=="dat":
				predict_label=model.predict(np.array(body))
				predict_prob=model.predict_proba(np.array(body))[:,1]
				predict_label=predict_label.astype(float).astype('int')
				#pdb.set_trace()
				m6a_num=sum(predict_label==1)
			else:
				with torch.no_grad():
					RNN=torch.tensor(np.array(body).reshape(len(X2[site]),5,-1)).permute(0,2,1).to(torch.float32)
					predict_label = model(RNN)
					_, predicted = torch.max(predict_label, 1)
					predict_prob=predict_label[:,1].numpy()
					m6a_num=sum(predicted==1).item()
			f3.write("\t".join(site.split("_")+[str(x) for x in [m6a_num,predict_label.shape[0],m6a_num/predict_label.shape[0]]]))
			f3.write("\n")
			if args.d:
				f4.write("\t".join(site.split("_"))+"\n")
				pd.DataFrame([X2[site],predict_prob]).T.to_csv(f4,sep="\t",mode="a",header=0,index=0)
	f3.close()
	if args.d:
		f4.close()
	os.remove(tmp_fn)
if __name__ == '__main__':
	main()
