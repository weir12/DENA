from torch.utils.data import Dataset
from torch import tensor
import pandas as pd
import h5py,os
import numpy as np
from utils.util import polish_signal

feature_table="/home/weir/m6a_model/DENA/step3_train_model/pytorch/feature_selection/RNN_features_motif.txt"
class DENA_dataset(Dataset):
	def __init__(self,hdf5_dir,file_list,features=False,windows_len=[7,7],balance_signal=False,nonsense_labels=False):
		super(DENA_dataset, self).__init__()
		self.index = pd.read_csv(file_list,header=None)
		self.windows_len=windows_len
		self.hdf5_dir=hdf5_dir
		self.features=['length', 'mean', 'med', 'motif', 'sd']		
		if nonsense_labels:
			self.target=np.random.choice([0,1],self.index.shape[0])
		else:
			self.target=self.index[0].str.split("_").str[0].to_numpy(dtype=np.int)	
	def __len__(self):
		return self.index.shape[0]
	def __getitem__(self,idx):
		with h5py.File(os.path.join(self.hdf5_dir,self.index.at[idx,0]),'r') as h5:
			RNN=tensor(np.vstack([h5[key][10-self.windows_len[0]:11+self.windows_len[1]] if key != 'motif' else h5['motif'][:].reshape(4,-1)[:,10-self.windows_len[0]:11+self.windows_len[1]] for key in self.features]).T).float()
			seg=np.insert(h5['seg'],0,0)
			signal=h5['r_sig']
			CNN=tensor([polish_signal(signal[seg[10-self.windows_len[0]:11+self.windows_len[1]]].tolist(),length=36*(sum(self.windows_len)+1))]) .float()
			target=tensor(self.target[idx])
			return RNN,CNN,target
			
			
class DENA_dataset_lstm(Dataset):
			"""
			motif="AGACT"
			"""
			def __init__(self, motif,windows=[-2,2],features=['base_q','length','mean','med','sd']):
				super(DENA_dataset_lstm, self).__init__()
				index = pd.read_csv(feature_table,sep="\t")
				colunms=[i+"_"+str(j) for i in features for j in range(windows[0],windows[1]+1)]
				
				colunms.append('label')
				self.index=index[index.motif.str.contains("[ATCG]{5}"+motif+"[ATCG]{5}")][colunms]
				self.index=self.index.reset_index(drop=True) 
				self.target=self.index['label']
				self.inter=windows[1]-windows[0]+1
				print(self.index)
				print(self.index.shape)
			def __len__(self):
				return self.index.shape[0]
			def __getitem__(self,idx):
				RNN=tensor(self.index.iloc[idx,:-1]).reshape(-1,self.inter).T.float()
				target=tensor(self.target[idx]).long()
				return RNN,target
			
		#return tensor(np.concatenate([np.expand_dims(x,1) for x in list(data[0:4]+[data[6]])]+[data[4].reshape(4,-1).T],axis=1)).float(),tensor([polish_signal(data[5].tolist(),length=768)]).float(),tensor(int(self.index.at[idx,0].split("_")[0]))