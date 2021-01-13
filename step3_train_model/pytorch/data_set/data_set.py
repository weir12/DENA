from torch.utils.data import Dataset
from torch import tensor
import pandas as pd
import h5py,os
import numpy as np
from utils.util import polish_signal
class DENA_dataset(Dataset):
	def __init__(self,hdf5_dir,file_list,features=False,windows_len=[10,10],balance_signal=False):
		super(DENA_dataset, self).__init__()
		uni_set={'sd','mean','med','motif','length'}
		msg="Please select a valid feature name from {}".format(" ".join(uni_set))
		self.index = pd.read_csv(file_list,header=None)
		self.windows_len=windows_len
		self.hdf5_dir=hdf5_dir
		self.features=set(features) if features else uni_set		
		assert self.features.issubset(uni_set),msg
	def __len__(self):
		return self.index.shape[0]
	def __getitem__(self,idx):
		with h5py.File(os.path.join(self.hdf5_dir,self.index.at[idx,0]),'r') as h5:
			RNN=tensor(np.vstack([h5[key][10-self.windows_len[0]:11+self.windows_len[1]] if key != 'motif' else h5['motif'][:].reshape(4,-1)[:,10-self.windows_len[0]:11+self.windows_len[1]] for key in self.features]).T).float()
			seg=np.insert(h5['seg'],0,0)
			signal=h5['r_sig']
			CNN=tensor([polish_signal(signal[seg[10-self.windows_len[0]:11+self.windows_len[1]]].tolist(),length=36*(sum(self.windows_len)+1))]) .float()
			target=tensor(int(self.index.at[idx,0].split("_")[0]))
			return RNN,CNN,target
			
			
			
			
			
		#return tensor(np.concatenate([np.expand_dims(x,1) for x in list(data[0:4]+[data[6]])]+[data[4].reshape(4,-1).T],axis=1)).float(),tensor([polish_signal(data[5].tolist(),length=768)]).float(),tensor(int(self.index.at[idx,0].split("_")[0]))