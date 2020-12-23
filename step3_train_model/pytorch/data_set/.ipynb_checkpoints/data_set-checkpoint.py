from torch.utils.data import Dataset
from torch import tensor
import pandas as pd
import h5py,os
import numpy as np
from utils.util import polish_signal
class DENA_dataset(Dataset):
	def __init__(self,hdf5_dir,file_list):
		super(DENA_dataset, self).__init__()
		self.index = pd.read_csv(file_list,header=None)
		self.hdf5_dir=hdf5_dir
	def __len__(self):
		return self.index.shape[0]
	def __getitem__(self,idx):
		try:
			with h5py.File(os.path.join(self.hdf5_dir,self.index.at[idx,0]),'r') as h5:
				data = [h5[key][:] for key in h5.keys()]
		except:
			#print("remove invalid file...",self.index.at[idx,0])
			print(self.index.at[idx,0])
			return tensor([1,2,3])
		return tensor(np.concatenate([np.expand_dims(x,1) for x in list(data[1:4]+[data[6]])]+[data[4].reshape(4,-1).T],axis=1)).float(),tensor([polish_signal(data[5].tolist(),length=768)]).float(),tensor(int(self.index.at[idx,0].split("_")[0]))