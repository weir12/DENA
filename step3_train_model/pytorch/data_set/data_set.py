from torch.utils.data import Dataset
import pandas as pd
import h5py,os
class DENA_dataset(Dataset):
	def __init__(self,hdf5_dir,file_list):
		super(DENA_dataset, self).__init__()
		self.index = pd.read_csv(file_list)
		self.hdf5_dir=hdf5_dir
	def __len__(self):
		return self.index.shape[0]
	def __gettime__(self,idx):
		with h5py.File(os.path.join(hdf5_dir,index.at[idx,0]),'r') as h5:
			data = [h5[key][:] for key in h5.keys()]
		return data