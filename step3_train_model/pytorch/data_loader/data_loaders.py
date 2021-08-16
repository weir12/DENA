import data_set.data_set as datasets
from base import BaseDataLoader
class DENADataLoader(BaseDataLoader):
	"""
	DENA data loader
	hdf5_dir:folder(abs_path) contains total dataset (*.hdf5)
	file_list:list of all files in the hdf5_dir
	"""
	def __init__(self, hdf5_dir,file_list,dataset,batch_size=64, shuffle=True, validation_split=0.3, num_workers=24, training=True,features=None,windows_len=[10,10],sampler=None,balance_signal=False,nonsense_labels=False):
		self.hdf5_dir = hdf5_dir
		self.dataset = getattr(datasets,dataset)(self.hdf5_dir, file_list,features,windows_len,balance_signal,nonsense_labels)
		super().__init__(self.dataset,batch_size,shuffle,validation_split,num_workers,sampler)


class DENA_LSTM_DataLoader(BaseDataLoader):
	"""
	DENA LSTM data loader
	"""
	def __init__(self,motif,dataset,batch_size=64, shuffle=True, validation_split=0.3, num_workers=24, training=True,features=None,windows_len=[-2,2],sampler=None,nonsense_labels=False):

		self.dataset = getattr(datasets,dataset)(motif,windows_len)
		super().__init__(self.dataset,batch_size,shuffle,validation_split,num_workers,sampler)