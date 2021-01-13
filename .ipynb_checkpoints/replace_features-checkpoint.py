import h5py
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True,use_memory_fs=False)
def handle(x):
	file_name=x[0]
	new_data=[int(y) for y in x[1].split(",")]
	try:
		with h5py.File("./features/"+file_name,"r+") as file:
			del file['base_q']  
			file['base_q'] = new_data
		return file_name
	except:
		return None




def main():
	rawdata=pd.read_csv("new_features",sep="\t",index_col=0)
	result=rawdata.parallel_apply(handle,axis=1)













if __name__ == '__main__':
	main()