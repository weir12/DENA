from pandarallel import pandarallel
pandarallel.initialize(nb_workers=32,progress_bar=True,use_memory_fs=False)
import pandas as pd
import h5py,os
def get_exter_args_parser():
	 parser = argparse.ArgumentParser(
	 	description='Check the integrity of the HDF5 file and output a usable file name',epilog="@Author: liang ou\t@e-mail:liangou@ips.ac.cn")
	 req_args = parser.add_argument_group('Required Arguments')
	 req_args.add_argument('--hdf5_basedir', required=True,metavar="./fast5",type=unicode,help="re-squiggled Directory containing fast5 files")
	 req_args.add_argument('--file_list', type=str, required=True,metavar="./my_sites",help="A text file containing hdf5_file_names need to checked")
	 opt_args = parser.add_argument_group('Miscellaneous Arguments')
	 opt_args.add_argument('--processes',default=32,type=int,metavar="64",help="Number of CPU cores used")
	 opt_args.add_argument('--progress_bar',default=True,type=bool,help="display progress bar")
	 opt_args.add_argument('--output_prefix',default='filterd_',metavar="filterd",type=unicode,help="The prefix of the output file")
	 args = parser.parse_args()
	 print(args)
	 return args
def check_hdf5(file_name,folder):
	try:
		with h5py.File(os.path.join(folder,file_name[0]),'r') as h5:
			data = [h5[key][:] for key in h5.keys()]
			return file_name
	except:
		return None



if __name__ == "__main__":	
	args=get_exter_args_parser()
	pandarallel.initialize(nb_workers=args.processes,progress_bar=args.progress_bar)
	rawdata=pd.read_csv(args.file_list,header=None)
	rawdata.apply(check_hdf5,args=(args.hdf5_basedir,),axis=1)
	
