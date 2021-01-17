# -*- coding: utf-8 -*-
# @Author: liangou
# @Date:   2021-01-15 17:00:49
# @Last Modified by:   liangou
# @Last Modified time: 2021-01-17 15:46:19



###	import modules
import argparse
import multiprocessing
###

###	def func
def get_pos(fasta,finder,motif,output):
	'''
	a samll C++ program used to find positions in fasta which matching motif pattern 
	'''
	
if __name__ == "__main__":
	description = textwrap.dedent("""
		Neural network model was used to predict the m6a RNA modification status with single reads resolution
		We recommend that you execute the following sub-commands in order:
			get_pos	  get candidate position(based on motif_find) from the reference fasta(transcriptome is recommended) of the species 
			Extract   Extract dataset(feature+lables in hdf5 format) from a directory(has been re-squiggled by tombo) that contains the FAST5 files
			predict   dataset is feed into the model and obtain the predicted results   
		See 'python ./main.py {sub-command} -h' to read about a options details.
		""")
	parser = argparse.ArgumentParser(
			description=description,
			formatter_class=argparse.RawDescriptionHelpFormatter)
	subparsers = parser.add_subparsers(title="Sub-command",dest='command')
	parser_a = subparsers.add_parser('get_pos',formatter_class=argparse.RawDescriptionHelpFormatter,help='get candidate position')
	parser_a.add_argument('--fasta',required=True, default=None,
						help=("reference fasta"))
	parser_a.add_argument('--finder',required=True,default='./motiflocation/motiflocation',
						help=("path to the executable file after compiled of (source code [https://github.com/auton1/motiflocation])"))
	parser_a.add_argument('--motif',  default='RRACH',
						help=("specifies a motif pattern"))
	parser_a.add_argument('--output', default='./candidate_predict_pos.txt',
						help=("output file"))
	parser_a.set_defaults(func=get_pos)
	
	args = parser.parse_args()
	args.func(args)