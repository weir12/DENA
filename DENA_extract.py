# -*- coding: utf-8 -*-
# @Author: liangou
# @Date:   2021-02-25 12:39:08
# @Last Modified by:   liangou
# @Last Modified time: 2021-04-18 16:32:48

# -*- coding: utf-8 -*-
# @Author: liangou
# @Date:   2021-01-25 22:08:08
# @Last Modified by:   liangou
# @Last Modified time: 2021-02-22 20:53:50
import pandas as pd
import math
import numpy as np
from sklearn.utils import shuffle
import multiprocessing as mp
from multiprocessing import Process,Queue
import sys,os
from concurrent.futures import ThreadPoolExecutor
from tombo.tombo_helper import intervalData,get_raw_signal
from tombo import tombo_helper, tombo_stats, resquiggle
from tombo.tombo_helper import Fasta
from tombo.tombo_helper import TomboMotif
import torch.nn.functional as F
from joblib import Parallel, delayed
import textwrap
import argparse
from gettext import gettext
import pdb,re
from tqdm.auto import tqdm
import pysam
import time
from functools import partial
import logging 
import h5py
try:
	import bripy
	from bripy import *
except :
	sys.stderr.write("Warning: The BRI module could not be loaded")

import random
import threading

def get_raw_read_slot(fast5_data):
	raw_read_slot = next(iter(fast5_data['/Raw/Reads'].values()))
	return raw_read_slot

def find_event(mid_path,name):
	if "Events" in name:
		return mid_path+"/"+name

def get_single_slot_genome_centric(r_data, *slot_name):
	r_slot_values=[]
	mid_path='/'.join(('/Analyses', r_data.corr_group))
	find_e=partial(find_event, mid_path)
	h5=h5py.File(r_data.fn, 'r')
	mid_path=h5[mid_path].visit(find_e)
	for sigle_slot in slot_name:
		r_slot_values.append(h5[mid_path][sigle_slot])
	h5.close()
	return r_slot_values


 
'''

Multiple processes + coroutines was Used to optimize program execution speed
py ./tair_test.py predict --fast5 001"
  363  python ./tair_test.py predict --fast5 "/home/shihan/qinh_NCBI/elife_NDRS/col0_nanopore_drs/elif_col0_guppy324_bc/col0_guppy324_allpassf5_tombo/col0_all" 
--bam "/home/shihan/qinh_NCBI/elife_NDRS/col0_nanopore_drs/elif_col0_guppy324_bc_fq/trans_align/col0_drs_all_guppy324_bc_cdna.bam" --label VIRc
'''
def polish_signal(signal_list,padding_num=0,length=256):
	if len(signal_list)==length:
		return signal_list
	elif len(signal_list) < length:
		signal_list.extend([padding_num]*abs(len(signal_list)-length))
		return signal_list
	else:
		return signal_list[int((len(signal_list)//2)-(length/2)):int((len(signal_list)//2)+(length/2))]
def get_pos(args):
	'''
	
	'''
	fasta=Fasta(args.fasta)
	motif=TomboMotif(args.motif)
	#Cycle each chromosome/transcript
	with open(args.output,'w') as out_hl:
		for chrm in fasta.iter_chrms():
			chrm_seq=fasta.get_seq(chrm)
			for hit in motif.motif_pat.finditer(chrm_seq):
				 out_hl.write("\t".join([chrm,str(hit.start()),str(hit.end()),'+',hit.group()])+"\n")
def ont_hot(base_list):
	'''
	encode base symbol with one-hot 
	'''
	dict={"A":[1,0,0,0],"T":[0,1,0,0],"C":[0,0,1,0],"G":[0,0,0,1]}
	return np.array([dict[base.decode('UTF-8')] for base in base_list]).reshape(4,-1)
class ColoredArgParser(argparse.ArgumentParser):

	# color_dict is a class attribute, here we avoid compatibility
	# issues by attempting to override the __init__ method
	# RED : Error, GREEN : Okay, YELLOW : Warning, Blue: Help/Info 
	color_dict = {'RED' : '1;31', 'GREEN' : '1;32', 
				  'YELLOW' : '1;33', 'BLUE' : '1;36'}

	def print_usage(self, file = None):
		if file is None:
			file = sys.stdout
		self._print_message(self.format_usage()[0].upper() + 
							self.format_usage()[1:],
							file, self.color_dict['YELLOW'])

	def print_help(self, file = None):
		if file is None:
			file = sys.stdout
		self._print_message(self.format_help()[0].upper() +
							self.format_help()[1:],
							file, self.color_dict['BLUE'])

	def _print_message(self, message, file = None, color = None):
		if message:
			if file is None:
				file = sys.stderr
			# Print messages in bold, colored text if color is given.
			if color is None:
				file.write(message)
			else:
				# \x1b[ is the ANSI Control Sequence Introducer (CSI)
				file.write('\x1b[' + color + 'm' + message.strip() + '\x1b[0m\n')

	def exit(self, status = 0, message = None):
		if message:
			self._print_message(message, sys.stderr, self.color_dict['RED'])
		sys.exit(status)

	def error(self, message):
		self.print_usage(sys.stderr)
		args = {'prog' : self.prog, 'message': message}
		self.exit(2, gettext('%(prog)s: Error: %(message)s\n') % args)


def extract_features(read,start,end,bamfile,wins,debug=False,bri=None,header=None):
	sequ_fea=[x[start - read.start:end - read.start] for x in get_single_slot_genome_centric(read,*['length','norm_mean','norm_stdev'])]
	r_sig,seg,*_,scale =get_raw_signal(read,start,end)
	seg=seg[1:]-seg[0]
	r_sig=(r_sig-scale.shift)/scale.scale
	try:
		if isinstance(bamfile,str):
			sam_data = bri.get_alignments(read.read_id)
			read_bam=pysam.AlignedSegment.fromstring(sam_data,header)
		else:
			read_bam=bamfile[read.read_id]
		mapping_relation=pd.DataFrame(read_bam.get_aligned_pairs(matches_only=False),dtype=pd.Int64Dtype())
		m6a_pos_inread=mapping_relation[(mapping_relation[1]==start+wins[0])].index.tolist()[0]
		selected_base=mapping_relation.iloc[m6a_pos_inread-wins[0]:m6a_pos_inread+wins[1]+1,0].tolist()
		read_baseQ=[read_bam.query_qualities[x] if isinstance(x,np.int64) else 0 for x in selected_base]
		assert len(read_baseQ)==(end-start),"illegal length base-quality!"
		sequ_fea.insert(2,[np.median(x) for x in np.split(r_sig,seg)[:-1]])
		sequ_fea.insert(0,read_baseQ)
	except Exception as e:
		if debug:
			sys.stderr.write(str(e)+'\n')
		return None
	return read.read_id,",".join([str(x) for x in np.vstack(sequ_fea).flatten()])

def not_empty(s):
	return s

def get_logger():
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)
	rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
	log_path = os.path.dirname(os.getcwd()) + '/Logs/'
	if not os.path.exists(log_path):
		os.mkdir(log_path) 
	log_name = log_path + rq + '.log'
	logfile = log_name
	fh = logging.FileHandler(logfile, mode='w')
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return logger
#c7d8d026-1266-4e90-9537-42b82abd89e6

def predict_worker(region,args,p_num):
	reads_index = tombo_helper.TomboReads([args.fast5,],corrected_group=args.corr_grp)
	#If BRI-MODE, pass the file path, else write it into memory
	if not args.bri:
		bamfile={read.query_name:read for read in pysam.AlignmentFile(args.bam,"rb")}
		bri,header=None,None
	else:
		bamfile=args.bam
		bri = BamReadIndex(bamfile)
		header=pysam.AlignmentFile(bamfile, "rb").header


	tmp_fl=open("{0}_tmp".format(p_num),mode="w")
	for site in tqdm(list(region.itertuples()),position=p_num+1,desc="_".join(["processes",str(p_num)])):
		start,end=site.start,site.end
		reg_data = intervalData(chrm=site.ref,start=start,end=end,strand='+').add_reads(reads_index,require_full_span=True)
		reg_data.add_seq()
		reg_seq=reg_data.seq
		reg_data = intervalData(chrm=site.ref,start=start+(2-args.windows[0]),end=end+(args.windows[1]-2),strand='+')
		reg_data.add_reads(reads_index,require_full_span=True)
		if(len(reg_data.reads))<1:
			continue
		if args.debug:
			sys.stderr.write("{}\t{}found {:,} reads in fast5\n".format(site.ref,str(site.start)+"-"+str(site.end),len(reg_data.reads)))
		if(len(reg_data.reads))>10000:
			reg_data.reads=random.sample(reg_data.reads, 10000)
		site_fea=[]
		for read in reg_data.reads:
			site_fea.append(extract_features(read,start+(2-args.windows[0]),end+(args.windows[1]-2),bamfile,args.windows,args.debug,bri,header))
		site_fea=list(filter(not_empty,site_fea))
		if args.debug:
			sys.stderr.write("{}\t{}found {:,} reads in bamfile\n".format(site.ref,str(site.start)+"-"+str(site.end),len(site_fea)))
		if len(site_fea) >0:
			tmp_fl.write(">{}_{}_{}\n".format(site.ref,str(start+3),reg_seq))
			pd.DataFrame(site_fea).to_csv(tmp_fl,sep="\t",mode="a",header=0,index=0)
			tmp_fl.flush()  
	tmp_fl.close()
		#	q.put((site._asdict(),site_fea))
		 

def predict(args):

	'''
	fast5:a directory(has been re-squiggled by tombo) that contains the FAST5 files
	label:The label of the dataset You must choose one from 'a', 'm6a' and 'unknown"
	bam:BAM file used to extract base-quality(feature)
	sites:candidate position are used to extract features of mapped reads
	'''
	def split_df(df,N):
		'''
		Optimize multi-core usage
		Divide a df evenly into n pieces
		return a generator
		'''
		assert N >=1 and isinstance(N,int),"you should input a positive interge!"
		df = shuffle(df.values)
		return (i for i in np.array_split(site_df, N))
	if args.bri:
		status=os.system(" bri index "+args.bam)
		if status==0 :
			print("Bri index success") 
		else:
			args.bri=False
			print("Bri index fail,switch to normal mode")
	if args.debug:
		mylog=get_logger()
	site_df = pd.read_csv(args.sites,sep="\t",header=None,names=["ref", "start", "end", "strand","motif"])
	site_df = site_df[site_df['strand']=="+"] #We consider only motifs located in the plus-strand of the transcriptome
	if args.debug:
		mylog.debug('found {0} candidate sites'.format(site_df.shape[0]))
	region_batchs=split_df(site_df,args.processes)
	pool = mp.Pool(processes=args.processes, initargs=(mp.RLock(),), initializer=tqdm.set_lock)
	jobs =[]
	for p_num,region in enumerate(region_batchs):
		jobs.append(pool.apply_async(predict_worker, (region,args,p_num,)))
	pool.close()
	pool.join()
	for j in jobs:
		print(j.get())
	print("all tasks has done!")
def argument_parsing():
	description = textwrap.dedent("""
		=====================================================================================================
		using Neural network model to predict the RNA m6a modification status with single reads resolution
		We strongly recommend that you execute the following sub-commands in order:
			get_pos	  get candidate position from the reference fasta of the species 
			predict   dataset is feed into the model and obtain the predicted results   
		See 'python ./main.py {sub-command} -h' to read about a options details.
		author:https://github.com/weir12
		=====================================================================================================
		""")
	parser = ColoredArgParser(
			description=description,
			formatter_class=argparse.RawDescriptionHelpFormatter)

	subparsers = parser.add_subparsers(title="Sub-command",dest='command')
	parser_a = subparsers.add_parser('get_pos',formatter_class=argparse.RawDescriptionHelpFormatter,help='get candidate position')
	parser_a.add_argument('--fasta',required=True, default=None,
						help=("reference fasta"))
	parser_a.add_argument('--motif',  default='RRACH',
						help=("specifies a motif pattern"))
	parser_a.add_argument('--output', default='./candidate_predict_pos.txt',
						help=("output file"))
	parser_a.set_defaults(func=get_pos)
	parser_b = subparsers.add_parser('predict',formatter_class=argparse.RawDescriptionHelpFormatter,help='predict')
	parser_b.add_argument('--fast5',required=True, default=None,
						help=("a directory(has been re-squiggled by tombo) that contains the FAST5 files"))	
	parser_b.add_argument('--corr_grp',default="RawGenomeCorrected_000",
						help=("Analysis slot containing the re-squiggle information"))		
	parser_b.add_argument('--bam',required=True, default=None,	
						help=("BAM file used to extract base-quality(feature)"))
	parser_b.add_argument('--sites',default='./candidate_predict_pos.txt',	
						help=("candidate position are used to extract features of mapped reads"))
	parser_b.add_argument('--label',required=True,	
						help=("The string used to distinguish the sample"))	
	parser_b.add_argument('--windows',required=True,nargs=2,metavar='3',type=int,
						help=("Window drift away from the center of m6A"))
	parser_b.add_argument('--debug',action='store_true',default=False,
						help=("Enable debug mode (output more detailed run log)"))
	parser_b.add_argument('--bri',action='store_true',default=False,
						help=("Enable BRI mode (Reduce RAM consumption of BAM files)"))	
	parser_b.add_argument('--processes',default=8,type=int,
						help=("Number of processes allocated"))																		
	parser_b.set_defaults(func=predict)
	args = parser.parse_args()
	args.func(args)
if __name__ == '__main__':
	argument_parsing()
