# -*- coding: utf-8 -*-
# @Author: liangou
# @Date:   2020-12-14 15:44:53
# @Last Modified by:   liangou
# @Last Modified time: 2020-12-16 23:01:42
#import modules
from __future__ import unicode_literals, absolute_import
from builtins import map
import sys
import argparse
import pandas as pd
import numpy as np
from tombo.tombo_helper import intervalData,get_raw_signal
import h5py
from tombo import tombo_helper, tombo_stats, resquiggle
import multiprocessing
import joblib,os
import pdb,time
import datetime
import tqdm,contextlib

if sys.version_info[0] > 2:
    unicode = str
class my_intervalData(intervalData):
	def __setattr__(self, name, value):
		self.__dict__[name] = value
	def __init__(self, chrm, start, end, strand=None, reg_id=None,\
				reg_text='', reads=None, seq=None,white_list=None):
		super().__init__(chrm, start, end, strand, reg_id,\
											reg_text, reads, seq)
		self.white_list=white_list
	def get_base_levels(self, output_path,label,site,group):
		def get_read_reg_events(r_data,start,end):
			r_means = tombo_helper.get_single_slot_genome_centric(r_data, 'norm_mean')
			r_stdev = tombo_helper.get_single_slot_genome_centric(r_data, 'norm_stdev')
			r_length = tombo_helper.get_single_slot_genome_centric(r_data, 'length')
			r_motif = tombo_helper.get_single_slot_genome_centric(r_data, 'base')
			r_sig,seg,*_,scale =get_raw_signal(r_data,start,end)
			seg=seg[1:]-seg[0]
			r_sig=(r_sig-scale.shift)/scale.scale
			r_med=[np.median(x) for x in np.split(r_sig,seg)[:-1]]
			if r_means is None:
				return None
			if r_data.start > self.start and r_data.end < self.end:
				# handle reads that are contained in a region
				# create region with nan values
				r_reg_means = np.full(self.end - self.start, np.NAN)
				r_reg_means[r_data.start - self.start:
							r_data.end - self.start] = r_means
			elif r_data.start > self.start:
				# handle reads that start in middle of region
				start_overlap = self.end - r_data.start
				# create region with nan values
				r_reg_means = np.full(self.end - self.start, np.NAN)
				r_reg_means[-start_overlap:] = r_means[:start_overlap]
			elif r_data.end < self.end:
				# handle reads that end inside region
				end_overlap = r_data.end - self.start
				# create region with nan values
				r_reg_means = np.full(self.end - self.start, np.NAN)
				r_reg_means[:end_overlap] = r_means[-end_overlap:]
			else:
				r_reg_means = r_means[
					self.start - r_data.start:self.end - r_data.start]
				r_reg_stdev = r_stdev[
					self.start - r_data.start:self.end - r_data.start]
				r_reg_length = r_length[
					self.start - r_data.start:self.end - r_data.start]
				r_reg_motif = r_motif[
					self.start - r_data.start:self.end - r_data.start]
				r_reg_motif=[x.decode('UTF-8') for x in r_reg_motif]
				#transform to one hot for a base sequence
				motif_onehot=[]
				for base in "ATCG":
					motif_onehot.extend([1 if x==base else 0 for x in r_reg_motif])
			return r_reg_means,r_med, r_reg_stdev, r_reg_length,motif_onehot,r_sig,seg

		if self.reads is None or len(self.reads) == 0:
			raise TomboError(
				'Must annotate region with reads ' +
				'(see `TomboInterval.add_reads`) to extract base levels.')


		for r_data in self.reads:
			if r_data.read_id not in self.white_list:
				continue
			if self.strand is not None and r_data.strand != self.strand:
				continue
			r_means, r_med,r_stdev, r_length,r_motif,r_sig,seg= get_read_reg_events(r_data,self.start,self.end)
			if r_means is None:
				continue
			with h5py.File(os.path.join(output_path,"_".join([str(label),str(site[0]),str(site[1]),r_data.read_id,".hdf5"])),'w') as h5:
				h5['mean']=r_means
				h5['sd']=r_stdev
				h5['med']=r_med
				h5['length']=r_length
				h5['motif']=r_motif
				h5['r_sig']=r_sig
				h5['seg']=seg
				tmp=[int(x) for x in group.loc[group[2]==r_data.read_id,3].str.split(",").values[0]]
				h5['base_q']=tmp



		return None
def docs():
	"""
	@Author: liang ou
	Brief Description:
		This script is used to extract electrical signal related features 
	Input file format requested:
	1.TAB delimited file without header and row-name(str):
		{chromosome}\t{position}\t{uuid}\t{extra columns}\n
			- {chromosome}:contig name,Must be consistent with the TOMBO Index
			- {position}:1-base coordinate
			- {uuid}:The UUID contained in the FAST5 file,NOT run_id or file name 
			- {extra columns}optional,retained in the output file without any processing)
	2.fast5_basedir(str):
		contians a set of reads and Tombo index file
	3.window_length(int):
		default:21
		determines the length of the window extending from the center point of M6A,to extract features
	4.features(list):
		default:all
		a set is picked from ['current_mean','current_std','current_median','num_samples','ont-hot_base','norm_singnal']
	5.output_directory:
		Path to save the output files
	If you have any questions,please discuss them in the issues section of https://github.com/weir12/DENA
	"""
def get_exter_args_parser():
	 parser = argparse.ArgumentParser(
	 	description='Extraction of current signal-related features from labeled reads of known loci  \
	 	with k-mer windows length',epilog="@Author: liang ou\t@e-mail:liangou@ips.ac.cn")
	 req_args = parser.add_argument_group('Required Arguments')
	 req_args.add_argument('--fast5s_basedir', required=True,metavar="./fast5",type=unicode,help="re-squiggled Directory containing fast5 files")
	 req_args.add_argument('--label', type=int, required=True,metavar="1",help="set 0 for negative dataset,1 for positive dataset")
	 req_args.add_argument('--site_files', type=str, required=True,metavar="./my_sites",help="A text file containing the extracted sites and reads_ids")
	 opt_args = parser.add_argument_group('Miscellaneous Arguments')
	 opt_args.add_argument('--processes',default=64,type=int,metavar="64",help="Number of CPU cores used")
	 opt_args.add_argument('--windows_length',default=[10,10],nargs = '+',type=int,metavar="[10,10]",help="The kmer length windows used to extract the feature")
	 opt_args.add_argument('--output_directory',default='./features',metavar="./output",type=unicode,help="Path to save the output files")
	 parser.add_argument('--print_readme', action='store_true',help="Print the detailed document and exit the program")
	 args = parser.parse_args()
	 print(args)
	 return args


def handle(reads_index,site,group,windows,outpath,label):
	with contextlib.redirect_stderr(open(os.devnull, "w")):
		reads_index = tombo_helper.TomboReads([reads_index,])
	transid=site[0]
	pos=site[1]
	reg_data = my_intervalData(chrm=transid,start=pos-1-windows[0],end=pos+windows[1],strand='+',white_list=group[2].values)
	reg_base_levels = reg_data.add_reads(reads_index,require_full_span=True).get_base_levels(outpath,label,site,group)


	return None
if __name__ == "__main__":	
	args=get_exter_args_parser()
	if args.print_readme is True:
		print(docs.__doc__)
		sys.exit(0)
	site_files=pd.read_csv(args.site_files,header=None,sep='\s+')
	if not os.path.exists(args.output_directory):
		os.makedirs(args.output_directory)
	pool = multiprocessing.Pool(processes = args.processes)

	overall_progressbar = tqdm.tqdm(total=site_files.drop_duplicates([0,1]).shape[0], desc='extracting current features...', position=0)
	def callback(result):
		overall_progressbar.update(1)
	for site,group in site_files.groupby([0,1]):
		result=pool.apply_async(handle,(args.fast5s_basedir,site,group,args.windows_length,args.output_directory,args.label), callback=callback)
	pool.close()
	pool.join()
	print("Done")