# -*- coding: utf-8 -*-
# @Author: liangou
# @Date:   2020-11-23 19:02:25
# @Last Modified by:   liangou
# @Last Modified time: 2020-11-23 19:41:54

'''
import modules
'''
import multiprocessing
import numpy as np
import pandas as pd
import os,sys,glob
import argparse
import joblib
import pdb,time
import datetime
import tqdm,contextlib


'''
define functions
'''

#function that resolves external(CMD) parameters
def my_argparse():
	parser = argparse.ArgumentParser(
	description="Label each read of the candidate site (based on the alignment error)",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--bam",help="indexed bam file",required=True)
	parser.add_argument("--label",type=int, choices=[1,0], help="The m6a reads(1) or non-m6a(0) class.",required=True)
	parser.add_argument("--sites_file",type=str,help="The file containing the site needed to be extracted",required=True)
	parser.add_argument("--outFolder", default='./feature', help="The default folder for outputing the results. Default: ./labels")
	parser.add_argument("--processes",type=int,default=64,help="number of processes for utilize")
	parser.add_argument("--plain_text_output",action='store_true',default=False,help="whether Output the generated labels in plain text format")	
	myargs=parser.parse_args()
	return myargs



if __name__ =='__main__':
	stime = time.time()
	myargs=my_argparse() 
	if not os.path.exists(myargs.outFolder):
		os.makedirs(myargs.outFolder)
	candidate_sites=pd.read_csv(myargs.sites_file,sep="\t",usecols=[1,3])
	pool = multiprocessing.Pool(processes = myargs.processes)
	q = multiprocessing.Manager().Queue()
	overall_progressbar = tqdm.tqdm(total=candidate_sites.shape[0], desc='extracting label...', position=0)
	def callback(result):
		overall_progressbar.update(1)
	for tup in zip(candidate_sites['chrom'], candidate_sites['error_geno_loci']):
		result=pool.apply_async(handle,(i,all_sites,myargs.windows,q,white_list,myargs.fast5_fn,myargs.uuid_limit,ambiguity,), callback=callback)	
