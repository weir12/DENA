# -*- coding: utf-8 -*-
# @Author: liangou
# @Date:   2020-11-23 19:02:25
# @Last Modified by:   liangou
# @Last Modified time: 2020-12-09 14:53:26

'''
import modules
'''
import multiprocessing
import numpy as np
import pandas as pd
import pysam
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
#The features of the alignment level of reads were extracted
#base-quality & mapping-quanlyt
def extract_align_features(reads,pos,win_len=10):
	result=[]
	uuids=[]
	for read in reads:
		try:
			quality=read.query_qualities #base-quality of read
			mapping_relation=pd.DataFrame(read.get_aligned_pairs(matches_only=False),dtype=pd.Int64Dtype()) #mapping_relation between read(query) and reference 
			m6a_pos_inread=mapping_relation[(mapping_relation[1]==int(pos))].index.tolist()[0]
			selected_base=mapping_relation.iloc[m6a_pos_inread-win_len:m6a_pos_inread+win_len+1,0].tolist()
			result.append([quality[x] if isinstance(x,np.int64) else 0 for x in selected_base])
			uuids.append(read.query_name)
		except:
			pass
	return uuids,result
#Reads were classified into positive and negative classes according to the characteristics of the align results
def handle(site,group,q,bamdir,label):
	bamfile = pysam.AlignmentFile(bamdir, "rb")
	error_reads=[]
	match_reads=[]
	for i,r in group.iterrows():
		iter = bamfile.fetch(r[0], r[1], r[2])
		for read in iter:
			if read.is_supplementary or read.is_secondary or read.is_unmapped :
				continue
			if r[1] in read.get_reference_positions():
				match_reads.append(read)
			else:
				error_reads.append(read)

	positive_reads=set(error_reads)
	negative_reads=set(match_reads)-set(error_reads)
	if label==0:
		q.put({site:extract_align_features(negative_reads,site[1],10)})
	q.put({site:extract_align_features(positive_reads,site[1],10)})
	bamfile.close()

'''
main func
'''
if __name__ =='__main__':
	stime = time.time()
	myargs=my_argparse() 
	if not os.path.exists(myargs.outFolder):
		os.makedirs(myargs.outFolder)
	candidate_sites=pd.read_csv(myargs.sites_file,sep="\t",header=None)
	pool = multiprocessing.Pool(processes = myargs.processes)
	q = multiprocessing.Manager().Queue()
	overall_progressbar = tqdm.tqdm(total=candidate_sites.drop_duplicates([0,14]).shape[0], desc='extracting label...', position=0)
	def callback(result):
		overall_progressbar.update(1)
	for site,group in candidate_sites.groupby([0,14]):
		result=pool.apply_async(handle,(site,group,q,myargs.bam,myargs.label,), callback=callback)
	pool.close()
	pool.join()
	res=[]
	mapping_dict={1:"positive_dataset",0:"negative_dataset"}
	print("Summarizing results and writing to file...\n")
	for i in range(q.qsize()) :
		values=q.get()
		res.append(values)
	with open(os.path.join(myargs.outFolder,mapping_dict[myargs.label]),'wb') as file:
		joblib.dump(res,file)
	if myargs.plain_text_output :
		with open(os.path.join(myargs.outFolder,mapping_dict[myargs.label]+".txt"),'w') as file:
			for site in res:
				for a,b in zip(list(site.values())[0][0],list(site.values())[0][1]):
					file.write('\t'.join([list(site.keys())[0][0],str(list(site.keys())[0][1]),a,",".join([str(x) for x in b]),"\n"]))

	etime = time.time()
	delta=etime - stime
	print('\ntime taken by program\t',datetime.timedelta(seconds=delta))




'''
sites_file="/home/weir/output/vir_fip37_mtb/result_11_28/filter_final_result"
bamfile = pysam.AlignmentFile("/home/weir/m6a_model/result/align_result/VIRc/VIRc.sort.bam", "rb")
'''
