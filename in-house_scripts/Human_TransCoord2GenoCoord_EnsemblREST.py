# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 21:20:14 2021
@refence:
    @http://rest.ensembl.org
    @http://rest.ensembl.org/documentation/info/assembly_cdna
    @https://github.com/Ensembl/ensembl-rest/wiki
    @https://github.com/Ensembl/ensembl-rest/wiki/Getting-Started
@author: hqin
@discripution:
    @using the rest.ensembl tools/API and python to covert the transcriptome coordinate to genomic coordinate
"""
import requests, sys, json
import argparse, os
import time

def GenoLoci2TransLoci(inf, ouf):
	f2 = open(ouf, 'w')
	server = "http://rest.ensembl.org"
	for line in open(inf, 'r'):
		time.sleep(2)
		ID, Loci = line.rstrip().split()[0:2]
		ext = "/map/cdna/"+ str(ID)+"/"+ str(Loci)+ ".." + str(Loci)
#		print(ext)
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			if not r.ok:
				r.raise_for_status()
				sys.exit()
			decoded = r.json()
			dic=json.loads(str(decoded).replace("'","\""))
			genoLoci, strand, chr = dic["mappings"][0]["end"], dic["mappings"][0]["strand"], dic["mappings"][0]["seq_region_name"]
			strand = int(strand)
#			print(strand)
			if strand == 1:
				strand = "+"
			elif strand == -1:
				strand = "-"
			else:
				strand = "NA"
#			print(strand)
#			print(genoLoci)
			f2.write(line.rstrip()+"\t"+ str(chr) + "\t"+str(genoLoci)+"\t" + str(strand) + "\n")
			f2.flush()
		except:
			genoLoci, strand, chr  = "NA", "NA", "NA"
			f2.write(line.rstrip()+"\t"+ str(chr) + "\t"+str(genoLoci)+"\t" + str(strand) + "\n")
			f2.flush()
			continue
	f2.close()
	return f2

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='covert trans coordinate to geno coordinate')
	parser.add_argument('-i', '--input', required = True,help="trans coordinate files builded on reference cdna file from ensembl")
	parser.add_argument('-o', '--output', required = True, help="Output file")
#	parser.add_argument('--cpu', default=8,help='cpu number usage,default=8')
	args = parser.parse_args(sys.argv[1:])
	global FLAGS
	FLAGS = args
	GenoLoci2TransLoci(FLAGS.input, FLAGS.output)