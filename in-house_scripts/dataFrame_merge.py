# -*- coding: utf-8 -*-
import pandas as pd
import argparse, sys

def decodeInputList(flist):
    labels, dfs, files, sampleNames = list(), list(), list(), list()
    for i, line in enumerate(open(flist,'r',encoding="utf-8")):
        try:
            file, sampleName, label = line.strip().split()
            dfs.append ("".join(["df",str(i)]))
            files.append(file)
            sampleNames.append(sampleName)
            dfs[i] = pd.read_csv(file,sep="\t", dtype = object)
            dfs[i].fillna("NA")
            labels.append(label.split(","))
        except:
            sys.stderr.write("Please check the input fileList whether is as: fileName  sampleName  feture1,feture2,feture3,...\n")
    return dfs, labels, files, sampleNames

def get_intersect(dfs, labels, select):
    # dfs, labels, files, sampleNames = decodeInputList(flist)
    #find the intersect lines from four files
    if len(dfs) <= 1:
        sys.exit("Please check filelist contains two files at least.\n")
    if len(dfs) == 2:
        ins = pd.merge(dfs[0],dfs[1],on=labels[0], how=str(select))
    elif len(dfs) == 3:
        ins1 = pd.merge(dfs[0],dfs[1],on=labels[0], how=str(select))
        ins = pd.merge(ins1,dfs[2],on=labels[0], how=str(select))
    elif len(dfs) == 4:
        ins1 = pd.merge(dfs[0],dfs[1],on=labels[0], how=str(select))
        ins2 = pd.merge(dfs[2],dfs[3],on=labels[0], how=str(select))
        ins = pd.merge(ins1,ins2,on=labels[0], how=str(select))
    elif len(dfs) > 4:
        sys.exit("Please check filelist containing files must be less than 4.\n")
    label = []
    for index, row in ins.iterrows():
        label.append( "_".join([str(row[i]) for i in labels[0]]))
    ins["label"] = label
    return ins
  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='covert trans coordinate to geno coordinate, inputlist as eg: fileName    sampleName    Chr,genoLoci,strand,geneid')
    parser.add_argument('-i', '--filelist', required = True, help="Input the files list with path and the columns used to merge")
    parser.add_argument('-s', '--select', default="inner", help="choose the inner (intersect) or outer (union) ")
    args = parser.parse_args(sys.argv[1:])
    ###############
    dfs, labels, files, sampleNames = decodeInputList(args.filelist)
    ins= get_intersect(dfs, labels, args.select)
    #use for DENA results
    slt = ""
    if args.select == "outer":
        slt = "union"
    elif args.select == "inner":
        slt = "intersect"
    sampleNames.append("%s_DFmerge.txt"%(slt))
    ins.to_csv("_".join(sampleNames), sep='\t',index=False) #write reults to csv output file
    #use for Nanom6A results
    # outf = ins[["chr", "loci","ids","modN","totN","ratio", "motif", "label"]]
    # outf.columns=["Chr", "genoloci","geneid","modRN","totalRN","ratio", "motif", "label"]
    # outf.to_csv("_".join([files[0].split(".txt")[0], "motif.txt"]), sep='\t',index=False) #write reults to csv output file
    ###############