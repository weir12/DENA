# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 14:28:48 2021

@author: hqin

Description:
    This script main used to extract the intersect sites with one liminted depth(eg: totalRN>=10) in all reps, 
    in order to used varify the accuracy of our LSTM models among reps with the same sites with same depth.
    This can avoid the effect of different depth on predicted m6A sites among reps.
"""

import pandas as pd
import argparse, sys, os
import dataFrame_merge as dFm
import Python_matlibPlot_V2 as pmp
import matplotlib.pyplot as plt
import matplotlib_venn as mv
import venn
# import seaborn as sns

def get_path(flist):
    """
    path/file1   name1   Chr,genoLoci,strand,geneid,...
    path/file2   name2   Chr,genoLoci,strand,geneid,...
    """
    files = list()
    for i, line in enumerate(open(flist,'r',encoding="utf-8")):
        file = line.strip().split("\t")[0]
        files.append(file)
    path = os.path.dirname(files[0])
    return path

def get_label(df, labels):
    """
    labels->label:
    Chr,genoLoci,strand -> Chr_genoLoci_strand
    """
    labs =list()
    for item in df.iterrows():# 按行遍历df
        item = item[1] #返回该行所有值，item[0]为返回改行index
        ele = "_".join([str(item[i]) for i in labels])
        labs.append(ele)
    return labs

def Extract_DepthSites(dfsList, labelList, sampleNameList, depth, ratio):
    dfs = list()
    label = labelList[0]
    for i, df in enumerate(dfsList):
        df["label"] = get_label(df, label)
        #df["totalRN"] = df["totalRN"].astype("int64")
        df["totalRN"] = df["totalRN"].astype("float64").astype("int64")
        df["ratio"] = df["ratio"].astype("float64")
        # df = df.loc[df["totalRN"] >= df["totalRN"].mean() & df["ratio" >= 0.1]]
        df = df[(df["totalRN"] >= float(depth)) & (df["ratio"] >= float(ratio))] 
        sampleName = sampleNameList[i]
        tmp = "LSTM_genoLoci_T%sR%s.txt"%(str(depth), str(ratio))
        ouf = "_".join([sampleName, tmp])
        oufPath = "/".join([path, ouf])
        print(oufPath)
        df.to_csv(oufPath, sep = "\t", index = False)
        dfs.append(df)
    return dfs

def extract_label(dfsList):
    labelValList = list()
    for i, df in enumerate(dfsList):
        labelVal = df["label"].drop_duplicates().values.tolist()
        labelValList.append(labelVal)
    return labelValList

def Veen(dfsList, sampleNames):
    ValueSet = list()
    labelValList = extract_label(dfsList)
    sampleNametup = tuple(sampleNames)
    for i in range(len(labelValList)):
        ValueSet.append(set(labelValList[i])) #covert each element of valuelist from list() to set()
    my_dpi=150
   
    plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)
    if len(labelValList) > 4:
        sys.exit("Do not plot over 3 datasets")
    elif len(labelValList) == 2:
        g = mv.venn2(subsets= ValueSet,
                    set_labels = sampleNametup,   #try ste_labels whether need a tuple() but ont list()
                    set_colors = ('#3d9a9b', '#c59a38'),
                    alpha = 0.8,
                    normalize_to = 1.0
                    )
    elif len(labelValList) == 3:
        g = mv.venn3(subsets= ValueSet,
                    set_labels = sampleNametup,
                    set_colors = ('#3d9a9b', '#8d4c4e', '#c59a38'),
                    alpha = 0.8,
                    normalize_to = 1.0
                    )
    elif len(labelValList) == 4:
        labels = venn.get_labels(list(ValueSet), fill=['number','percent'])
        g, ax = venn.venn4(labels, names=list(sampleNametup),fontsize=8)
    else:
        """
        pyvenn plot over four-dimensional venn figures
        https://zhuanlan.zhihu.com/p/195541937
        """
        pass
    # for text in g.set_labels:
        # text.set_fontsize(8)
    # for text in g.subset_labels:
        # text.set_fontsize(6)

def Df_merge(dfs, labelList, sampleNames, depth, ratio):
    label = labelList[0]
    if len(dfs) <= 1:
        sys.exit("Please check filelist contains two files at least.\n")
    if len(dfs) == 2:
        ins = pd.merge(dfs[0],dfs[1],on=label)
    elif len(dfs) == 3:
        ins1 = pd.merge(dfs[0],dfs[1],on=label)
        ins = pd.merge(ins1,dfs[2],on=label)
    elif len(dfs) == 4:
        ins1 = pd.merge(dfs[0],dfs[1],on=label)
        ins2 = pd.merge(dfs[2],dfs[3],on=label)
        ins = pd.merge(ins1,ins2,on=label)
    elif len(dfs) > 4:
        sys.exit("Please check filelist containing files must be less than 4.\n")
    lab = []
    for index, row in ins.iterrows():
        lab.append( "_".join([str(row[i]) for i in label]))
    ins["label"] = lab
    tmp = "LSTM_genoLoci_T%sR%s_intersectSites.txt"%(str(depth), str(ratio))
    Names = [str(i) for i in sampleNames]
    Names.append(tmp)
    ouf = "_".join(Names)
    #oufPath = "/".join([path, ouf])
    print(os.path.join(path, ouf))
    ins.to_csv(os.path.join(path, ouf), sep="\t", index = False)
    return ins

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='deal LSTM212 ori output files between reps')
    parser.add_argument('-l', '--filelist', required = True, help="input the filelist containing: fileName  sampleName  feture1,feture2,feture3,...")
    parser.add_argument('-d', '--depth', default = "10,20,30,40,50", help="Set the list of depth, eg.: 10,20,30,40,50...")
    parser.add_argument('-r', '--ratio', default = "0", help="Set the list of depth, eg.: 10,20,30,40,50...")
    #parser.add_argument('-p', '--path', help="set the path of output files")
    args = parser.parse_args(sys.argv[1:])
    
    global path 
    flist = args.filelist
    path = get_path(flist)
    # path = eval(repr(path).replace('\\', '/'))
    # path = path.replace('\\', '/')
    depth = args.depth
    ratio = args.ratio
    depthes = depth.strip().split(",")
    dfsList, labelList, fileList, sampleNameList = dFm.decodeInputList(flist) #for Df_merge
    for depth in depthes:
        dfs = Extract_DepthSites(dfsList, labelList, sampleNameList, depth, ratio)
        Df_merge(dfs, labelList, sampleNameList, depth, ratio)      
#for veen plot     
        Veen(dfs, sampleNameList)
        tmp = "_".join(["T%sR%s"%(depth, ratio),"VeenPlot.pdf"])
        Names = [str(i) for i in sampleNameList]
        Names.append(tmp)
        ouf = "_".join(Names)
        #oufPath = "/".join([path, ouf])
        print(os.path.join(path, ouf))
        plt.savefig(os.path.join(path, ouf), bbox_inches='tight') #save as the .pdf figure
