# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import argparse, sys, re, os
import matplotlib.pyplot as plt
import matplotlib_venn as mv
import seaborn as sns
from matplotlib.pyplot import MultipleLocator
import scipy.stats as sci
import venn

# from scipy.stats import spearmastat_funcnr
# import pylib
# import time
# from tqdm import tqdm, trange

def get_path(flist):
    """
    path/file1   name1   Chr,genoLoci,strand,geneid,...
    path/file2   name2   Chr,genoLoci,strand,geneid,...
    """
    files = list()
    for i, line in enumerate(open(flist,'r',encoding="utf-8")):
        file = line.strip().split()[0]
        files.append(file)
    path = os.path.dirname(files[0])
    return path

def data_group(df, axX, axY, sampleNames):
    """
    axY: ["ratio_x","ratio_y"]
    axX: "motif"
    From:
    motif    ratio_x    ratio_y
    AAACA    1    2
    To:
    motif    ratio    sample
    AAACA    1    mod
    AAACA    2    unmod
    """
    dfs = list()
    for index, ax in enumerate(axY):
        dfs.append("".join(["df", str(index)]))
        y = axY[index]
        dfs[index] = df[[axX, y]]
        dfs[index].rename(columns = {y:"ratio"}, inplace=True)
        dfs[index].loc[:, "sample"] = [sampleNames[index]]*len(df[axX])  #generate the list of sample name euqual to the length of df[axX]: ["mod",",mod",....]
    df = pd.DataFrame(columns=[axX, "ratio", "sample"]) #build null df
    for dfe in dfs:
        df = pd.concat([df, dfe], join = "inner", ignore_index = False, axis = 0) #merge all df in dfs
    print(df.columns.tolist())
    return df

def get_label(df, labels):
    """
    labels->label:
    Chr,genoLoci,strand -> Chr_genoLoci_strand
    """
    labs =list()
    lab = [str(i) for i in labels.strip().split(",")]
    for item in df.iterrows():# 按行遍历df
        item = item[1] #返回该行所有值，item[0]为返回改行index
        ele = "_".join([str(item[i]) for i in lab])
        labs.append(ele)
    return labs

def get_SpecialColnames(df,pattern):
    ColNames = df.columns.str.strip()
    Patterns = list()
    pat = "".join(["^", "(", pattern, ".*", ")"])
    for i in ColNames:
        ele = re.match(pat, i)
        if ele is not None:
            Patterns.append(ele.group())
    return Patterns

def decodeInputlist(flist):
    """
    file1   name1   Chr,genoLoci,strand,geneid,...
    file2   name2   Chr,genoLoci,strand,geneid,...
    """
    valueList, sampleNames, dfs, RatioValues = list(), list(), list(), list()
    for i, line in enumerate(open(flist,'r',encoding="utf-8")):
        file, SampleName, label = line.strip().split()[0:3]
        dfs.append ("".join(["df",str(i)]))
        lab = "_".join([str(i) for i in label.split(",")])
        dfs[i] = pd.read_csv(file,sep="\t", dtype = object)
        dfs[i].fillna("NA")
        RatioValues.append(dfs[i]["ratio"])
        if str(lab) not in dfs[i].columns.tolist():
            dfs[i]["label"] = get_label(dfs[i], label)
        value = dfs[i]["label"].drop_duplicates().values.tolist() #drop out duplicates of value
        valueList.append(value)
        sampleNames.append(SampleName)
    return dfs, sampleNames, valueList, RatioValues

def get_intersect(flist):
    dfs, sampleNames, valueList, RatioValues= decodeInputlist(flist)
    #find the intersect lines from four files
    if len(dfs) == 1:
        insdf = dfs[0]
    if len(dfs) == 2:
        insdf = pd.merge(dfs[0],dfs[1],on=["label"])
    elif len(dfs) == 3:
        insdf1 = pd.merge(dfs[0],dfs[1],on=["label"])
        insdf = pd.merge(insdf1,dfs[2],on=["label"])
    elif len(dfs) == 4:
        insdf1 = pd.merge(dfs[0],dfs[1],on=["label"])
        insdf2 = pd.merge(dfs[2],dfs[3],on=["label"])
        insdf = pd.merge(insdf1,insdf2,on=["label"])
    elif len(dfs) > 4:
        sys.exit("Please check filelist containing files must be less than 4.\n")
    Names = [str(i) for i in sampleNames]
    Names.append("intersectSites.txt")
    ouf = "_".join(Names)
    insdf.to_csv(ouf, sep='\t',index=False)
    return insdf, sampleNames

def Veen(flist):
    """
    file1   name1   Chr,genoLoci,strand,geneid,...
    file2   name2   Chr,genoLoci,strand,geneid,...
    """
    valueList, sampleNames, dfs = list(), list(), list()
    for i, line in enumerate(open(flist,'r',encoding="utf-8")):
        file, sampleName, label = line.strip().split()[0:3]
        dfs.append ("".join(["df",str(i)]))
        lab = "_".join([str(i) for i in label.split(",")])
        dfs[i] = pd.read_csv(file,sep="\t", dtype = object)
        dfs[i].fillna("NA")
        if str(lab) not in dfs[i].columns.tolist():
            dfs[i]["label"] = get_label(dfs[i], label)
        value = dfs[i]["label"].drop_duplicates().values.tolist() #drop out duplicates of value
        valueList.append(value)
        sampleNames.append(sampleName)
    ValueSet = list()
    sampleNametup = tuple(sampleNames)
    for i in range(len(valueList)):
        ValueSet.append(set(valueList[i])) #covert each element of valuelist from list() to set()
    my_dpi=150
    #控制图尺寸的同时，使图高分辨率（高清）显示
    plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)
    if len(valueList) > 4:
        sys.exit("Do not plot over 3 datasets")
    elif len(valueList) == 2:
        g = mv.venn2(subsets= ValueSet,
                    set_labels = sampleNametup,   #try ste_labels whether need a tuple() but ont list()
                    set_colors = ('#3d9a9b', '#c59a38'),
                    alpha = 0.8,
                    normalize_to = 1.0
                    )
    elif len(valueList) == 3:
        g = mv.venn3(subsets= ValueSet,
                    set_labels = sampleNametup,
                    set_colors = ('#3d9a9b', '#8d4c4e', '#c59a38'),
                    alpha = 0.8,
                    normalize_to = 1.0
                    )
    elif len(valueList) == 4:
        labels = venn.get_labels(list(ValueSet), fill=['number','percent'])
        g, ax = venn.venn4(labels, names=list(sampleNametup),fontsize=8)
    else:
        """
        pyvenn plot over three-dimensional venn figures
        https://zhuanlan.zhihu.com/p/195541937
        """
        pass
    if len(valueList) < 4:
        for text in g.set_labels:
            text.set_fontsize(8)
        for text in g.subset_labels:
            text.set_fontsize(6)
    Names = [str(i) for i in sampleNames]
    Names.append("VeenPlot.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight') #save as the .pdf figure

def Box_all(flist):
    dfs, sampleNames, valueList, RatioValues= decodeInputlist(flist)
    my_dpi=150 
    plt.figure(figsize=(150*len(valueList)/my_dpi, 450/my_dpi), dpi=my_dpi) #控制图尺寸的同时，使图高分辨率（高清）显示
    # sns.set(style="white", palette = "muted", color_codes=True) #背景白色，没有标线线,使用颜色就是传递参数给palette
    # sns.set_style("ticks") #xy轴都有非常短的小刻度
    current_palette = sns.color_palette("hls",len(sampleNames)) #就是这里穿参数，一般使用“hls”,l是亮度，s是饱和度,12表示12种颜色
    sns.boxplot(data = RatioValues, palette = current_palette, width = 0.2, fliersize= 0.6, linewidth= 0.8, showfliers=False)
    sns.despine(offset=0.3)
    ax=plt.gca();#获得坐标轴的句柄
    ax.spines['bottom'].set_linewidth(0.8);#设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(0.8);#设置左边坐标轴的粗细
    ax.set_ylim([-0.05, 1]) 
    ax.set_xlabel("Samples", fontsize=7) #设置坐标标签字体大小
    ax.set_ylabel("Ratio", fontsize=7)
    plt.xticks(range(0,len(sampleNames),1),labels = [sampleNames[0],sampleNames[1]],fontsize=6) #设置坐标刻度字体的大小
    plt.yticks(fontsize=6)
    Names = [str(i) for i in sampleNames]
    Names.append("BoxPlot_allSites.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight') #save as the .pdf figure

def Box_intersect(flist):
    df, sampleNames = get_intersect(flist)
    RatioValues = list()
    RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
    for val in RatioLabs:
        ratiovalue = df[val].values.tolist()
        RatioValues.append(ratiovalue)
    my_dpi=150 
    plt.figure(figsize=(150*len(RatioValues)/my_dpi, 450/my_dpi), dpi=my_dpi) 
    # sns.set(style="white", palette = "muted", color_codes=True)
    # sns.set_style("ticks")
    current_palette = sns.color_palette("hls",len(sampleNames))
    sns.boxplot(data = RatioValues, palette = current_palette, width = 0.2, fliersize= 0.6, linewidth= 0.8, showfliers=False)
    sns.despine(offset=0.3)
    ax=plt.gca();#获得坐标轴的句柄
    ax.spines['bottom'].set_linewidth(0.8);#设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(0.8);#设置左边坐标轴的粗细
    # ax.set_ylim([-0.05, 1]) 
    ax.set_xlabel("Samples", fontsize=7) #设置坐标标签字体大小
    ax.set_ylabel("Ratio", fontsize=7)
    plt.xticks(range(0,len(sampleNames),1),labels = sampleNames,fontsize=6) #设置坐标刻度字体的大小
    plt.yticks(fontsize=6)
    Names = [str(i) for i in sampleNames]
    Names.append("BoxPlot_IntersectSites.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight') #save as the .pdf figure

def classify_motif(flist):
    """
    file1   name1   Chr,genoLoci,strand,geneid,...
    file2   name2   Chr,genoLoci,strand,geneid,...
    """
    df, sampleNames = get_intersect(flist)
    RatioLabs = get_SpecialColnames(df,"ratio")
    motifColname = get_SpecialColnames(df,"motif")[0]
    motifs = df[motifColname].drop_duplicates().values.tolist()
    dfs = list()
    for index, motif in enumerate(motifs):
        dfs.append("".join(["df", str(index)]))
        dfs[index] = df[df[motifColname].str.match(motif)]
    return dfs, sampleNames, RatioLabs, motifs

def EachMotifRatio_Box(flist):
    dfs, sampleNames, RatioLabs, motifs = classify_motif(flist)
    for index, df in enumerate(dfs):
        motif = motifs[index]
        valueList2 = list()
        valueList = [df[str(RatioLabs[0])].tolist(), df[str(RatioLabs[1])].tolist()]
        for index, value in enumerate(valueList):
            valueList2.append(list(map(float, value)))
        my_dpi=150 
        plt.figure(figsize=(150*len(valueList)/my_dpi, 450/my_dpi), dpi=my_dpi) #控制图尺寸的同时，使图高分辨率（高清）显示
        current_palette = sns.color_palette("hls",len(valueList))
        sns.boxplot(data = valueList2, palette = current_palette, width = 0.2, fliersize= 0.6, linewidth= 0.8, showfliers=False)
        sns.despine(offset=0.3)
        ax=plt.gca();#获得坐标轴的句柄
        ax.spines['bottom'].set_linewidth(0.8);#设置底部坐标轴的粗细
        ax.spines['left'].set_linewidth(0.8);#设置左边坐标轴的粗细
        ax.set_ylim([-0.05, 1]) 
        ax.set_xlabel(motif, fontsize=7) #设置坐标标签字体大小
        ax.set_ylabel("Ratio", fontsize=7)
        plt.xticks(range(0,len(sampleNames),1),labels = [sampleNames[0],sampleNames[1]],fontsize=6) #设置坐标刻度字体的大小
        plt.yticks(fontsize=6)
        Names = [str(i) for i in sampleNames]
        Names.extend([motif, "BoxPlot.pdf"])
        ouf = "_".join(Names)
        plt.savefig(os.path.join(path, ouf), bbox_inches='tight') #save as the .pdf figure

def Density_all(flist): #plot all ratio of each sample although the numer of ratio in each sample
    """
    file1   name1   Chr,genoLoci,strand,geneid,...
    file2   name2   Chr,genoLoci,strand,geneid,...
    """
    dfs, SampleNames, valueList, RatioValues= decodeInputlist(flist)
    SampleNametup = tuple(SampleNames)
    current_palette = sns.color_palette("hls",len(SampleNames)) #就是这里穿参数，一般使用“hls”,l是亮度，s是饱和度,12表示12种颜色
    sns.set(style="white", palette = current_palette, color_codes=True) #背景白色，没有标线线,使用颜色就是传递参数给palette
    sns.set_style("ticks") #xy轴都有非常短的小刻度
    for index, val in enumerate(RatioValues):
        value = (list(map(float, list(val))))
        sns.distplot(value, label=SampleNametup[index], norm_hist=True, bins=20, axlabel = "Ratio")
    # plt.legend()
    sns.despine(offset=0.3)
    ax=plt.gca();#获得坐标轴的句柄
    ax.spines['bottom'].set_linewidth(0.8);#设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(0.8);#设置左边坐标轴的粗细 
    ax.set_xlabel(fontsize=7) #设置坐标标签字体大小
    ax.set_ylabel(fontsize=7)
    plt.xticks(fontsize=6) #设置坐标刻度字体的大小
    plt.yticks(fontsize=6)
    Names = [str(i) for i in SampleNames]
    Names.append("DensityPlot_allSites.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight') #save as the .pdf figure
    
def Density_intersect(flist): #plot the ratio of intersect sites among all samples although the numer of ratio in each sample
    df, sampleNames = get_intersect(flist)
    current_palette = sns.color_palette("hls",len(sampleNames))
    sns.set(style="white", palette = current_palette, color_codes=True)
    sns.set_style("ticks")
    try:
        RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
    except:
        sys.stderr.write("please select colnames of values in %s and change the codes.\n"%(df.columns.tolist()))
        sys.exit()
    for index, val in enumerate(RatioLabs):
        value = (list(map(float, list(df[val].values.tolist()))))
        sns.distplot(value, label=sampleNames[index], norm_hist=True, bins=20, axlabel = "Ratio")
    sns.despine(offset=0.3)
    ax=plt.gca();#获得坐标轴的句柄
    ax.spines['bottom'].set_linewidth(0.8);#设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(0.8);#设置左边坐标轴的粗细 
    ax.set_xlabel(fontsize=7) #设置坐标标签字体大小
    ax.set_ylabel(fontsize=7)
    plt.xticks(fontsize=6) #设置坐标刻度字体的大小
    plt.yticks(fontsize=6)
    Names = [str(i) for i in sampleNames]
    Names.append("DensityPlot_intersectSites.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight') #save as the .pdf figure

def get_corr(df,meth):
    RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
    # for ele in RatioLabs:
    #     df[ele] = df[ele].astype("float64")
    cor = df[RatioLabs].corr(method= meth)
    return cor

def Jointplot(flist):
    """
    file1   name1   Chr,genoLoci,strand,geneid,...
    file2   name2   Chr,genoLoci,strand,geneid,...
    """
    df, sampleNames = get_intersect(flist)
    current_palette = sns.color_palette("hls",len(sampleNames))
    try:
        RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
    except:
        sys.stderr.write("please select colnames of values in %s and change the codes.\n"%(df.columns.tolist()))
        sys.exit()
    for ele in RatioLabs:
        df[ele] = df[ele].astype("float64") #Notes: Values of ratio are not string but float64, so they must be covert to float64.
    cor = get_corr(df, "spearman") 
    cor = "%.3f"%(cor.iat[1,0]) #extract the value at second row and first column
    sns.set(style="white", palette = "muted", color_codes=True) #背景白色，没有标线线,使用颜色就是传递参数给palette
    sns.set_style("ticks") #xy轴都有非常短的小刻度
    # sns.set(font_scale = 1.5)#坐标轴刻度字体放大倍数
    p = sns.jointplot(data = df, x = RatioLabs[0], y = RatioLabs[1], kind="reg", marker=".", marginal_kws=dict(bins=25, fill=True),marginal_ticks=True, truncate=False, xlim=(0, 1), ylim=(0, 1), palette=current_palette, height=6)
    ####
    #add y=x line
    x1 = df[RatioLabs[0]].tolist()
    x2 = df[RatioLabs[1]].tolist()
    x0, x1 = p.ax_joint.get_xlim()
    y0, y1 = p.ax_joint.get_ylim()
    lims = [max(x0, y0), min(x1, y1)]
    p.ax_joint.plot(lims, lims, ':k', color = 'r')
    ####
    # p = sns.jointplot(data = df, x = RatioLabs[0], y = RatioLabs[1],hue="motif_x",height=5, xlim=(0, 1), ylim=(0, 1))
    p.set_axis_labels(sampleNames[0], sampleNames[1], fontsize=16)
    ax=plt.gca()
    ax.text(0,0.9, cor, fontsize=11,horizontalalignment='center')
    plt.xticks(fontsize=11) #设置右边柱状图纵坐标刻度字体的大小
    plt.yticks(fontsize=11)
    ####
    #add spearman and pvalue
    x1 = df[RatioLabs[0]].tolist()
    x2 = df[RatioLabs[1]].tolist()
    spearman, pvalue = sci.spearmanr(x1,x2) #pvalue=0.0 represent very small
    spearman, pvalue = ":".join(["spearman", str(spearman)]), ":".join(["pvalue", str(pvalue)])
    plt.title("; ".join([spearman, pvalue]))
    # plt.text(-10, 0.9,"; ".join([spearman, pvalue]))
    ####
    # x_major_locator=MultipleLocator(0.1)
    # y_major_locator=MultipleLocator(0.1)
    # ax.xaxis.set_major_locator(x_major_locator)
    # ax.yaxis.set_major_locator(y_major_locator)
    # plt.xlim(-0.05,1)
    # plt.ylim(0,1)
    Names = [str(i) for i in sampleNames]
    Names.append("JointPlot_intersectSites.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight')

def Box_motif(flist):
    """
    file1   name1   Chr,genoLoci,strand,geneid,...
    file2   name2   Chr,genoLoci,strand,geneid,...
    """
    df, sampleNames = get_intersect(flist)
    try:
        RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
        MotifLabs = get_SpecialColnames(df,"motif")
    except:
        sys.stderr.write("please select colnames of values in %s and change the codes.\n"%(df.columns.tolist()))
        sys.exit()
    # current_palette = sns.color_palette("hls",len(MotifLabs))
    for ele in RatioLabs:
        df[ele] = df[ele].astype("float64") #Notes: Values of ratio are not string but float64, so they must be covert to float64.
    for ele in MotifLabs:
        df[ele] = df[ele].astype("object") #Notes: covert to string.
    df = data_group(df, MotifLabs[0], RatioLabs, sampleNames) #re-construct the df
    sns.set(style="white", palette = "muted", color_codes=True) #背景白色，没有标线线,使用颜色就是传递参数给palette
    sns.set_style("ticks") #xy轴都有非常短的小刻度
    p = sns.boxplot(data = df, x = MotifLabs[0], y = "ratio", hue = "sample", palette="Set2")
    p.set_xticklabels(p.get_xticklabels(),rotation=30) 
    ax=plt.gca()
    plt.xticks(fontsize=11) #设置坐标刻度字体的大小
    plt.yticks(fontsize=11)
    ax.set_xlabel("Motifs",fontsize=15) #设置坐标标签字体大小
    ax.set_ylabel("Ratio",fontsize=15)
    Names = [str(i) for i in sampleNames]
    Names.append("boxPlot_intersectMotifRatio.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight')

def Catplot(flist):
    """
    file1   name1   Chr,genoLoci,strand,geneid,...
    file2   name2   Chr,genoLoci,strand,geneid,...
    """
    df, sampleNames = get_intersect(flist)
    # current_palette = sns.color_palette("hls",len(sampleNames))
    try:
        RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
        MotifLabs = get_SpecialColnames(df,"motif")
    except:
        sys.stderr.write("please select colnames of values in %s and change the codes.\n"%(df.columns.tolist()))
        sys.exit()
    for ele in RatioLabs:
        df[ele] = df[ele].astype("float64") #Notes: Values of ratio are not string but float64, so they must be covert to float64.
    for ele in MotifLabs:
        df[ele] = df[ele].astype("object") #Notes: covert to string.
    df = data_group(df, MotifLabs[0], RatioLabs, sampleNames) #re-construct the df
    sns.set(style="white", palette = "muted", color_codes=True) #背景白色，没有标线线,使用颜色就是传递参数给palette
    sns.set_style("ticks") #xy轴都有非常短的小刻度
    p = sns.catplot(data = df, x = MotifLabs[0], y = "ratio", hue = "sample", kind="swarm", palette="Set2") #multiple samples
    # p = sns.catplot(data = df, x = MotifLabs[0], y = "ratio", hue = "sample", kind="swarm", palette=['seagreen','peru']) #two sample for change colors
    # p = sns.catplot(data = df, x = MotifLabs[0], y = "ratio", hue = "sample", kind="swarm", palette=['peru']) #single sample
    # p.set_axis_labels(sampleNames[0], "Ratio", fontsize=16) ##设置坐标轴标签及字体大小
    p.set_xticklabels(rotation=30) 
    ax=plt.gca()
    plt.xticks(fontsize=11) #设置坐标刻度字体的大小
    plt.yticks(fontsize=11)
    ax.set_xlabel("Motifs",fontsize=15) #设置坐标标签字体大小
    ax.set_ylabel("Ratio",fontsize=15)
    plt.axhline(y=0.1,c='b',ls='--',lw=1)
    # plt.axvline(x=0.1,c='red',ls='--',lw=0.5)
    plt.ylim(0,1)
    Names = [str(i) for i in sampleNames]
    Names.append("CatPlot_intersectMotifRatio.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight')
    
def vioplot(flist):
    """
    plot violinplot
    """
    df, sampleNames = get_intersect(flist)
    try:
        RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
        MotifLabs = get_SpecialColnames(df,"motif")
    except:
        sys.stderr.write("please select colnames of values in %s and change the codes.\n"%(df.columns.tolist()))
        sys.exit()
    # current_palette = sns.color_palette("hls",len(MotifLabs))
    for ele in RatioLabs:
        df[ele] = df[ele].astype("float64") #Notes: Values of ratio are not string but float64, so they must be covert to float64.
    for ele in MotifLabs:
        df[ele] = df[ele].astype("object") #Notes: covert to string.
    df = data_group(df, MotifLabs[0], RatioLabs, sampleNames) #re-construct the df
    sns.set(style="white", palette = "muted", color_codes=True) #背景白色，没有标线线,使用颜色就是传递参数给palette
    sns.set_style("ticks") #xy轴都有非常短的小刻度
    p = sns.violinplot(data = df, x = "sample", y = "ratio", hue = "sample")
    p.set_xticklabels(p.get_xticklabels(),rotation=30) 
    ax=plt.gca()
    plt.xticks(fontsize=11) #设置坐标刻度字体的大小
    plt.yticks(fontsize=11)
    ax.set_xlabel("Motifs",fontsize=15) #设置坐标标签字体大小
    ax.set_ylabel("Ratio",fontsize=15)
    Names = [str(i) for i in sampleNames]
    Names.append("violinPlot_intersectRatio.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight')
    
def vioplot_Motif(flist):
    """
    plot violinplot
    """
    df, sampleNames = get_intersect(flist)
    try:
        RatioLabs = get_SpecialColnames(df,"ratio") #use pandas str.contains method to extract all colnames containing keyword "ratio"
        MotifLabs = get_SpecialColnames(df,"motif")
    except:
        sys.stderr.write("please select colnames of values in %s and change the codes.\n"%(df.columns.tolist()))
        sys.exit()
    # current_palette = sns.color_palette("hls",len(MotifLabs))
    for ele in RatioLabs:
        df[ele] = df[ele].astype("float64") #Notes: Values of ratio are not string but float64, so they must be covert to float64.
    for ele in MotifLabs:
        df[ele] = df[ele].astype("object") #Notes: covert to string.
    df = data_group(df, MotifLabs[0], RatioLabs, sampleNames) #re-construct the df
    sns.set(style="white", palette = "muted", color_codes=True) #背景白色，没有标线线,使用颜色就是传递参数给palette
    sns.set_style("ticks") #xy轴都有非常短的小刻度
    if len(sampleNames ) != 2:
        sys.stderr.write("violinplot_motif/beaplot only need two samples, suach as WT and KO.\n")
        sys.exit()
    p = sns.violinplot(data = df, x = df.columns[0], y = "ratio", hue = "sample", split=True)
    p.set_xticklabels(p.get_xticklabels(),rotation=30) 
    ax=plt.gca()
    plt.xticks(fontsize=11) #设置坐标刻度字体的大小
    plt.yticks(fontsize=11)
    ax.set_xlabel("Motifs",fontsize=15) #设置坐标标签字体大小
    ax.set_ylabel("Ratio",fontsize=15)
    Names = [str(i) for i in sampleNames]
    Names.append("violinPlot_intersectMotifRatio.pdf")
    ouf = "_".join(Names)
    plt.savefig(os.path.join(path, ouf), bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='python plots <"Veen", "Box_all","Box_intersect", "EachMotifRatio_Box", "Density_all", "Density_intersect", "Jointplot","Catplot", "Box_motif", "vioplot", "vioplot_Motif">,\
    flist:\
    fileName1    sampleName1    Chr,genoLoci,strand,geneid,...\
    fileName2    sampleName2    Chr,genoLoci,strand,geneid,...\
    ')
    parser.add_argument('-l', '--flist', required = True, help="Input the file list containing <fileName1    sampleName1    Chr,genoLoci,strand,geneid,...>")
    parser.add_argument('-c', '--choose', default="Box",help="select plot type, default plot Box")
    args = parser.parse_args(sys.argv[1:])
    
    flist = args.flist
    global path
    path = get_path(flist)
    
    typ = args.choose
    Typ = ("Veen", "Box_all", "Box_intersect", "EachMotifRatio_Box", "Density_all", "Density_intersect", "Jointplot", "Catplot", "Box_motif", "vioplot", "vioplot_Motif")
    if typ not in Typ:
        sys.exit("Do not plot %s, please re-choose figure type from %s"%(typ, str(Typ)))
    else:
        sys.stderr.write("You choose to plot %s.\n"%(typ))
        eval(typ)(flist)
        # plt.show() #show the figure
        plt.close()
    sys.stderr.write("%s figure has save as pdf.\n"%(typ))