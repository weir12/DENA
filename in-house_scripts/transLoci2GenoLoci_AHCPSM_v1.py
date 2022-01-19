# -*- coding: utf-8 -*-
"""
Created on Tue Mar 9 21:20:14 2021
@convert transLoci to genoLoci base on refrence of gtf annotation
@author: Hang Qin
Platform: Linux, python3.6
Notes: the gtf file must be consistent with the transcriptome reference used for mapping in running DENA, such as download the cdna.fa and .gtf files of the same release version from Ensembl database
"""
import argparse, os, sys
import pandas as pd

def CalDistance(dic_exonLen, line):
    transID, transLoci, motif, modRN, totalRN, ratio = line[:]
    #<dic_exonLen> {'AT1G01010.1': [283, 281, 610]}
    exonLen = dic_exonLen[transID]
    length, exonID, distance = 0, 0, 0
    res = dict()
    for ele in range(len(exonLen)):
        length += int(exonLen[ele])
        if length <= int(transLoci):
            continue
        else:
            distance = length - int(transLoci)
            exonID = ele
            res={"transID":transID, "distance":distance, "exonID":exonID}
            break
    return res

def get_GenoLoci(df, inf, ouf):
    df, inf, ouf = df, inf, ouf
    """
    <df>:
    Chr    start    end    strand    geneid    transid    exonid    geneName    exonLength
    1    3631    3913    +    AT1G01010    AT1G01010.1    1    NAC001    283
    1    3996    4276    +    AT1G01010    AT1G01010.1    2    NAC001    281
    1    4486    5095    +    AT1G01010    AT1G01010.1    3    NAC001    610
    """
    dic_exonLen = df.groupby('transid').exonLength.apply(list).to_dict() #<dic_exonLen> {'AT1G01010.1': [283, 281, 610]}
    dic_leftBoundary = df.groupby('transid').start.apply(list).to_dict() #<dic_exonLen> {'AT1G01010.1': [3631, 3996, 4486]}
    dic_rightBoundary = df.groupby('transid').end.apply(list).to_dict() #<dic_exonLen> {'AT1G01010.1': [3913, 4276, 5095]}
    all_transids = set(df['transid'].tolist()) #get all transids in the GTF
    dic_strand = df.groupby('transid').strand.apply(set).apply(list).to_dict() #<dic_strand> {'AT1G01010.1': ["+"]}
    #dic_strand = df[['transid','strand']].set_index('transid').to_dict()['strand'] #<dic_strand> {'AT1G01010.1': [+, +, +]}
    dic_geneid = df.groupby('transid').geneid.apply(set).apply(list).to_dict() #<dic_geneid> {'AT1G01010.1': ["AT1G01010"]}
    dic_Chr = df.groupby('transid').Chr.apply(set).apply(list).to_dict() #<dic_Chr> {'AT1G01010.1': ["1"]}
    dic_geneName = df.groupby('transid').geneName.apply(set).apply(list).to_dict() #<dic_geneName> {'AT1G01010.1': ["NAC001"]}

    df_inf = pd.read_csv(inf, sep="\t", header = None, low_memory=False, error_bad_lines=False) #read the input files, and "error_bad_lines=False" can remove these lines more than 7 columns
    df_inf.columns=['transid','transLoci','motif','modRN','totalRN','ratio'] #add colume names
    """
    transid    transLoci    modRN    totalRN    ratio
    AT1G01620.1    606    AAACA    13    376    0.034574468085106384
    AT1G01820.1    1203    AAACA    10    26    0.38461538461538464
    AT1G02130.1    514    AAACA    11    56    0.19642857142857142
    """
    Strand, GenoLoci, Chrom, GeneID, GeneName = [], [], [], [], []
    GenoLoci = []
    for index, row in df_inf.iterrows():
        genoLoci = 0
        transID, transLoci, motif, modRN, totalRN, ratio = row[:]
        try:
            if str(transID) in all_transids:
                Chr = str(dic_Chr[transID][0])
                geneid = str(dic_geneid[transID][0])
                geneName = str(dic_geneName[transID][0])
                strand = str(dic_strand[transID][0])
                line = list(row)
                dicLoci = CalDistance(dic_exonLen, line) #{"transID":transID, "distance":distance, "exonID":exonID}
                exonNum = int(dicLoci["exonID"])
                if strand == "+":
                    genoLoci = int(dic_rightBoundary[transID][exonNum]) - int(dicLoci["distance"])
                elif strand == "-":
                    genoLoci = int(dic_leftBoundary[transID][exonNum]) + int(dicLoci["distance"])
                else:
                    strand, genoLoci = str(strand), str("NA")
            else:
                strand, genoLoci, Chr, geneid, geneName = str("NA"), str("NA"), str("NA"), str("NA"), str("NA")
        except:
            strand, genoLoci, Chr, geneid, geneName = str("Outlier"), str("Outlier"), str("Outlier"), str("Outlier"), str("Outlier")
            GenoLoci.append(genoLoci)
            Strand.append(strand)
            Chrom.append(Chr)
            GeneID.append(geneid)
            GeneName.append(geneName)
            continue
        GenoLoci.append(genoLoci)
        Strand.append(strand)
        Chrom.append(Chr)
        GeneID.append(geneid)
        GeneName.append(geneName)
    df_inf.insert(loc = 0, column = "Chr", value = Chrom)
    df_inf.insert(loc = 1, column = "genoLoci", value = GenoLoci)
    df_inf.insert(loc = 2, column = "strand", value= Strand)
    df_inf.insert(loc = 3, column = "geneid", value = GeneID)
    df_inf["geneName"] = GeneName
    df_inf.to_csv(ouf, sep='\t',index=False)
    """
    Chr    genoLoci    strand    geneid    transid    transLoci    modRN    totalRN    ratio    geneName
    """

def Tair10(GTF): #Arabidopsis_thaliana.TAIR10.50.gtf
    Gtf = GTF
    """
    deal with the gtf files with shell
    """
    cmd = str('''awk '{if($3=="exon"){print $0}}' %s > gtf.deal'''%(Gtf))
    os.system(cmd)
    os.system(r'''sed 's/gene_id //g' gtf.deal |sed 's/transcript_id //g'|sed 's/exon_number //g'|sed 's/gene_name //g'|sed 's/"//g'|sed 's/;//g'|sed 's/ /\t/g' > gtf.deal1''')
    os.system(r"cut -f1,4,5,7,9,10,11,12 gtf.deal1 > gtf.deal")
    """
    #gtf.deal:
    Chr    start    end    strand    geneid    transid    exonid    geneName
    1    3631    3913    +    AT1G01010    AT1G01010.1    1    NAC001
    """
    df = pd.read_csv("gtf.deal", sep="\t", header = None, low_memory=False)
    df[8] = df[2] - df[1] + 1
    df.columns=['Chr','start','end','strand','geneid','transid','exonid','geneName','exonLength']
    """
    <df>
    Chr    start    end    strand    geneid    transid    exonid    geneName    exonLength
    1    3631    3913    +    AT1G01010    AT1G01010.1    1    NAC001    283
    """
    tmpf = "_".join([Gtf, "deal"])
    df.to_csv(tmpf, sep='\t',index=False)
    os.system(r"rm gtf.deal gtf.deal1")

def S_cerevisiae(GTF): #Saccharomyces_cerevisiae.R64-1-1.50.gtf
    Tair10(GTF)

def S_pombe(GTF): # Schizosaccharomyces_pombe.ASM294v2.50.gtf
    Tair10(GTF)

def Celegans101(GTF): #Caenorhabditis_elegans.WBcel235.103.gtf
    Tair10(GTF)

def Human38(GTF): #Homo_sapiens.GRCh38.103.gtf
    Gtf = GTF #Homo_sapiens.GRCh38.103.gtf
    cmd = str('''awk '{if($3=="exon"){print $0}}' %s > gtf.deal'''%(Gtf))
    os.system(cmd)
    os.system(r'''sed 's/gene_id //g' gtf.deal |sed 's/gene_version //g'|sed 's/transcript_id //g' |sed 's/"; transcript_version "/./g' |sed 's/exon_number //g'|sed 's/gene_name //g'|sed 's/"//g'|sed 's/;//g'|sed 's/ /\t/g' > gtf.deal1''')
    os.system(r"cut -f1,4,5,7,9,11,12,13 gtf.deal1 > gtf.deal")
    df = pd.read_csv("gtf.deal", sep="\t", header = None, low_memory=False)
    df[8] = df[2] - df[1] + 1
    df.columns=['Chr','start','end','strand','geneid','transid','exonid','geneName','exonLength']
    tmpf = "_".join([Gtf, "deal"])
    df.to_csv(tmpf, sep='\t',index=False)
    os.system(r"rm gtf.deal gtf.deal1")

def Mouse(GTF): #Mus_musculus.GRCm39.103.gtf
    Human38(GTF)
    
def Zebrefish(GTF): #Danio_rerio.GRCz11.103_Zebrafish.gtf
    Human38(GTF)

def Pop_tri(GTF): #Populus_trichocarpa.Pop_tri_v3.50.gtf  
    Gtf = GTF 
    cmd = str('''awk '{if($3=="exon"){print $0}}' %s > gtf.deal'''%(Gtf))
    os.system(cmd)
    os.system(r'''sed 's/gene_id //g' gtf.deal |sed 's/transcript_id //g'|sed 's/exon_number //g'|sed 's/gene_source //g'|sed 's/"//g'|sed 's/;//g'|sed 's/ /\t/g' > gtf.deal1''')
    os.system(r"cut -f1,4,5,7,9,10,11,12 gtf.deal1 > gtf.deal")
    df = pd.read_csv("gtf.deal", sep="\t", header = None, low_memory=False)
    df[8] = df[2] - df[1] + 1
    df.columns=['Chr','start','end','strand','geneid','transid','exonid','geneName','exonLength']
    tmpf = "_".join([Gtf, "deal"])
    df.to_csv(tmpf, sep='\t',index=False)
    os.system(r"rm gtf.deal gtf.deal1")

def TransLoci2GenoLoci():
    species = ["Tair10", "Human38", "Celegans101", "Pop_tri", "S_cerevisiae", "S_pombe", "Mouse", "Zebrefish"]
    slt = FLAGS.select
    gtf, inf, ouf = FLAGS.gtf, FLAGS.input, FLAGS.output
    if slt not in species:
        sys.exit("please select the right species(Tair10, Human38, Celegans101, Pop_tri, S_cerevisiae, S_pombe), default: Tair10")
    else:
        sys.stderr.write("You select the %s species!\n"%(slt))
        tmpf = str("_".join([gtf, "deal"]))
        if os.path.isfile(tmpf):
            sys.stderr.write("%s has existed.\n"%(tmpf))
        else:
            eval(slt)(gtf)
        df = pd.read_csv(tmpf, sep="\t", low_memory=False)
        get_GenoLoci(df, inf, ouf)
        sys.stderr.write("%s has finished the conversion of transLoci to genoLoci.\n"%(inf))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='covert trans coordinate to geno coordinate. Notes: the gtf file must be consistent with the transcriptome reference used for mapping in running DENA, such as download the cdna.fa and .gtf files of the same release version from Ensembl database')
    parser.add_argument('-s', '--select', default="Tair10", help="select the species <Tair10, Human38, Celegans101, Pop_tri, S_cerevisiae, S_pombe, Mouse, Zebrefish>, default: Tair10")
    parser.add_argument('-g', '--gtf', required = True,help="the gtf file download from ensembl database")
    parser.add_argument('-i', '--input', required = True,help="trans coordinate files from DENA predict, containing six columns <transid,transLoci,motif,modRN,totalRN,ratio>")
    parser.add_argument('-o', '--output', required = True, help="name of Output file")
    args = parser.parse_args(sys.argv[1:])
    global FLAGS
    FLAGS = args
    TransLoci2GenoLoci()
