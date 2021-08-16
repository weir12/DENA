import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# 设置风格，seaborn有5种基本风格，context表示环境
sns.set(style="white", context="notebook")
# 处理中文问题
sns.set_style('whitegrid', {'font.sans-serif':['Arial']})
data = pd.read_csv("./motif_dis.txt",sep=" ",name=['freq','motif'])
data=data.sort_values(by='freq',ascending=[False])
fig = plt.figure(figsize=(10,5))
sns.barplot(x="motif",y="freq",data=data, palette="BuPu_r")
plt.savefig("motif_dist.pdf")
plt.close()
plt.close()
fig = plt.figure(figsize=(10,10))
plt.pie(data['freq'], labels=data['motif'],palette="BuPu_r")
plt.savefig("motif_dist_pie.pdf")
vali_res=pd.read_csv("/home/haopei/DENA/step3_train_model/ml/valid_performance.txt",sep="\t"
fig = plt.figure(figsize=(10,5))
sns.barplot(x="motif",y="values",hue="metrics",data=data, palette="BuPu_r")
sns.barplot(x="motif",y="values",hue="metrics",data=vali_res, palette="BuPu_r")
plt.savefig("XGboost_vali.pdf")