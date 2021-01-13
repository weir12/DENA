import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
#from sklearn.metrics import precision_recall_curve
#from sklearn.metrics import plot_precision_recall_curve
FONT_SIZE=12
plt.rc('font', size=FONT_SIZE)          # controls default text sizes

plt.rc('axes', titlesize=FONT_SIZE)     # fontsize of the axes title

plt.rc('axes', labelsize=FONT_SIZE)    # fontsize of the x and y labels

plt.rc('xtick', labelsize=FONT_SIZE)    # fontsize of the tick labels

plt.rc('ytick', labelsize=FONT_SIZE)    # fontsize of the tick labels

plt.rc('legend', fontsize=FONT_SIZE)    # legend fontsize

plt.rc('figure', titlesize=FONT_SIZE)  # fontsize of the figure title
plt.rc('font',family='Arial') 
def ROC(outfile="ROC_curve.pdf",*files):
	'''
	files:multi files that store the predicted results and the true labels of validation set 
	'''
	plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim(0,1)
	plt.ylim(0,1)
	for file in files:
		label=file.split("/")[-2]
		data=pd.read_csv(file,sep="\t")
		fpr, tpr, thresholds = metrics.roc_curve(data['target'], data['output'], pos_label=1)
		auc=metrics.auc(fpr,tpr)
		plt.plot(fpr, tpr, alpha=0.8,label="{0} ( AUC={1:.3f} )".format(label, auc))
	plt.xlabel('False positive rate')
	plt.ylabel('True positive rate')
	plt.title('ROC curve')
	plt.legend(loc='best')
	plt.savefig(outfile,bbox_inches='tight',format='pdf')
	plt.close()    


def PRC(outfile="PR_curve.pdf",*files):
	'''
	files:multi files that store the predicted results and the true labels of validation set 
	'''
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.plot([0, 1], [1, 0], 'k--')
	for file in files:
		label=file.split("/")[-2]
		data=pd.read_csv(file,sep="\t")
		precision, recall, thresholds = metrics.precision_recall_curve(data['target'], data['output'])
		ap=metrics.average_precision_score(data['target'], data['output'])
		plt.plot(recall, precision, label="{0} ( AP={1:.3f} )".format(label, ap))
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.title('Precision-Recall curve')
	plt.legend(loc='best')
	plt.savefig(outfile,bbox_inches='tight',format='pdf')
	plt.close()
	
	