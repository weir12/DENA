import torch
from sklearn.metrics import auc

def accuracy(output, target):
    with torch.no_grad():
        pred = torch.argmax(output, dim=1)
        assert pred.shape[0] == len(target)
        correct = 0
        correct += torch.sum(pred == target).item()
    return correct / len(target)


def top_k_acc(output, target, k=3):
    with torch.no_grad():
        pred = torch.topk(output, k, dim=1)[1]
        assert pred.shape[0] == len(target)
        correct = 0
        for i in range(k):
            correct += torch.sum(pred[:, i] == target).item()
    return correct / len(target)

def auc(output, target):
	output,target=output.numpy(),target.numpy()
	fpr, tpr, thresholds = metrics.roc_curve(target, output, pos_label=1)
	return metrics.auc(fpr, tpr)
