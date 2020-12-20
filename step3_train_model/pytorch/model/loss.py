import torch.nn.functional as F
import torch.nn as nn


def CE_loss(output, target):
    return nn.CrossEntropyLoss()(output, target)
