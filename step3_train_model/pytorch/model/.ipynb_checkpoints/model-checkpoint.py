import torch.nn as nn
import torch
import torch.nn.functional as F
from base import BaseModel
from .layers import Inception_ResNetv2

class DENAModel(BaseModel):
	def __init__(self, in_dim, hidden_dim, n_layer, n_classes,drop_out=False):
		super(DENAModel, self).__init__()
		self.lstm = nn.LSTM(in_dim, hidden_dim, n_layer, batch_first=True,bidirectional=True)
		self.cnn=Inception_ResNetv2(in_channels=1, classes=2)
		self.classifier=nn.Sequential(
			nn.Linear(hidden_dim*2+769,1024),
			nn.ReLU6(inplace=True),
			nn.Linear(1024,2))

	def forward(self, x_event,x_signal):
		out, (h_n, c_n) = self.lstm(x_event)
		# 此时可以从out中获得最终输出的状态h
		x_event = out[:, -1, :]
		x_signal=self.cnn(x_signal)
		x = self.classifier(torch.cat((x_signal,x_event),1))
		return x
