import torch
import torch.nn as nn

class INCEPTION(nn.Module):
	def __init__(self,in_channels,outputs_channels,scale=1):
		super(INCEPTION,self).__init__()
		self.scale=scale
		self.branch1=nn.Sequential(
			Conv1d(in_channels=in_channels, out_channels=outputs_channels, padding=0,kernel_size=1,stride=1)
			)
		self.branch2=nn.Sequential(
		Conv1d(in_channels=in_channels, out_channels=outputs_channels, padding=0,kernel_size=1,stride=1),
		Conv1d(in_channels=outputs_channels, out_channels=outputs_channels,kernel_size=3,padding=1,stride=1)
		)
		self.branch3=nn.Sequential(
			Conv1d(in_channels=in_channels, out_channels=outputs_channels, padding=0,kernel_size=1,stride=1),
			Conv1d(in_channels=outputs_channels, out_channels=outputs_channels, kernel_size=5,padding=2,stride=1)
			)
		self.branch4=nn.Sequential(
			nn.MaxPool1d(kernel_size=3,padding=1,stride=1),
			Conv1d(in_channels=in_channels, out_channels=outputs_channels, padding=0,kernel_size=1,stride=1)
			)
		self.shortcut=Conv1d(in_channels=in_channels, out_channels=outputs_channels*4,kernel_size=1,stride=1,padding=0)
	def forward(self, x):
		out1 = self.branch1(x)
		out2 = self.branch2(x)
		out3 = self.branch3(x)
		out4 = self.branch4(x)
		out = torch.cat([out1, out2, out3, out4], dim=1)
		out=out*self.scale+self.shortcut(x)
		out=F.relu6(out)
		return out
class Conv1d(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, padding, stride=1, bias=False):
        super(Conv1d, self).__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size, stride=stride, padding=padding, bias=bias)
        self.bn = nn.BatchNorm1d(out_channels, eps=0.0001, momentum=0.1)
        self.relu = nn.ReLU(inplace=True)
    def forward(self, x):
        x = self.conv(x)
        x = self.bn(x)
        x = self.relu(x)
        return x
class Reduction_A(nn.Module):
    # 35 -> 17
    def __init__(self, in_channels, k, l, m, n):
        super(Reduction_A, self).__init__()
        self.branch_0 = Conv1d(in_channels, n, 3, stride=2, padding=0, bias=True)
        self.branch_1 = nn.Sequential(
            Conv1d(in_channels, k, 1, stride=1, padding=0, bias=True),
            Conv1d(k, l, 3, stride=1, padding=1, bias=True),
            Conv1d(l, m, 3, stride=2, padding=0, bias=True),
        )
        self.branch_2 = nn.MaxPool1d(3, stride=2, padding=0)

    def forward(self, x):
        x0 = self.branch_0(x)
        x1 = self.branch_1(x)
        x2 = self.branch_2(x)
        return torch.cat((x0, x1, x2), dim=1) # [b,1088,17]

class Stem(nn.Module):
    def __init__(self, input_size):
        super(Stem, self).__init__()
        self.features = nn.Sequential(
            Conv1d(1, 32, 3, stride=2, padding=0, bias=True), #  [b, 32, 149]
            Conv1d(32, 32, 3, stride=1, padding=0, bias=True), # 147 x 32 
            Conv1d(32, 64, 3, stride=1, padding=1, bias=True), # 147  x 64
            nn.MaxPool1d(3, stride=2, padding=0), # 73 x 64
            Conv1d(64, 80, 1, stride=1, padding=0, bias=True), # 73 x  80
            Conv1d(80, 192, 3, stride=1, padding=0, bias=True), # 71 x 192
            nn.MaxPool1d(3, stride=2, padding=0), # 35 x 192
        )
        self.branch_0 = Conv1d(192, 96, 1, stride=1, padding=0, bias=True)
        self.branch_1 = nn.Sequential(
            Conv1d(192, 48, 1, stride=1, padding=0, bias=True),
            Conv1d(48, 64, 5, stride=1, padding=2, bias=True),
        )
        self.branch_2 = nn.Sequential(
            Conv1d(192, 64, 1, stride=1, padding=0, bias=True),
            Conv1d(64, 96, 3, stride=1, padding=1, bias=True),
            Conv1d(96, 96, 3, stride=1, padding=1, bias=True),
        )
        self.branch_3 = nn.Sequential(
            nn.MaxPool1d(3, stride=1, padding=1),
            Conv1d(192, 64, 1, stride=1, padding=0, bias=True)
        )
    def forward(self, x):
        x = self.features(x)
        x0 = self.branch_0(x)
        x1 = self.branch_1(x)
        x2 = self.branch_2(x)
        x3 = self.branch_3(x)
        return torch.cat((x0, x1, x2, x3), dim=1)


class Inception_ResNet_A(nn.Module):
    def __init__(self, in_channels, scale=1.0):
        super(Inception_ResNet_A, self).__init__()
        self.scale = scale
        self.branch_0 = Conv1d(in_channels, 32, 1, stride=1, padding=0, bias=True)
        self.branch_1 = nn.Sequential(
            Conv1d(in_channels, 32, 1, stride=1, padding=0, bias=True),
            Conv1d(32, 32, 3, stride=1, padding=1, bias=True)
        )
        self.branch_2 = nn.Sequential(
            Conv1d(in_channels, 32, 1, stride=1, padding=0, bias=True),
            Conv1d(32, 48, 3, stride=1, padding=1, bias=True),
            Conv1d(48, 64, 3, stride=1, padding=1, bias=True)
        )
        self.conv = nn.Conv1d(128, 320, 1, stride=1, padding=0, bias=True)
        self.relu = nn.ReLU(inplace=True)
    def forward(self, x):
        x0 = self.branch_0(x)
        x1 = self.branch_1(x)
        x2 = self.branch_2(x)
        x_res = torch.cat((x0, x1, x2), dim=1)
        x_res = self.conv(x_res)
        return self.relu(x + self.scale * x_res)


class Inception_ResNet_B(nn.Module):
    def __init__(self, in_channels, scale=1.0):
        super(Inception_ResNet_B, self).__init__()
        self.scale = scale
        self.branch_0 = Conv1d(in_channels, 192, 1, stride=1, padding=0, bias=True)
        self.branch_1 = nn.Sequential(
            Conv1d(in_channels, 128, 1, stride=1, padding=0, bias=True),
            Conv1d(128, 160, 7, stride=1, padding=3, bias=True),
            Conv1d(160, 192, 1, stride=1, padding=0, bias=True)
        )
        self.conv = nn.Conv1d(384, 1088, 1, stride=1, padding=0, bias=True)
        self.relu = nn.ReLU(inplace=True)
    def forward(self, x):
        x0 = self.branch_0(x)
        x1 = self.branch_1(x)
        x_res = torch.cat((x0, x1), dim=1)
        x_res = self.conv(x_res)
        return self.relu(x + self.scale * x_res)


class Reduciton_B(nn.Module):
    def __init__(self, in_channels):
        super(Reduciton_B, self).__init__()
        self.branch_0 = nn.Sequential(
            Conv1d(in_channels, 256, 1, stride=1, padding=0, bias=True),
            Conv1d(256, 384, 3, stride=2, padding=0, bias=True)
        )
        self.branch_1 = nn.Sequential(
            Conv1d(in_channels, 256, 1, stride=1, padding=0, bias=True),
            Conv1d(256, 288, 3, stride=2, padding=0, bias=True),
        )
        self.branch_2 = nn.Sequential(
            Conv1d(in_channels, 256, 1, stride=1, padding=0, bias=True),
            Conv1d(256, 288, 3, stride=1, padding=1, bias=True),
            Conv1d(288, 320, 3, stride=2, padding=0, bias=True)
        )
        self.branch_3 = nn.MaxPool1d(3, stride=2, padding=0)

    def forward(self, x):
        x0 = self.branch_0(x)
        x1 = self.branch_1(x)
        x2 = self.branch_2(x)
        x3 = self.branch_3(x)
        return torch.cat((x0, x1, x2, x3), dim=1)


class Inception_ResNet_C(nn.Module):
    def __init__(self, in_channels, scale=1.0, activation=True):
        super(Inception_ResNet_C, self).__init__()
        self.scale = scale
        self.activation = activation
        self.branch_0 = Conv1d(in_channels, 192, 1, stride=1, padding=0, bias=True)
        self.branch_1 = nn.Sequential(
            Conv1d(in_channels, 192, 1, stride=1, padding=0, bias=True),
            Conv1d(192, 224, 3, stride=1, padding=1, bias=True),
            Conv1d(224, 256, 1, stride=1, padding=0, bias=True)
        )
        self.conv = nn.Conv1d(448, 2080, 1, stride=1, padding=0, bias=True)
        self.relu = nn.ReLU(inplace=True)
    def forward(self, x):
        x0 = self.branch_0(x)
        x1 = self.branch_1(x)
        x_res = torch.cat((x0, x1), dim=1)
        x_res = self.conv(x_res)
        if self.activation:
            return self.relu(x + self.scale * x_res)
        return x + self.scale * x_res


class Inception_ResNetv2(nn.Module):
    def __init__(self, in_channels=3, classes=2, k=256, l=256, m=384, n=384):
        super(Inception_ResNetv2, self).__init__()
        blocks = []
        blocks.append(Stem(in_channels))
        for i in range(3):
            blocks.append(Inception_ResNet_A(320, 0.1)) #[b, 320, 35]
        blocks.append(Reduction_A(320, k, l, m, n)) #[b, 1088, 17]
        for i in range(6):
            blocks.append(Inception_ResNet_B(1088, 0.1))
        blocks.append(Reduciton_B(1088))
        for i in range(4):
            blocks.append(Inception_ResNet_C(2080, 0.1))
        blocks.append(Inception_ResNet_C(2080, activation=False))
        self.features = nn.Sequential(*blocks)
        self.conv = Conv1d(2080, 1536, 1, stride=1, padding=0, bias=True)
        self.global_average_pooling = nn.AdaptiveAvgPool1d(1)
        self.linear = nn.Linear(1536, 769)

    def forward(self, x):
        x = self.features(x)
        x = self.conv(x)
        x = self.global_average_pooling(x)
        x = x.view(x.size(0), -1)
        x = self.linear(x)
        return x