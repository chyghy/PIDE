import torch.nn as nn
import torch.nn.functional as F
import torch
import esm

class ESMORF_binary(nn.Module):
    def __init__(self, n_class, lis = ["32","31","30","29"]): # 调整4层
        super(ESMORF_binary, self).__init__()
        self.esm = esm.pretrained.esm2_t33_650M_UR50D()[0]
        for name, param in self.named_parameters():
            sp = name.split('.')
            if len(set(lis) & set(sp)) != 0:
                param.requires_grad =True
            else:
                param.requires_grad =False
        self.fc1 = nn.Linear(1280, 960)
        self.fc2 = nn.Linear(960, 480)
        self.fc3 = nn.Linear(480, 120)
        self.fc4 = nn.Linear(120, 30)
        self.fc5 = nn.Linear(30, n_class)
    def forward(self, batch_tokens):
        x = self.esm(batch_tokens, repr_layers=[33], return_contacts=False)["representations"][33]
        batch_tokens = batch_tokens.unsqueeze(-1)
        x = x.masked_fill(batch_tokens==2, 0)
        x = x.masked_fill(batch_tokens==1, 0)[:, 1:, :]
        num = torch.sum(batch_tokens>2, axis=1)
        x = x.sum(axis=1) / num
        ret = []
        ret.append(x)
        x = F.relu(self.fc1(x))
        ret.append(x)
        x = F.relu(self.fc2(x))
        ret.append(x)
        x = F.relu(self.fc3(x))
        ret.append(x)
        x = F.relu(self.fc4(x))
        ret.append(x)
        ret.append(self.fc5(x))        
        return ret



class MLP_esm650(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden1 = nn.Linear(1280,960)
        self.hidden2 = nn.Linear(960,400)
        self.hidden3 = nn.Linear(400,100)
        self.hidden4 = nn.Linear(100,20)
        self.out = nn.Linear(20,2)
    
    def forward(self, X):
        x = F.relu(self.hidden1(X))
        x = F.relu(self.hidden2(x))
        x = F.relu(self.hidden3(x))
        x = F.relu(self.hidden4(x))
        net = self.out(x)
        return net



