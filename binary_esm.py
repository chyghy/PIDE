import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
import torch.nn.functional as F
import random
import model
import sys

def binary_test(input_data,esmorf,device,output_data,output_data2,batch_size):

    # 定义一个数据加载器
    class ORFData(Dataset):
        def __init__(self):
            super(ORFData, self).__init__()
            self.x = np.array(input_data)       
        def __getitem__(self,idx):
            x = torch.tensor(self.x[idx], dtype=torch.long)
            return x
        def __len__(self):
            return len(self.x)


    test_dataset = ORFData()
    test_dataloader = DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=False)
    all_predict = np.array([])
    all_max = np.array([])

    with torch.no_grad():
        esmorf = esmorf.eval()
        for batch_id, batch_data in enumerate(test_dataloader):
            inputs = batch_data.to(device)
            predicts = F.softmax(esmorf(inputs)[5],1).cpu().detach().numpy()
            all_predict = np.concatenate((all_predict, predicts[:, 1]))
            all_max = np.concatenate((all_max, np.argmax(predicts, axis=1)))
        Predict = pd.DataFrame(all_predict)
        Max = pd.DataFrame(all_max)
        Predict.to_csv(output_data)
        Max.to_csv(output_data2)

