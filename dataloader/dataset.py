import torch
from torch.utils.data import Dataset, DataLoader

class RobotData(Dataset):
    def __init__(self, np_data, mean, std):
        self.data = []
        for i in range(np_data.shape[0]):
            alpha = (np_data[i, 0] - mean[0])/std[0]
            compressL = (np_data[i, 1] - mean[1])/std[1]
            mu = (np_data[i, 2] - mean[2])/std[2]
            rho = (np_data[i, 4] - mean[4])/std[4]

            y = np_data[i, 6]/0.02
            x = np_data[i, 7]/0.02

            features = [alpha, compressL, mu, rho]
            output = [x, y]

            self.data.append((torch.Tensor(features), torch.tensor(output)))

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]


