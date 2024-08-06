import torch
import torch.nn as nn
from torch.nn.init import orthogonal_

class RNN(nn.Module):
    def __init__(self, args):
        super(RNN, self).__init__()
        input_dim = args["model_param"]["input_dim"]
        output_dim = args["model_param"]["output_dim"]
        layer_sizes = args["model_param"]["layer_sizes"]
        dropout = args["model_param"]["dropout"]
        # define the input module
        self.input_block = nn.Sequential(
            nn.Linear(input_dim, layer_sizes[0]),
            nn.ReLU()
        )
        # the hidden layers
        self.hidden_block = nn.Sequential()
        for i in range(len(layer_sizes) - 1):
            temp_layers = nn.Sequential(
                nn.Linear(layer_sizes[i], layer_sizes[i+1]),
                nn.ReLU()
            )

            self.hidden_block.add_module(f"hidden({i})",
                                         temp_layers)
        # output block
        self.output = nn.Sequential(
            nn.Linear(layer_sizes[-1], output_dim)
        )

    def get_num_params(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad) 

    def forward(self, x):
        out1 = self.input_block(x)
        out2 = self.hidden_block(out1)
        out3 = self.output(out2)

        return out3
