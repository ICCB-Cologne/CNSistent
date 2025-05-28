import torch
import torch.nn as nn
import numpy as np
from torchview import draw_graph

# Your model definition
class MixedClassifier(nn.Module):
    def _calc_exponent(self, in_size, out_size):
        return (np.log(out_size) / np.log(in_size)) ** (1 / 2)

    def __init__(self, in_c, in_f, out_size):
        super(MixedClassifier, self).__init__()
        kernel = 5
        pad = kernel // 2
        pool = 2
        channels1 = 10
        channels2 = 5
        self.conv1 = nn.Conv1d(in_channels=in_c, out_channels=in_c*channels1, kernel_size=kernel, padding=pad)
        self.pool1 = nn.MaxPool1d(kernel_size=pool)
        self.bnc = nn.BatchNorm1d(in_c*channels1)
        self.conv2 = nn.Conv1d(in_channels=in_c*channels1, out_channels=in_c*channels2, kernel_size=kernel, padding=pad)
        self.pool2 = nn.MaxPool1d(kernel_size=pool)  
        self.bn1 = nn.BatchNorm1d(in_c*channels2)
        self.dropout = nn.Dropout() 
        flat_size = in_c*channels2 * (in_f // (pool ** 2))
        exp = self._calc_exponent(flat_size, out_size)
        next_size = int(flat_size ** exp)
        print(f"Layer sizes: {flat_size} -> {next_size} -> {out_size}")
        self.fc1 = nn.Linear(in_features=flat_size, out_features=next_size)
        self.bn2 = nn.BatchNorm1d(next_size)
        self.fc2 = nn.Linear(in_features=next_size, out_features=out_size)

    def forward(self, x): 
        x = self.pool1(torch.relu(self.conv1(x)))
        x = self.pool2(torch.relu(self.conv2(x)))
        x = self.bn1(x)  
        x = self.dropout(x) 
        x = torch.flatten(x, 1)
        x = torch.relu(self.fc1(x))
        x = self.bn2(x)
        x = self.dropout(x)
        x = self.fc2(x)
        return x

# Example model and input shape
in_c = 2        # e.g. number of channels
in_f = 128      # e.g. sequence length
out_size = 4    # e.g. number of classes

model = MixedClassifier(in_c, in_f, out_size)

# TorchView visualization
input_tensor = torch.zeros(1, in_c, in_f)
graph = draw_graph(
    model,
    input_data=input_tensor,
    graph_name="MixedClassifier",
    save_graph=True,              # Set to True to save as file
    directory="./",               # Output directory
    roll=True                     # Expand submodules inline
)
# To display in Jupyter, just:
# graph.visual_graph

# To save manually:
graph.visual_graph.render("mixed_classifier", format="png")
