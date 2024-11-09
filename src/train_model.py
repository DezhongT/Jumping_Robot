import sys
import os
import numpy as np
import yaml
from config import configer
import logging
from dataloader import dataset
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch
import torch.optim as optim
import matplotlib.pyplot as plt
from os import path as osp

def save_model(path, epoch, network, optimizer):
    model_path = osp.join(path, "checkpoints", "checkpoint_%d.pt" % epoch)
    if not osp.isdir(osp.join(path, "checkpoints")):
        os.makedirs(osp.join(path, "checkpoints"))

    state_dict = {
        "model_state_dict": network.state_dict(),
        "epoch": epoch,
        "optimizer_state_dict": optimizer.state_dict(),
    }
    torch.save(state_dict, model_path)
    logging.info(f"Model saved to {model_path}")

def create_output_dir():
    out_dir = "./output"
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.isdir(os.path.join(out_dir, "checkpoints")):
        os.makedirs(os.path.join(out_dir, "checkpoints"))
    if not os.path.isdir(os.path.join(out_dir, "logs")):
        os.makedirs(os.path.join(out_dir, "logs"))
    logging.info(f"Creating the training output folder: {out_dir}")

def train_load_data(batch_size, train_ratio, data_name):
    np_data = np.loadtxt(data_name)
    np_data = np_data.astype(np.float32)
    idx = (np_data[:, 6] != 0)
    np_data = np_data[idx]
    mean = np.mean(np_data, axis = 0)
    std = np.std(np_data, axis = 0)

    indices = list(range(np_data.shape[0]))
    np.random.shuffle(indices)
    index = int(train_ratio * np_data.shape[0])
    train_data = dataset.RobotData(np_data[:index], mean, std)
    test_data = dataset.RobotData(np_data[index:], mean, std)

    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=batch_size, shuffle=False)

    return train_loader, test_loader

yaml_path = "model.yaml"
with open(yaml_path, 'r') as file:
    args = yaml.safe_load(file)

# build model
model = configer.build_model(args)
# create a output dir
create_output_dir()

logging.info(f"Loading training data")
data_path = "./data/train_data.txt"
batch_size = 1024
train_ratio = 0.8

learning_rate = 1e-4
weight_decay = 0.0
criterion = nn.L1Loss()
optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

num_epochs = 200
out_dir = "./output/"

train_losses = []
test_losses = []

# train with gradually increased batch size
batch_size = 64 / 2
for epoch in range(num_epochs):
    model.train()
    running_loss = 0.0
    if epoch % 50 == 0:
        batch_size = int(batch_size * 2)
        train_loader, val_loader = train_load_data(batch_size, train_ratio, data_path)

    for i, (x_train, y_train) in enumerate(train_loader):
        x_train = x_train.to(device)
        y_train = y_train.to(device)
        y_pred = model(x_train)
        loss = criterion(y_pred, y_train)
        optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)  # Adjust max_norm as needed

        optimizer.step()
        running_loss += loss.item()

    epoch_train_loss = running_loss / len(train_loader)
    train_losses.append(epoch_train_loss)

    model.eval()
    running_loss  = 0.0

    with torch.no_grad():
        for x_val, y_val in val_loader:
            x_val = x_val.to(device)
            y_val = y_val.to(device)
            y_pred = model(x_val)
            loss = criterion(y_val, y_pred)
            running_loss += loss.item()
    
    # Calculate test loss
    epoch_test_loss = running_loss / len(val_loader)
    test_losses.append(epoch_test_loss)

    if (epoch + 1)%10 == 0:
        save_model(out_dir, epoch+1, model, optimizer)

    print(f'Epoch [{epoch+1}/{num_epochs}], '
          f'Train Loss: {epoch_train_loss:.7f}, '
          f'Test Loss: {epoch_test_loss:.7f}')
    plt.clf()
    plt.plot(train_losses, label='Training Loss')
    plt.plot(test_losses, label='Test Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Loss Curve')
    plt.legend()
    plt.yscale('log')
    plt.draw()
    plt.pause(0.1)

plt.show()
