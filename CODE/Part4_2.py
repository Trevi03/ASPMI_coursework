import numpy as np
import torch
import torch.nn as nn
from sklearn.utils import shuffle
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.io
import random

import os

# Get the directory of the current script
script_directory = os.path.dirname(os.path.abspath(__file__))
# Change to that directory
os.chdir(script_directory)

data_dct = scipy.io.loadmat('example_training_n2500_1freq_rescale_k8.mat')

#extract training data and labels
trainX = data_dct['trainX']
trainy = data_dct['trainy']
trainX, trainy = shuffle(trainX, trainy)

#define training and model parameters
num_epochs = 200
batch_size = 30
n_in, n_out = 1, 8
kernel_size = 500
padding = 0
dilation = 1
stride = 1
learning_rate = 0.000001
L_in = trainX.shape[-1]
trainX = trainX.reshape((trainX.shape[0], 1, L_in))

# calculate the output length after convolution
L_out = int((L_in+2*padding-dilation*(kernel_size-1)-1)/stride+1)

print("Defining model")
# define the model using pytorch
class ConvNet1D(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer1 = nn.Conv1d(n_in, n_out, kernel_size=kernel_size,padding=padding)
        self.layer2 = nn.MaxPool1d(L_out)
        # self.layer2 = nn.AvgPool1d(L_out)
        self.act = nn.Softmax()

    def forward(self, x):
        out = self.layer1(x)
        out = torch.abs(out)
        out = self.layer2(out)
        out = self.act(out)
        return out
    
model = ConvNet1D()
# Loss and optimizer
criterion = nn.L1Loss()#nn.BCEWithLogitsLoss()#nn.L1Loss()#nn.MSELoss()#.BCEWithLogitsLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

print('Training model')
# Train the model
total_step = trainX.shape[0]

# transformation of data into torch tensors
trainXT = torch.from_numpy(trainX.astype('float32'))
trainyT = torch.from_numpy(trainy.astype('float32'))

##### create initialised weights for each of the 8 kernels
Nt = kernel_size
t = np.linspace(0, 2 * np.pi, kernel_size)
layer1_init = np.reshape(np.zeros((1, kernel_size)),(1,1,kernel_size))
layer2_init = np.reshape(np.cos(0.5 * t), (1,1,kernel_size))   # cos(0.5t)
layer3_init = np.reshape(np.cos(1.0 * t), (1,1,kernel_size))   # cos(t)
layer4_init = np.reshape(np.cos(1.5 * t), (1,1,kernel_size))   # cos(1.5t)
layer5_init = np.reshape(np.cos(2.0 * t), (1,1,kernel_size))   # cos(2.0t)
layer6_init = np.reshape(np.cos(2.5 * t), (1,1,kernel_size))   # cos(2.5t)
layer7_init = np.reshape(np.cos(3.0 * t), (1,1,kernel_size))   # cos(3.0t)
layer8_init = np.reshape(np.cos(3.5 * t), (1,1,kernel_size))   # cos(3.5t)
weights_init = np.concatenate((layer1_init,layer2_init,layer3_init,layer4_init,layer5_init,layer6_init,layer7_init,layer8_init),axis=0)
weights_initT = torch.from_numpy(weights_init.astype('float32'))
model.layer1.weight.data = weights_initT

loss_list = []
for epoch in range(num_epochs):
    correct_sum = 0
    correct_sum_test = 0

    for i in range(total_step//batch_size): # split data into batches
        trainXT_seg = trainXT[i*batch_size:(i+1)*batch_size, :, :]
        trainyT_seg = trainyT[i*batch_size:(i+1)*batch_size, :]

        # Run the forward pass
        outputs = model(trainXT_seg)

        #ensure outputs are the correct shape
        outputs = outputs.reshape((outputs.shape[0], 8))

        #calculate loss for the batch
        loss = criterion(trainyT_seg,outputs)

        # Backprop and perform Adam optimisation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    output_max = torch.argmax(outputs, 1, keepdim=False)
    reference_max = torch.argmax(trainyT_seg, 1, keepdim=False)
    count = torch.count_nonzero((output_max-reference_max)).cpu().detach()
    accuracy = 100*((batch_size-count)/batch_size)

    print("Epoch")
    print(epoch)
    print("Train loss")
    print(loss)
    print("Accuracy (%)")
    print(accuracy)
    loss_list.append(loss.item())

#plot the training accuracy per epoch
# plt.title('Training loss: AvgPool1d')
plt.title('Training loss: MaxPool1d')
plt.ylabel('L1 Loss')
plt.xlabel('Epoch')
plt.grid(True)
plt.autoscale(axis='x', tight=True)
plt.plot(np.log(loss_list), label = "Training loss")
plt.legend()
plt.show()

#### plot the trained kernel weights
kernels = model.layer1.weight.T.cpu().detach().numpy()[:,0,:]
fig, axs = plt.subplots(kernels.shape[1],1, figsize=(5,8), squeeze=True)
for i, ax in enumerate(axs):
    ax.plot(kernels[:,i])
plt.show()