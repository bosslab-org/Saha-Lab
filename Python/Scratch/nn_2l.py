import numpy as np
import sys

sys.path.insert(0,'/Users/alexanderfarnum/Documents/Code/Python/ML')
from read_data import read_mnist, one_hot

x_train, y_train, x_val, y_val, x_test, y_test = read_mnist()
one_hot_y = one_hot(y_train)
# x_train = [784, num_samples=48000]
# y_train = [num_samples=48000,1]
# one_hot_y = [10, num_samples=480000]

alpha = 0.1

w1 = np.random.rand(10,784)  # w = [10, 784]
b1 = np.random.rand(10,1)    # b = [10, 1]
w2 = np.random.rand(10,10)   # w = [10, 10]
b2 = np.random.rand(10,1)    # b = [10, 1]

for i in range(500):
    # forward- h1
    z1 = np.dot(w1, x_train) + b1
    a1 = np.maximum(z1, 0)
    # forward- h2
    z2 = np.dot(w2, a1) + b2
    a2 = np.exp(z2) / sum(np.exp(z2))
    # compute total loss

    # compute h2 gradients- z, w, b
    dz2 = a2 - one_hot_y
    dw2 = 1 / x_train.shape[1] * dz2.dot(a1.T)
    db2 = 1 / x_train.shape[1] * np.sum(dz2)
    # compute h1 gradients- z, w, b
    dz1 = w2.T.dot(dz2) * (z1 > 0)
    dw1 = 1 / x_train.shape[1] * dz1.dot(x_train.T)
    db1 = 1 / x_train.shape[1] * np.sum(dz1)
    # update parameters
    w1 -= alpha * dw1
    b1 -= alpha * db1    
    w2 -= alpha * dw2  
    b2 -= alpha * db2  
    if i % 10 == 0:
        y_hat = np.argmax(a2, 0)
        print("Total accuracy: ", round(np.sum(y_hat == y_train) / y_train.size * 100,2), "%")
