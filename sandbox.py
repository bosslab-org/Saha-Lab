import numpy as np
import random
from matplotlib import pyplot as plt
from Read_Data import read_mnist, read_cifar, spiral_data, read_locust_data


# Design matrix = [m data samples , n features]
# X_train, y_train, X_val, y_val, X_test, y_test = read_mnist()
# X_train, y_train, X_val, y_val, X_test, y_test = read_cifar()
# X_train, y_train, spiral_fig = spiral_data(100, 2, 3)

data, labels = read_locust_data()

print(data.shape)
print(labels.shape)
print(data.shape[2])
print(data[3].shape)

# for i in range(data.shape[2]):
#     X_train = data(:, :, i)
#     X_test =


alpha = 1e-0
reg = 1e-3  # regularization strength


def single_layer_perceptron(X, y, num_classes, alpha, reg, num_epochs=10000, fig=0):
    W = 0.01 * np.random.randn(X.shape[1], num_classes)
    b = np.zeros((1, num_classes))
    for i in range(num_epochs):
        # Forward propagation
        scores = np.dot(X, W) + b
        # Normalize scores to avoid overflow error (min = 0, max = 1)
        # scores = (scores - np.min(scores)) / (np.max(scores) - np.min(scores))
        # Each row of probs contains the class probabilities
        probs = np.exp(scores) / np.sum(np.exp(scores), axis=1, keepdims=True)

        correct_logprobs = -np.log(probs[range(X.shape[0]), y])
        # compute data loss: average cross-entropy loss
        data_loss = np.sum(correct_logprobs)/X.shape[0]
        # compute regularization loss
        reg_loss = 0.5*reg*np.sum(W*W)
        # compute total loss
        loss = data_loss + reg_loss
        if i % (num_epochs/10) == 0:
            print("iteration %d: loss %f" % (i, loss))

        probs[range(X.shape[0]), y] -= 1
        probs /= X.shape[0]

        dW = np.dot(X.T, probs)
        db = np.sum(probs, axis=0, keepdims=True)
        # Regularization gradient
        dW += reg*W

        # Update parameters
        W += -alpha * dW
        b += -alpha * db

    # evaluate training set accuracy
    scores = np.dot(X, W) + b
    predicted_class = np.argmax(scores, axis=1)
    print('training accuracy: %.2f' % (np.mean(predicted_class == y)))

    if fig == 1:
        h = 0.02
        x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
        y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
        xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
        z = np.dot(np.c_[xx.ravel(), yy.ravel()], W) + b
        z = np.argmax(z, axis=1)
        z = z.reshape(xx.shape)
        plt.contourf(xx, yy, z, cmap=plt.cm.Spectral, alpha=0.3)
        plt.scatter(X[:, 0], X[:, 1], c=y, s=30, cmap=plt.cm.Spectral, edgecolors=[0, 0, 0, 1])
        plt.xlim(xx.min(), xx.max())
        plt.ylim(yy.min(), yy.max())
        plt.show()
    return X, y, W, b


def single_layer_neural_network(X, y, num_classes, alpha, reg, h_size=100, batch_size=256, num_epochs=10000, fig=0):
    W1 = 0.01 * np.random.randn(X.shape[1], h_size)
    b1 = np.zeros((1, h_size))
    W2 = 0.01 * np.random.randn(h_size, num_classes)
    b2 = np.zeros((1, num_classes))
    for i in range(num_epochs):
        # Delineate minibatch samples
        batch_idx = random.sample(range(X.shape[0]), batch_size)
        X_batch = X[batch_idx, :]
        y_batch = y[batch_idx]

        # Forward propagation with ReLU
        hidden_layer = np.maximum(0, np.dot(X_batch, W1) + b1)
        scores = np.dot(hidden_layer, W2) + b2

        # Each row of probs contains the class probabilities
        probs = np.exp(scores) / np.sum(np.exp(scores), axis=1, keepdims=True)

        correct_logprobs = -np.log(probs[range(X_batch.shape[0]), y_batch])
        # compute data loss: average cross-entropy loss
        data_loss = np.sum(correct_logprobs)/X_batch.shape[0]
        # compute regularization loss
        reg_loss = 0.5 * reg * np.sum(W1 * W1) + 0.5 * reg * np.sum(W2 * W2)
        # compute total loss
        loss = data_loss + reg_loss
        if i % (num_epochs/10) == 0:
            # print("iteration ", i)
            print("iteration %d: loss %f" % (i, loss))

        probs[range(X_batch.shape[0]), y_batch] -= 1
        probs /= X_batch.shape[0]

        dW2 = np.dot(hidden_layer.T, probs)
        db2 = np.sum(probs, axis=0, keepdims=True)
        dhidden = np.dot(probs, W2.T)
        dhidden[hidden_layer <= 0] = 0
        dW1 = np.dot(X_batch.T, dhidden)
        db1 = np.sum(dhidden, axis=0, keepdims=True)

        # Regularization gradient
        dW2 += reg * W2
        dW1 += reg * W1

        # Update parameters
        W1 += -alpha * dW1
        b1 += -alpha * db1
        W2 += -alpha * dW2
        b2 += -alpha * db2

    # evaluate training set accuracy
    hidden_layer = np.maximum(0, np.dot(X, W1) + b1)
    scores = np.dot(hidden_layer, W2) + b2
    predicted_class = np.argmax(scores, axis=1)
    print('training accuracy: %.2f' % (np.mean(predicted_class == y)))

    if fig == 1:
        h = 0.02
        x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
        y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
        xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
        z = np.dot(np.maximum(0, np.dot(np.c_[xx.ravel(), yy.ravel()], W1) + b1), W2) + b2
        z = np.argmax(z, axis=1)
        z = z.reshape(xx.shape)
        plt.contourf(xx, yy, z, cmap=plt.cm.Spectral, alpha=0.3)
        plt.scatter(X[:, 0], X[:, 1], c=y, s=30, cmap=plt.cm.Spectral, edgecolors=[0, 0, 0, 1])
        plt.xlim(xx.min(), xx.max())
        plt.ylim(yy.min(), yy.max())
        plt.show()
    return X, y, W1, b1, W2, b2


# X, y, W, b = single_layer_perceptron(X_train, y_train, 10, alpha, reg, num_epochs=10, fig=1)
# X, y, W1, b1, W2, b2 = single_layer_neural_network(X_train, y_train, 10, alpha, reg, batch_size=300, num_epochs=10000, fig=1)

# int(X_train.shape[0]/10)