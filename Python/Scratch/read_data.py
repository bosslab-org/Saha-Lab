import gzip
import os
import numpy as np
import matplotlib.pyplot as plt

readpath = '/Users/alexanderfarnum/Documents/Code/Datasets/'
mnist = 'mnist'
cifar = 'cifar-10-batches/'
np.random.seed(10)

def one_hot(y):
    one_hot_y = np.zeros((y.shape[1], np.max(y) + 1))
    one_hot_y[np.arange(y.shape[1]), y] = 1
    one_hot_y = one_hot_y.T
    return one_hot_y


def flatten(x):
    return np.reshape(x, (x.shape[0], x.shape[1] * x.shape[2]))


def read_mnist(val_size=0.2):
    """Import and preprocess MNIST dataset.
    Constructs one hot vector for labels.
    Accepts optional validation set size, else val_size = 0.2
    """
    X_train = 'train-images-idx3-ubyte.gz'
    y_train = 'train-labels-idx1-ubyte.gz'
    X_test = 't10k-images-idx3-ubyte.gz'
    y_test = 't10k-labels-idx1-ubyte.gz'

    img_size = 28  # 28 x 28 pixel images
    train_samples = 60000
    test_samples = 10000

    X_train_path = readpath + mnist + '/' + X_train
    y_train_path = readpath + mnist + '/' + y_train
    X_test_path = readpath + mnist + '/' + X_test
    y_test_path = readpath + mnist + '/' + y_test

    print('Extracting X_train')
    with gzip.open(X_train_path) as temp:
        temp.read(16)
        temp = temp.read(img_size * img_size * train_samples)
        X_train = np.frombuffer(temp, dtype=np.uint8).astype(np.int64)
        X_train = X_train.reshape(train_samples, img_size*img_size).T
    print('Extracting y_train')
    with gzip.open(y_train_path) as temp:
        temp.read(8)
        temp = temp.read(1 * train_samples)
        y_train = np.frombuffer(temp, dtype=np.uint8).astype(np.int64)
        y_train = y_train.reshape(train_samples, 1).T
    print('Extracting X_test')
    with gzip.open(X_test_path) as temp:
        temp.read(16)
        temp = temp.read(img_size * img_size * test_samples)
        X_test = np.frombuffer(temp, dtype=np.uint8).astype(np.float32)
        X_test = X_test.reshape(test_samples, img_size*img_size).T
    print('Extracting y_test')
    with gzip.open(y_test_path) as temp:
        temp.read(8)
        temp = temp.read(1 * test_samples)
        y_test = np.frombuffer(temp, dtype=np.uint8).astype(np.int64)

    X_val, X_train = np.split(X_train, [int(val_size*X_train.shape[1])], axis=1)
    y_val, y_train = np.split(y_train, [int(val_size*y_train.shape[1])], axis=1)

    return X_train/255, y_train, X_val/255, y_val, X_test/255, y_test


def read_cifar(val_size=0.2):
    """Import and preprocess CIFAR-10 dataset.
    Accepts optional validation size, else val_size = 0.2
    """
    source = 'https://www.cs.toronto.edu/~kriz/cifar.html'
    test_readpath = 'testdata'
    def unpickle(file):
        import pickle
        with open(file, 'rb') as fo:
            dict = pickle.load(fo, encoding='latin1')
        return dict
    X_train = []
    y_train = []
    for train_filepath in os.listdir(readpath + cifar):
        if train_filepath.startswith("traindata"):
            data_temp = unpickle(readpath + cifar + train_filepath)
            X_train.append(data_temp['data'])
            y_train.append(data_temp['labels'])
    X_train = np.array(X_train)
    X_train = X_train.reshape(X_train.shape[0]*X_train.shape[1], X_train.shape[2])
    y_train = np.array(y_train)
    y_train = y_train.reshape(y_train.shape[0]*y_train.shape[1])

    np.random.shuffle(X_train)
    np.random.shuffle(y_train)

    X_val, X_train = np.split(X_train, [int(val_size*X_train.shape[0])])
    y_val, y_train = np.split(y_train, [int(val_size*y_train.shape[0])])

    data_temp = unpickle(readpath + cifar + test_readpath)
    X_test = np.array(data_temp['data'])
    y_test = np.array(data_temp['labels'])

    return X_train/255, y_train, X_val/255, y_val, X_test/255, y_test


def read_locust_data():
    import scipy.io  # v7
    data = scipy.io.loadmat('/Users/XanderFarnum/Documents/NN/Sandbox/Datasets/Locust/Cell_Culture_RMS_Master.mat')

    # import h5py  # v7.3
    # f = h5py.File('somefile.mat', 'r')
    # data = f.get('data/variable1')
    # data = np.array(data)  # For converting to a NumPy array

    values = [v for k, v in data.items() if k.endswith('RMS_data_filt')]
    # data_dict = {keys: values for (key, value) in data.items() if key.endswith('RMS_data_filt')}

    num_neurons = values[0][:, 0, 0, 0].size
    num_trials = values[0][0, 0, :, 0].size
    num_samples = values[0][0, 0, 0, :].size
    num_stimuli = len(values)

    data = np.zeros((num_neurons, num_trials, num_samples, num_stimuli))
    labels = np.zeros((num_neurons, num_stimuli))
    for i in range(len(values)):
        data[:, :, :, i] = values[i].mean(axis=1)
        labels[:, i] = i

    data = data.transpose((0, 3, 2, 1)).reshape((num_neurons * num_stimuli, num_samples, num_trials)) / data.max()

    # data = data.transpose((0, 1, 3, 2)).reshape((num_neurons * num_trials * num_stimuli, num_samples)) / data.max()
    labels = labels.reshape((num_neurons * num_stimuli, 1))

    # np.random.shuffle(data)
    # np.random.shuffle(labels)
    # X_test, X_train = np.split(data, [int(val_size*data.shape[0])])
    # y_test, y_train = np.split(labels, [int(val_size*labels.shape[0])])

    return data, labels


def spiral_data(n, m, k):
    """Constructs spiral dataset.
    Number of samples per class, number of dimensions, number of classes.
    """
    X_train = np.zeros((n*k, m))
    y_train = np.zeros(n*k, dtype='uint8')
    for j in range(k):
        ix = range(n*j, n*(j+1))
        r = np.linspace(0.0, 1, n)  # radius
        t = np.linspace(j*4, (j+1)*4, n) + np.random.randn(n)*0.2  # theta
        X_train[ix] = np.c_[r*np.sin(t), r*np.cos(t)]
        y_train[ix] = j
    fig = plt.figure()
    plt.scatter(X_train[:, 0], X_train[:, 1], c=y_train, s=40, cmap=plt.cm.Spectral)
    plt.xlim([-1, 1])
    plt.ylim([-1, 1])

    return X_train, y_train, fig
