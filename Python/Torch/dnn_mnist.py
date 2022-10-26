import torch
import torch.nn as nn
import torchvision
import torchvision.transforms as transforms


class Model(nn.Module):
    def __init__(self, input_size, h1, h2, num_classes):
        super().__init__()
        self.input_size = input_size
        self.l1 = nn.Linear(input_size, h1) 
        self.relu = nn.ReLU()
        self.l2 = nn.Linear(h1, h2)
        self.relu = nn.ReLU()  
        self.l3 = nn.Linear(h2, num_classes)  
    
    def forward(self, x):
        out = self.l1(x)
        out = self.relu(out)
        out = self.l2(out)
        out = self.relu(out)
        out = self.l3(out)
        return out


input_size = 784
h1 = 100
h2 = 50
num_classes = 10
num_epochs = 4
batch_size = 100
learning_rate = 0.001

train_dataset = torchvision.datasets.MNIST(root='/Users/alexanderfarnum/Documents/Code/Datasets', 
                                           train=True, 
                                           transform=transforms.ToTensor(),  
                                           download=True)

test_dataset = torchvision.datasets.MNIST(root='/Users/alexanderfarnum/Documents/Code/Datasets', 
                                          train=False, 
                                          transform=transforms.ToTensor())

train_loader = torch.utils.data.DataLoader(dataset=train_dataset, 
                                           batch_size=batch_size, 
                                           shuffle=True)

test_loader = torch.utils.data.DataLoader(dataset=test_dataset, 
                                          batch_size=batch_size, 
                                          shuffle=False)


# Model instantiation
model = Model(input_size, h1, h2, num_classes)

# Select loss
criterion = nn.CrossEntropyLoss()

# Select optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)  

# Training
for epoch in range(num_epochs):
    for i, (x_train, y_train) in enumerate(train_loader):  
        x_train = x_train.reshape(-1, 28*28)
        
        outputs = model(x_train)
        loss = criterion(outputs, y_train)
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if (i+1) % 50 == 0:
            print (f'Epoch {epoch+1}/{num_epochs}, Iteration {i+1}/{len(train_loader)}, Loss: {loss.item():.4f}')

# Testing
with torch.no_grad():
    test_correct = 0
    test_total = 0
    for x_test, y_test in test_loader:
        x_test = x_test.reshape(-1, 28*28)

        outputs = model(x_test)
        _, predicted = torch.max(outputs.data, 1)
        test_total += y_test.size(0)
        test_correct += (predicted == y_test).sum().item()

    print(f'Test accuracy: {100.0 * test_correct / test_total} %')
