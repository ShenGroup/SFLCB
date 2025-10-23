import numpy as np
import cvxpy as cp
import time
from torch.autograd import grad 


import matplotlib.pyplot as plt

import torch
from torch import nn 
from torch.nn import functional as F

import sys
sys.path.append('..')

device = torch.device("cpu")

from utils import load_diabetes, train_val_test_split


def Loss_upper(w, b, C_tensor_val, y_val, x_val):
    x = torch.reshape(y_val, (y_val.shape[0],1)) 
    x = x* F.linear(x_val, w, b)
    #x = x* (x_val @ w1 +b1)
    return torch.sum(torch.exp(1-x)) + 0.5 * torch.linalg.norm(C_tensor_val)**2

def gaffa(x_train, y_train, x_val, y_val, x_test, y_test, hparams, epochs, early_stopping_th = False,verbose=True):
    feature=x_train.shape[1] # = 8
    nTr=y_train.shape[0]

    # Dataset to tensor
    y_train = torch.tensor(y_train).float()
    x_train = torch.tensor(x_train).float()
    y_val = torch.tensor(y_val).float()
    x_val = torch.tensor(x_val).float()
    y_test = torch.tensor(y_test).float()
    x_test = torch.tensor(x_test).float()
    
    
    Loss_lower = lambda w, b : .5 * (torch.sum(torch.square(w)))
    con_lower = lambda w, b, c: 1 - c - y_train * (F.linear(x_train, w, b).squeeze())
    relu0 = nn.ReLU()
    relu1 = nn.LeakyReLU()


    ###### Parameters
    p = hparams['p']
    gam1 = hparams['gam1']
    gam2 = hparams['gam2']
    eta = hparams['eta']
    alpha = hparams['alpha']
    beta = hparams['beta']
    ck0 = hparams['ck0']
    R = hparams['R']

    # For storage
    val_loss_list=[]
    test_loss_list=[]
    val_acc_list=[]
    test_acc_list=[]
    time_computation=[]
    algorithm_start_time=time.time()

    metrics = []
    variables = []

    # Locate Variable
    w0 = torch.ones(feature, requires_grad=True, device=device)
    b0 = torch.zeros(1, requires_grad=True, device=device)
    w1 = torch.ones(feature, requires_grad=True, device=device)
    b1 = torch.zeros(1, requires_grad=True, device=device)
    lam0 = torch.ones(nTr, device=device)
    lam1 = torch.ones(nTr, device=device)
    #c0 = -6. + 1. * torch.rand(nTr, device=device, requires_grad=True)
    c0 =  torch.Tensor(x_train.shape[0]).uniform_(1.,5.)

    # Set the Initial Value
    svm_w = cp.Variable(feature)
    svm_b = cp.Variable()
    loss = 1/2 *(cp.sum_squares(svm_w))
    constraints = [1 - c0.detach().numpy() - cp.multiply(y_train.cpu().numpy(), (x_train.cpu().numpy() @ svm_w + svm_b)) <= 0]
    prob = cp.Problem(cp.Minimize(loss), constraints)

    algorithm_start_time = time.time()
    # prob.solve(solver='ECOS')
    prob.solve(solver='ECOS', abstol=2e-3,reltol=2e-3,max_iters=1000000000, warm_start=True)  

    w0.data = torch.tensor(svm_w.value, device = device, dtype=torch.float32)
    w0=w0.reshape(1,feature)
    b0.data = torch.tensor(svm_b.value, device = device, dtype=torch.float32)
    w1.data = torch.clone(w0)
    b1.data = torch.tensor(svm_b.value, device = device, dtype=torch.float32)
    lam0.data = torch.tensor(constraints[0].dual_value, device = device, dtype=torch.float32)
    lam1.data = torch.tensor(constraints[0].dual_value, device = device, dtype=torch.float32)


    # print(f"{0:4d}-th iter: Val Loss = {ValLoss[-1]:.4f} Test Loss = {TestLoss[-1]: .4f} ValAcc = {ValAcc[-1]:.2f} Test Acc = {TestAcc[-1]:.2f}")


    for epoch in range(epochs):
        variables.append({
            'C': c0,
            'w': w0,
            'b': b0,
            'w_F': w1,
            'b_F': b1
        })

        x = torch.reshape(y_val, (y_val.shape[0],1)) 
        x = x* F.linear(x_val, w0, b0) 

        x1 = torch.reshape(y_test, (y_test.shape[0],1)) 
        x1 = x1 * F.linear(x_test, w0, b0) 
        test_loss_upper= torch.sum(torch.exp(1-x1))

        val_loss_F = (torch.sum(torch.exp(1-x))).detach().numpy()/y_val.shape[0]
        test_loss_F = test_loss_upper.detach().numpy()/y_test.shape[0]

        x = torch.reshape(y_val, (y_val.shape[0],1))
        x = x* F.linear(x_val, w1, b1) 


        x1 = torch.reshape(y_test, (y_test.shape[0],1)) 
        x1 = x1 * F.linear(x_test, w1, b1) 

        test_loss_upper= torch.sum(torch.exp(1-x1))

        val_loss = (torch.sum(torch.exp(1-x))).detach().numpy()/y_val.shape[0]
        test_loss = test_loss_upper.detach().numpy()/y_test.shape[0]

        loss_upper = val_loss
        loss_lower = (1/2) * (w0**2).sum()

        ###### Accuracy
        q = y_train * (w0 @ x_train.T + b0)
        train_acc = (q>0).sum() / len(y_train)

        q = y_val * (w0 @ x_val.T + b0)
        val_acc = (q>0).sum() / len(y_val)

        q = y_test * (w0 @ x_test.T + b0)
        test_acc = (q>0).sum() / len(y_test)

        q = y_train * (w1 @ x_train.T + b1)
        train_acc_F = (q>0).sum() / len(y_train)

        q = y_val * (w1 @ x_val.T + b1)
        val_acc_F = (q>0).sum() / len(y_val)

        q = y_test * (w1 @ x_test.T + b1)
        test_acc_F = (q>0).sum() / len(y_test)

        metrics.append({
            #'train_loss': train_loss,
            'train_acc': train_acc,
            'train_acc_F': train_acc_F,
            'val_loss': val_loss,
            'val_loss_F': val_loss_F,
            'val_acc': val_acc,
            'val_acc_F': val_acc_F,
            'test_loss': test_loss,
            'test_loss_F': test_loss_F,
            'test_acc': test_acc,
            'test_acc_F': test_acc_F,
            'loss_upper': loss_upper,
            'loss_lower': loss_lower,
            'time_computation': time.time()-algorithm_start_time
        })




        c0.requires_grad_(True)
        w0.requires_grad_(True)
        w1.requires_grad_(True)
        b0.requires_grad_(True)
        b1.requires_grad_(True)
        c0.grad=None
        w0.grad=None
        w1.grad=None
        b0.grad=None
        b1.grad=None



        ck = ck0 * (epoch + 1)**p
        f1k = Loss_lower(w1, b1)
        g1k = lam0 @ con_lower(w1, b1, c0)

        dw1 = grad(f1k, w1, retain_graph=True)[0] + grad(g1k, w1, retain_graph=True)[0] + 1 / gam1 * (w1 - w0)
        db1 = grad(g1k, b1, retain_graph=True)[0] + 1 / gam1 * (b1 - b0)

        w1p = w1 - eta * dw1
        b1p = b1 - eta * db1

        g0k = con_lower(w0, b0, c0)
        lam0 = relu0(lam1 + gam2 * g0k)

        F0k = Loss_upper(w0, b0, c0, y_val, x_val)
        f0k = Loss_lower(w0, b0)
        f1k = Loss_lower(w1p, b1p)
        g0k = lam0 @ con_lower(w0, b0, c0)

        dc0 = 1/ck * grad(F0k, c0, retain_graph=True)[0]+grad(g0k, c0, retain_graph=True)[0]-grad((lam1 @ con_lower(w1, b1, c0)), c0, retain_graph=True)[0]
        c0p = c0 - alpha * dc0

        dw0 = 1/ck * grad(F0k, w0, retain_graph=True)[0] + grad(f0k, w0, retain_graph=True)[0] + grad(g0k, w0, retain_graph=True)[0] - 1 / gam1 * (w0 - w1)
        db0 = 1/ck * grad(F0k, b0, retain_graph=True)[0] + grad(g0k, b0, retain_graph=True)[0] - 1 / gam1 * (b0 - b1)

        dlam1 = -(lam1 - lam0) / gam2 - con_lower(w1, b1, c0)
        with torch.no_grad():
            w0 = w0 - alpha * dw0
            b0 = b0 - alpha * db0
            lam1 = torch.minimum(torch.maximum(lam1 - beta * dlam1, torch.zeros(nTr)), R*torch.ones(nTr))

        w1, b1, c0 = w1p.detach().clone(), b1p.detach().clone(), c0p.detach().clone()




  
    


        #################
        if epoch%20==0 and verbose:
            print(f"Epoch [{epoch}/{epochs}]:",
              "val acc: {:.2f}".format(val_acc),
              "val loss: {:.2f}".format(val_loss),
              "test acc: {:.2f}".format(test_acc),
              "test loss: {:.2f}".format(test_loss))

        val_loss_list.append(val_loss) # length 150
        test_loss_list.append(test_loss) # length 118
        val_acc_list.append(val_acc)
        test_acc_list.append(test_acc)
        time_computation.append(time.time()-algorithm_start_time)

        if torch.linalg.norm(dc0) < early_stopping_th:
            break

    return metrics, variables


if __name__ == "__main__":
    ############ Load data code ###########

    data = load_diabetes()

    n_train = 500
    n_val = 150

    metrics = []
    variables = []

    hparams = {
        'p': 0.3,
        'gam1': 10,
        'gam2': 0.01,
        'eta': 0.01,
        'alpha': 0.001,
        'beta': 0.001,
        'ck0': 10,
        'R': 10,
    }

    epochs = 5000
    plot_results = True

    for seed in range(10):

        x_train, y_train, x_val, y_val, x_test, y_test = train_val_test_split(data, seed, n_train, n_val)

        metrics_seed, variables_seed = gaffa(x_train, y_train, x_val, y_val, x_test, y_test, hparams, epochs)
        metrics.append(metrics_seed)
        variables_seed.append(variables_seed)

    train_acc = np.array([[x['train_acc'] for x in metric] for metric in metrics])
    val_acc = np.array([[x['val_acc'] for x in metric] for metric in metrics])
    test_acc = np.array([[x['test_acc'] for x in metric] for metric in metrics])

    val_loss = np.array([[x['val_loss'] for x in metric] for metric in metrics])
    test_loss = np.array([[x['test_loss'] for x in metric] for metric in metrics])

    time_computation = np.array([[x['time_computation'] for x in metric] for metric in metrics])

    if plot_results:
        val_loss_mean=np.mean(val_loss,axis=0)
        val_loss_sd=np.std(val_loss,axis=0)/2.0
        test_loss_mean=np.mean(test_loss,axis=0)
        test_loss_sd=np.std(test_loss,axis=0)/2.0

        val_acc_mean=np.mean(val_acc,axis=0)
        val_acc_sd=np.std(val_acc,axis=0)/2.0
        test_acc_mean=np.mean(test_acc,axis=0)
        test_acc_sd=np.std(test_acc,axis=0)/2.0

        axis = np.mean(time_computation,axis=0)

        plt.rcParams.update({'font.size': 18})
        plt.rcParams['font.sans-serif']=['Arial']#如果要显示中文字体，则在此处设为：SimHei
        plt.rcParams['axes.unicode_minus']=False #显示负号
        axis=time_computation.mean(0)
        plt.figure(figsize=(8,6))
        #plt.grid(linestyle = "--") #设置背景网格线为虚线
        ax = plt.gca()
        plt.plot(axis,val_loss_mean,'-',label="Training loss")
        ax.fill_between(axis,val_loss_mean-val_loss_sd,val_loss_mean+val_loss_sd,alpha=0.2)
        plt.plot(axis,test_loss_mean,'--',label="Test loss")
        ax.fill_between(axis,test_loss_mean-test_loss_sd,test_loss_mean+test_loss_sd,alpha=0.2)
        #plt.xticks(np.arange(0,iterations,40))
        plt.title('Kernelized SVM')
        plt.xlabel('Running time /s')
        #plt.legend(loc=4)
        plt.ylabel("Loss")
        #plt.xlim(-0.5,3.5)
        #plt.ylim(0.5,1.0)
        plt.legend(loc=0, numpoints=1)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        #plt.setp(ltext, fontsize=18,fontweight='bold') #设置图例字体的大小和粗细
        plt.savefig('gaffa_loss.png') 
        #plt.show()

        plt.figure(figsize=(8,6))
        ax = plt.gca()
        plt.plot(axis,val_acc_mean,'-',label="Training accuracy")
        ax.fill_between(axis,val_acc_mean-val_acc_sd,val_acc_mean+val_acc_sd,alpha=0.2)
        plt.plot(axis,test_acc_mean,'--',label="Test accuracy")
        ax.fill_between(axis,test_acc_mean-test_acc_sd,test_acc_mean+test_acc_sd,alpha=0.2) 
        #plt.xticks(np.arange(0,iterations,40))
        plt.title('Kernelized SVM')
        plt.xlabel('Running time /s')
        plt.ylabel("Accuracy")
        # plt.ylim(0.64,0.8)
        #plt.legend(loc=4)
        plt.legend(loc=0, numpoints=1)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        #plt.setp(ltext, fontsize=18,fontweight='bold') #设置图例字体的大小和粗细
        plt.savefig('gaffa_acc.png') 
        plt.show()
