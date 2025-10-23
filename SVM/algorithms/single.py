import numpy as np
import cvxpy as cp
import time

import matplotlib.pyplot as plt

import torch
from torch.nn import functional as F

import sys
sys.path.append('..')

device = torch.device("cpu")

from utils import load_diabetes, train_val_test_split



def single(x_train, y_train, x_val, y_val, x_test, y_test, hparams, epochs, early_stopping_th = False,verbose=True):
    feature=x_train.shape[1] # = 8
    nTr=y_train.shape[0]
    mu=0.01

    # Dataset to tensor
    y_train = torch.tensor(y_train).float()
    x_train = torch.tensor(x_train).float()
    y_val = torch.tensor(y_val).float()
    x_val = torch.tensor(x_val).float()
    y_test = torch.tensor(y_test).float()
    x_test = torch.tensor(x_test).float()
    
    Loss_lower = lambda w, b : .5 * (torch.sum(torch.square(w)) + mu* b**2)
    con_lower = lambda w, b, c: 1 - c - y_train * (F.linear(x_train, w, b).squeeze())


    ###### Parameters
    eta_x = hparams['eta_x']
    eta_y = hparams['eta_x']
    eta_z = hparams['eta_x']
    eta_u = hparams['eta_x']
    eta_v = hparams['eta_x']
    delta = hparams['delta']
    rho_1 = hparams['rho_1']
    rho_2 = hparams['rho_2']

    # Initialization of upper and lower level variables
    C_tensor_val= torch.Tensor(x_train.shape[0]).uniform_(1.,5.)

    w_z = torch.zeros(1,feature)
    b_z = torch.tensor(1.)
    u = torch.zeros(nTr)
    alpha = torch.zeros(nTr)

    w_y = w_z.clone()
    b_y = b_z.clone()
    v = torch.zeros(nTr)
    beta = torch.zeros(nTr)
    
    # For storage
    val_loss_list=[]
    test_loss_list=[]
    val_acc_list=[]
    test_acc_list=[]
    time_computation=[]
    algorithm_start_time=time.time()

    metrics = []
    variables = []


    # Set the Initial Value
    svm_w = cp.Variable(feature)
    svm_b = cp.Variable()
    loss = 1/2 *(cp.sum_squares(svm_w) + mu * cp.square(svm_b))
    constraints = [1 - C_tensor_val.cpu().numpy() - cp.multiply(y_train.cpu().numpy(), (x_train.cpu().numpy() @ svm_w + svm_b)) <= 0]
    prob = cp.Problem(cp.Minimize(loss), constraints)

    algorithm_start_time = time.time()
    # prob.solve(solver='ECOS')
    prob.solve(solver='ECOS', abstol=2e-3,reltol=2e-3,max_iters=1000000000, warm_start=True)  

    w_z.data = torch.tensor(svm_w.value, device = device, dtype=torch.float32)
    w_z=w_z.reshape(1,feature)
    b_z.data = torch.tensor(svm_b.value, device = device, dtype=torch.float32)
    w_y.data = torch.clone(w_z)
    b_y.data = torch.tensor(svm_b.value, device = device, dtype=torch.float32)

    alpha.data = torch.clone(con_lower(w_y, b_y, C_tensor_val))
    beta.data = torch.clone(alpha)



    u.data = torch.tensor(constraints[0].dual_value, device = device, dtype=torch.float32)
    v.data = torch.tensor(constraints[0].dual_value, device = device, dtype=torch.float32)



    C_tensor_val.requires_grad_(False)
    w_z.requires_grad_(False)
    b_z.requires_grad_(False)
    w_y.requires_grad_(False)
    b_y.requires_grad_(False)
    u.requires_grad_(False)
    v.requires_grad_(False)
    alpha.requires_grad_(False)
    beta.requires_grad_(False)


    for epoch in range(epochs):

        variables.append({
            'C': C_tensor_val,
            'w': w_z,
            'b': b_z,
            'w_F': w_y,
            'b_F': b_y
        })

        x = torch.reshape(y_val, (y_val.shape[0],1)) 
        x = x* F.linear(x_val, w_y, b_y) # / torch.linalg.norm(w_tensor)
        #x = x* (x_val @ w_y+b_y)

        x1 = torch.reshape(y_test, (y_test.shape[0],1)) 
        x1 = x1 * F.linear(x_test, w_y, b_y) # / torch.linalg.norm(w_tensor)
        #x1= x1* (x_test @ w_y +b_y)
        # test_loss_upper= torch.sum(torch.sigmoid(x1))
        test_loss_upper= torch.sum(torch.exp(1-x1))

        val_loss_F = (torch.sum(torch.exp(1-x))).detach().numpy()/y_val.shape[0]
        test_loss_F = test_loss_upper.detach().numpy()/y_test.shape[0]

        x = torch.reshape(y_val, (y_val.shape[0],1))
        x = x* F.linear(x_val, w_z, b_z) # / torch.linalg.norm(w_z)

        #x= x* (x_val @ w_z +b_z)


        x1 = torch.reshape(y_test, (y_test.shape[0],1)) 
        x1 = x1 * F.linear(x_test, w_z, b_z) # / torch.linalg.norm(w_z)

        #x1= x1* (x_test @ w_z +b_z)

        test_loss_upper= torch.sum(torch.exp(1-x1))

        val_loss = (torch.sum(torch.exp(1-x))).detach().numpy()/y_val.shape[0]
        test_loss = test_loss_upper.detach().numpy()/y_test.shape[0]

        loss_upper = val_loss

        ###### Accuracy
        q = y_train * (w_z @ x_train.T + b_z)
        train_acc = (q>0).sum() / len(y_train)

        q = y_val * (w_z @ x_val.T + b_z)
        val_acc = (q>0).sum() / len(y_val)

        q = y_test * (w_z @ x_test.T + b_z)
        test_acc = (q>0).sum() / len(y_test)

        q = y_train * (w_y @ x_train.T + b_y)
        train_acc_F = (q>0).sum() / len(y_train)

        q = y_val * (w_y @ x_val.T + b_y)
        val_acc_F = (q>0).sum() / len(y_val)

        q = y_test * (w_y @ x_test.T + b_y)
        test_acc_F = (q>0).sum() / len(y_test)

        metrics.append({
            'train_acc_z': train_acc,
            'train_acc_y': train_acc_F,
            'val_loss_z': val_loss,
            'val_loss_y': val_loss_F,
            'val_acc_z': val_acc,
            'val_acc_y': val_acc_F,
            'test_loss_z': test_loss,
            'test_loss_y': test_loss_F,
            'test_acc_z': test_acc,
            'test_acc_y': test_acc_F,
            'time_computation': time.time()-algorithm_start_time
        })






        # Update
        C_tensor_val.requires_grad_(True)
        w_z.requires_grad_(True)
        b_z.requires_grad_(True)
        w_y.requires_grad_(True)
        b_y.requires_grad_(True)
        u.requires_grad_(True)
        v.requires_grad_(True)
        alpha.requires_grad_(True)
        beta.requires_grad_(True)

        C_tensor_val.grad = None # Reset gradients
        w_z.grad=None
        b_z.grad=None
        w_y.grad=None
        b_y.grad=None
        u.grad=None
        v.grad=None
        alpha.grad=None
        beta.grad=None

        x = torch.reshape(y_val, (y_val.shape[0],1)) 
        x = x* F.linear(x_val, w_y, b_y)
        #x = x* (x_val @ w_y +b_y)
        loss_upper= torch.sum(torch.exp(1-x)) + 0.5 * torch.linalg.norm(C_tensor_val)**2

        h_alpha=con_lower(w_y, b_y, C_tensor_val)-alpha
        h_beta=con_lower(w_z, b_z, C_tensor_val)-beta



        
        K= delta* loss_upper+Loss_lower(w_y, b_y)-Loss_lower(w_z,b_z) + torch.dot(u, h_alpha )-torch.dot(v, h_beta )+0.5 * rho_1 * torch.dot(h_alpha, h_alpha)-0.5 * rho_2 * torch.dot(h_beta, h_beta )
        
        # Upper level iteration


        K.backward()
        
        with torch.no_grad():
            u.data.add_(u.grad.data, alpha=eta_u)
            v.data.add_(v.grad.data, alpha=-eta_v)
            C_tensor_val.data.add_(C_tensor_val.grad.data, alpha=-eta_x)
            w_y.data.add_(w_y.grad.data, alpha=-eta_y)
            b_y.data.add_(b_y.grad.data, alpha=-eta_y)
            w_z.data.add_(w_z.grad.data, alpha=eta_z)
            b_z.data.add_(b_z.grad.data, alpha=eta_z)

            alpha.data.add_(alpha.grad.data, alpha=-eta_y)
            beta.data.add_(beta.grad.data, alpha=eta_z)
            alpha = torch.clamp(alpha, max=0.0)
            beta = torch.clamp(beta, max=0.0)

            #C_tensor_val.data = torch.maximum(C_tensor_val.data, torch.tensor(1e-4))


        C_tensor_val.requires_grad_(False)
        w_z.requires_grad_(False)
        b_z.requires_grad_(False)
        w_y.requires_grad_(False)
        b_y.requires_grad_(False)
        u.requires_grad_(False)
        v.requires_grad_(False)
        alpha.requires_grad_(False)
        beta.requires_grad_(False)

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

        if torch.linalg.norm(C_tensor_val.grad.detach()) < early_stopping_th:
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
        'delta': 0.1,
        'eta_x': 0.001,
        'eta_y': 0.001,
        'eta_z': 0.001,
        'eta_u': 0.001,
        'eta_v': 0.001,
        'rho_1': 0.1,
        'rho_2': 0.1,
    }

    epochs = 200
    plot_results = True

    for seed in range(10):

        x_train, y_train, x_val, y_val, x_test, y_test = train_val_test_split(data, seed, n_train, n_val)

        metrics_seed, variables_seed = single(x_train, y_train, x_val, y_val, x_test, y_test, hparams, epochs)
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
        plt.savefig('ho_svm_kernel_1.png') 
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
        plt.savefig('ho_svm_kernel_2.png') 
        plt.show()
