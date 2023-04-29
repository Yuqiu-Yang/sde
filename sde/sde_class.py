import numpy as np 
from matplotlib import pyplot as plt 

from tqdm import tqdm 
from typing import Callable, Optional

class sde_class:
    def __init__(self,
                 T: float, 
                 N: int,
                 M: int) -> None:
        self.T = T
        self.N = N
        self.M = M
        self.dt = T/N
        self.dW = np.sqrt(self.dt) * np.random.randn(M, N-1)
        self.dW = np.pad(self.dW, [(0,0),(1,0)], mode="constant", 
                         constant_values=0)
        self.W = np.cumsum(self.dW, axis=1)
        self.time = np.linspace(start = 0,
                                stop = T,
                                num=N)

    def transform_W(self, 
                    fun: Callable,
                    W: Optional[np.ndarray]=None) -> np.ndarray:
        if W is None:
            t_W = [fun(self.time, self.W[path, :]) for path in range(self.W.shape[0])]
        else:
            t_W = [fun(self.time, W[path, :]) for path in range(W.shape[0])]
        return np.array(t_W)
    
    def integrate(self, 
                  fun: Callable,
                  integral_type: str="Ito") -> np.ndarray:
        if integral_type == "Stratonovich":
            # Compute mid points
            W_array = np.zeros((self.W.shape[0], np.ceil(self.W.shape[1]/2)))
            for path in range(self.W.shape[0]):
                ma = np.convolve(self.W[path,:], np.ones(2), "valid")/2
                ma += np.random.randn(len(ma))*np.sqrt(self.dt)*0.5
                W_array[path,:] = ma
            t_W = self.transform_W(fun=fun,
                                   W=W_array) 
        elif integral_type == "Ito":
            t_W = self.transform_W(fun=fun)
        else:
            return NotImplementedError
        
        integral_result = np.zeros(self.W.shape[0])
        if integral_type == "Ito":
            for k in range(t_W.shape[0]):
                integral_result[k] = np.sum(self.dW[k, 1:] * t_W[k, :-1]) 
        elif integral_type == "Stratonovich":
            for k in range(t_W.shape[0]):
                integral_result[k] = np.sum(self.dW[k, 1:] * t_W[k,:])  
        else:
            raise NotImplementedError
        return integral_result
        
    def eular_maruyama(self,
                       mu_fun: Callable,
                       sigma_fun: Callable,
                       x0: float=0,
                       R: int=1) -> dict:
        dW = np.diff(self.W[:,::R], axis=1)
        dt = self.dt * R
        X = np.zeros((dW.shape[0], dW.shape[1]+1))
        X[:,0] = x0
        for j in tqdm(range(1, X.shape[1])):
            X[:,j] = X[:,j-1]
            X[:,j] += dt * mu_fun(X[:, j-1])
            X[:,j] += dW[:, j-1] * sigma_fun(X[:, j-1])
        return X
    
    def milstein(self,
                 mu_fun: Callable,
                 sigma_fun: Callable,
                 d_sigma_fun: Callable,
                 x0: float=0,
                 R: int=1) -> dict:
        dW = np.diff(self.W[:,::R], axis=1)
        dt = self.dt * R
        X = np.zeros((dW.shape[0], dW.shape[1]+1))
        X[:,0] = x0
        for j in tqdm(range(1, X.shape[1])):
            X[:,j] = X[:,j-1]
            X[:,j] += dt * mu_fun(X[:, j-1])
            X[:,j] += dW[:, j-1] * sigma_fun(X[:, j-1])
            X[:,j] += 0.5 * (dW[:,j-1]**2 - dt) * sigma_fun(X[:,j-1]) * d_sigma_fun(X[:,j-1])
        return X
        
    
    
    