# Numerical manipulation 
import numpy as np 

# Plotting 
from matplotlib import pyplot as plt 

# User entertainment
from tqdm import tqdm 

# Typing 
from typing import Callable, Optional


class sde_class:
    def __init__(self,
                 T: float, 
                 N: int,
                 M: int) -> None:
        """Initialize sde object 

        Parameters
        ----------
        T : float
            The end time point 
        N : int
            Number of time steps 
        M : int
            Number of BM paths 
        """
        self.T = T
        self.N = N
        self.M = M
        # Compute time steps 
        self.dt = T/N
        # Generate increments 
        # We pad 0 at the front so that 
        # the length is equal to that of 
        # the resulting BMs 
        self.dW = np.sqrt(self.dt) * np.random.randn(M, N-1)
        self.dW = np.pad(self.dW, [(0,0),(1,0)], mode="constant", 
                         constant_values=0)
        # Compute the resutling BMs
        self.W = np.cumsum(self.dW, axis=1)
        # Generate time for plotting purpose 
        self.time = np.linspace(start=0,
                                stop=T,
                                num=N)

    def transform_W(self, 
                    fun: Callable,
                    W: Optional[np.ndarray]=None) -> np.ndarray:
        """Apply transformation to BMs 

        Parameters
        ----------
        fun : Callable
            A funtion of time and brownian motions 
        W : Optional[np.ndarray], optional
            If None, we will use the stored BMs 
            Ohterwise, we will apply the transformation to 
            the user supplied arrays 

        Returns
        -------
        np.ndarray
            An array of transformed BMs 
        """
        if W is None:
            t_W = [fun(self.time, self.W[path, :]) for path in range(self.W.shape[0])]
        else:
            t_W = [fun(self.time, W[path, :]) for path in range(W.shape[0])]
        return np.array(t_W)
    
    def integrate(self, 
                  fun: Callable,
                  integral_type: str="Ito") -> np.ndarray:
        """Stochastic integration 

        Parameters
        ----------
        fun : Callable
            A function of time and BMs 
        integral_type : str, optional
            The type of stochastic integration to perform
            Can be "Ito" or "Stratonovich", by default "Ito".

        Returns
        -------
        np.ndarray
            An array of computed integral results 

        Raises
        ------
        NotImplementedError
            If type not of Ito or Stratonovich, the 
            algorithm is not implemented 
        """
        if integral_type == "Stratonovich":
            # Compute mid points
            W_array = np.zeros((self.W.shape[0], self.W.shape[1]-1))
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
        
    def euler_maruyama(self,
                       mu_fun: Callable,
                       sigma_fun: Callable,
                       x0: float=0,
                       R: int=1) -> np.ndarray:
        """Euler Maruyama method 

        Parameters
        ----------
        mu_fun : Callable
            A drift function of x
        sigma_fun : Callable
            A diffusion function of x 
        x0 : float, optional
            The initial value, by default 0
        R : int, optional
            The time steps, by default 1

        Returns
        -------
        np.ndarray
            An array of simulated sde values 
        """
        # Compute increments and time steps 
        dW = np.diff(self.W[:,::R], axis=1)
        dt = self.dt * R
        # Sequentially compute sde
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
                 R: int=1) -> np.ndarray:
        """Milstein's higher order method 

        Parameters
        ----------
        mu_fun : Callable
            A drift function of x
        sigma_fun : Callable
            A diffusion function of x 
        d_sigma_fun : Callable
            The derivative of the diffusion function of x 
        x0 : float, optional
            The initial value, by default 0
        R : int, optional
            The time steps, by default 1

        Returns
        -------
        np.ndarray
            An array of simulated sde values 
        """
        # Compute increments and time steps 
        dW = np.diff(self.W[:,::R], axis=1)
        dt = self.dt * R
        # Sequentially compute sde
        X = np.zeros((dW.shape[0], dW.shape[1]+1))
        X[:,0] = x0
        for j in tqdm(range(1, X.shape[1])):
            X[:,j] = X[:,j-1]
            X[:,j] += dt * mu_fun(X[:, j-1])
            X[:,j] += dW[:, j-1] * sigma_fun(X[:, j-1])
            X[:,j] += 0.5 * (dW[:,j-1]**2 - dt) * sigma_fun(X[:,j-1]) * d_sigma_fun(X[:,j-1])
        return X
        
    
    
    