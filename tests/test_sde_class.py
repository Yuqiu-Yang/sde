import pytest
import os
from copy import deepcopy
import numpy as np 
from sde.sde_class import sde_class

@pytest.mark.parametrize("T, N, M, expected", [
    (1, 1000, 2, (2, 1000))
])
def test_init(T, N, M, expected):
    sde = sde_class(T=T, N=N, M=M)
    assert (sde.W.shape == expected) and (sde.dW.shape == expected)
    
@pytest.mark.parametrize("T, N, M, W, expected", [
    (1,100,5, None, (5,100)),
    (1,100,5, np.random.randn(5,100), (5,100)),
])
def test_transform(T, N, M, W, expected):
    sde = sde_class(T=T, N=N, M=M)
    def u(t, w):
        return np.exp(t + 0.5 * w)
    transformed_W = sde.transform_W(fun=u,
                                    W=W)
    assert transformed_W.shape == expected
    
