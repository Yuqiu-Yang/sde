# Stochastic Differential Equations

![Logo](/assets/logo.png)

![forthebadge](/assets/numerical-sde.svg)
![forthebadge](/assets/fokker-planck.svg)

| Package | Documentation | Code Coverage |
| --- | --- | --- |
| sde |  [![Documentation Status](https://readthedocs.org/projects/sde/badge/?version=latest)](https://sde.readthedocs.io/en/latest/?badge=latest) | [![codecov](https://codecov.io/gh/Yuqiu-Yang/sde/branch/main/graph/badge.svg?token=KW3cp0XJky)](https://codecov.io/gh/Yuqiu-Yang/sde) |

So you say you are interested in stochastic differential equations? And you also subscribe to the idea of learning by example? Well, you are in luck!
<b>sde</b> provides basic tools for simulating brownian motions which is the basic ingredients to lots of SDE models, performing different types of stochastic integrations, Eular-Maruyama methods and so much more (to come...). 

We have also built a detailed [online documentation](https://sde.readthedocs.io/en/latest/) where we guide you step-by-step on how to use our package.

The origin of this project is this wonderful introductory [paper](https://epubs.siam.org/doi/pdf/10.1137/S0036144500378302) by Desmond J. Higham.

## Dependencies
- numpy==1.22.4
- tqdm==4.64.1
- matplotlib==3.7.1
- pandas==2.0.0
- seaborn==0.12.2

## Enviroment Setup
We highly recommend creating a virtual environment before proceeding to installing the package. For how to manage virtual
environment via conda, check out 
[their tutorial](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#).

```shell
pip install -r requirements.txt
```

## Installation 
```shell
pip install sde
```
