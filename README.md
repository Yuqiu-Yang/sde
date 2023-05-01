# Stochastic Differential Equations :cyclone:

![Logo](/assets/logo.png)

![forthebadge](/assets/numerical-sde.svg)
![forthebadge](/assets/fokker-planck.svg)

| Package | Documentation | Code Coverage |
| --- | --- | --- |
| sde |  [![Documentation Status](https://readthedocs.org/projects/sde/badge/?version=latest)](https://sde.readthedocs.io/en/latest/?badge=latest) | [![codecov](https://codecov.io/gh/Yuqiu-Yang/sde/branch/main/graph/badge.svg?token=KW3cp0XJky)](https://codecov.io/gh/Yuqiu-Yang/sde) |

So you say you are interested in stochastic differential equations? And you also subscribe to the idea of learning by example? Well, you are in luck! :smiley: 
<b>sde</b> provides basic tools for simulating brownian motions which is the basic ingredients to lots of SDE models, performing different types of stochastic integrations, Eular-Maruyama methods and so much more (to come... :stuck_out_tongue_winking_eye:). 

We have also built a detailed [online documentation :page_with_curl:](https://sde.readthedocs.io/en/latest/) where we guide you step-by-step on how to use our package.

The origin of this project is this wonderful introductory [paper :page_facing_up:](https://epubs.siam.org/doi/pdf/10.1137/S0036144500378302) by Desmond J. Higham.

## Dependencies :anchor:
- numpy==1.22.4
- tqdm==4.64.1
- matplotlib==3.7.1
- pandas==2.0.0
- seaborn==0.12.2

## Enviroment Setup :arrow_forward:
We highly recommend creating a virtual environment before proceeding to installing the package. For how to manage virtual
environment via conda, check out 
[their tutorial](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#).

```shell
pip install -r requirements.txt
```

## Installation :arrow_double_down:
```shell
pip install sdeIU
```

## Tutorials
We also included all the [jupyter notebooks :notebook_with_decorative_cover:](https://github.com/Yuqiu-Yang/sde/tree/main/notebooks). If you want, you can 
regenerate (of course not exactly due to the random nature of the problems) all the examples we went through in [our documentation :page_with_curl:](https://sde.readthedocs.io/en/latest/) by running these notebooks.
- [Brownian Motions :bread:](https://github.com/Yuqiu-Yang/sde/blob/main/notebooks/bm.ipynb)
- [Stochastic Integrals :egg:](https://github.com/Yuqiu-Yang/sde/blob/main/notebooks/integrals.ipynb)
- [Stochastic Differential Equations :cake:](https://github.com/Yuqiu-Yang/sde/blob/main/notebooks/sde.ipynb)
- [Convergence and Stability :custard:](https://github.com/Yuqiu-Yang/sde/blob/main/notebooks/convergence.ipynb)

