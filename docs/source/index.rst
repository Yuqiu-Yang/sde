.. sde documentation master file, created by
   sphinx-quickstart on Sat Apr 29 18:29:22 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ../../assets/logo.png
   :width: 600

**sde**: a simple tool for numerical stochastic differential equations
===========================================================================

**sde**provides basic tools for simulating brownian motions which is the basic 
ingredients to lots of SDE models, performing different types of stochastic integrations, 
Eular-Maruyama methods and so much more (to come...). 

The origin of this project is this wonderful introductory 
`paper <https://epubs.siam.org/doi/pdf/10.1137/S0036144500378302>`_:
by Desmond J. Higham.

In this document, we will demonstrate how to:

* Simulate brownian motion paths 
* Transform simulated bm paths 
* Perform stochastic integrals 
* Carry out numerical sde methods 

However, before we get started, we need to install the 
package first which is freely available on PyPI. 

Dependencies
---------------

.. list-table:: Dependencies
  :widths: 50 50
  :align: center
  :header-rows: 1

  * - Package 
    - Version 
  * - python
    - ``>=3.9``  
  * - numpy
    - ``==1.24.2``
  * - tqdm 
    - ``==4.64.1``
  * - matplotlib
    - ``==3.7.1``


Installation
--------------------------------

.. note:: 
   We highly recommend creating a virtual environment before 
   proceeding to installing the package. For how to manage virtual
   environment via ``conda``, check out 
   `their tutorial <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#>`_.

.. code:: bash
   
   pip install sde

To quickly test if it has been installed:

.. code:: bash

  python -m sde --version 


Tutorial
===================
.. toctree::
   :maxdepth: 1
   :caption: Brownian Motions

   bm

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api_reference/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
