# PCA after Diagonal Reduction

* This paper has been superseded by the follow-up work: E. Dobriban, W. Leeb, A. Singer: **Optimal prediction in the linearly transformed spiked model**. available at [arxiv](https://arxiv.org/abs/1709.03393), and to appear in The Annals of Statistics

Code to reproduce the computational results from the paper:

* E. Dobriban, W. Leeb, A. Singer: **PCA from noisy, linearly reduced data: the diagonal case**. available at [arxiv](https://arxiv.org/abs/1611.10333):

This directory contains work by W. Leeb and E. Dobriban.

# Acknowledgements
* Uses software from E. Dobriban's EigenEdge package, also available on [GitHub](github.com/dobriban/EigenEdge). 
The functions from EigenEdge used here are:
 ```general_spiked_forward.m ```.
 ```compute_esd_ode.m ```.
 ```compute_esd_fp.m ```.
 ```evaluate_inverse_ST.m ```.
 ```standard_spiked_forward.m ```.

* Includes for completeness the "OptSpace" method, from the paper 
Keshavan, Montanari, Oh: OptSpace: A Matrix Completion Algorithm. http://arxiv.org/abs/0901.3150.
This is used in the experiments in the  ```expt_shri2/ ``` folder.
