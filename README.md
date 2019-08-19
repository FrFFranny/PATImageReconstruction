# PAT Image Reconstruction
Matlab codes for PAT image reconstruction from subsampled data based on a novel regularisation term (Hessian Schatten-norm of the image filtered by Gaussian function), using k-Wave Matlab toolbox, FISTA and ADMM algorithms.

### Auther: Qing Fang
August 20, 2019

# Build Instruction
This project is developed based on the k-Wave matlab toolbox, hence the k-Wave toolbox should be installed first.   
[k-wave toolbox]http://www.k-wave.org   
  
Then the main_function.m script can be used for image reconstruction. By executing the main function, the PAT image cen be recovered for a specified subsampled data based on the FISTA algorithm with single Gaussian filter or with two Gaussian filters combined via ADMM algorithm.
