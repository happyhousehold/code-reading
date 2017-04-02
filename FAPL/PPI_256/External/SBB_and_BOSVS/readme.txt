    SBB algorithm by Xiaojing Ye

    WHAT IS THIS
This code solves isotropic image reconstruction with total variation (TV)
regularization model:
    minimize_u (wTV)*TV(u) + 0.5*|Au-f|^2   

where A is sensing matrix (using parallel MRI in this code for
demonstration in this code), f is input data, u is the image to be
reconstructed, wTV>0 is the weight of TV term.

    HOW TO USE
For demonstration, simply call 'demo_SBB' in MATLAB command window.
For general TV based image reconstruction, modify 'A_oper' and 'At_oper'
functions in 'utilities' folder. Detailed instructions are provided as
comments in the code.

    REFERENCE
This software implements the algorithm developed in this paper:
Computational Acceleration for MR Image Reconstruction in Partially Parallel Imaging
Xiaojing Ye, Yunmei Chen, and Feng Huang
IEEE Transactions on Medical Imaging (TMI), 30(5), pp. 1055-1063, 2011.
Modification with backtracking for convergence will be online soon.

    DISCLAIMER
This code is for academic (non-commercial) use only.
This code comes with NO warranty.