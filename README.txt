Extracting Sparse High-Dimensional Dynamics from Limited Data
Copyright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
Version 1, July 2017 
Reference: arxiv link

A. Lorenz96 dx/dt = F(x,t)
 
dx_{k}/dt = -x_{k-2} * x_{k-1} + x_{k-1} * x_{k+1} - x_k + F,  k=1,...,N

1. dictionary96.m   
	Include dictionary up to degree 2 built from monomials and from Legendre polynomials
2. lorenz96.m       
	Right-hand-side of the ODE Lorenz96
3. time_derivative.m
        Approximate velocity from data using 1st/2nd approximation
4. Lorenz96_XV.m
	Construct all-in-one including the data matrix, approximation
	velocity, and exact velocity given F, initialization, timestep, and
	length of the evolution.
5. Lorenz96_true_coefficients.m
  	Compute the true coefficients of the Lorenz 96 and their indices  
6. basisPursuit_Lorenz96.m
	Find the coefficient of component optEquation in Lorenz 96

7. test_strategy1.m
	Use strategy 1, where data is collected from K bursts starting from K random uniform initializations (see more in our paper), compute the recovered coefficients of any component of the Lorenz 96.
 
