Extracting Sparse High-Dimensional Dynamics from Limited Data
Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
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
7. test_all_NumIC.m
	phase transition 1D plot where number of initializations varies from the 		sparse level to N/(SizeOfBurst). For each NumIC, run 100 times and compute 		the probability of recovering exactly 4 monomial terms of the governing 		equation.
note: The recovered coefficients are close to the ground truth.
      To test with a specific number of initializations, assign NumICmin = NumICmax 
      The lower bound of number of initializations is suggested from the theorem in our paper.

8. test_Lorenz96_comparison.m
      Plot the basis pursuit solution, the least-square solution and the solution from the sequential thresholding algorithm.  

B. Fisher's Equation

 
