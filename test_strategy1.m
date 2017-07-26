% test_strategy1: 
% Description: test strategy 1: use K burst generated from K random initializations
% Tuning parameters:
%    n: number of variable (should be at least around 50 so that the lower bound on NumIC is smaller than the number of terms in the basis)
%    F: constant parameter of Lorenz 96 -- if F>1, the system is chaos; if F<1, the system eventually converges to a point.
%    NumIC: number of initializations

% Other parameters:
%    dt: timestep 
%    SizeofBurst: m, number of measurements for each initialization
%    N: number of terms in the monmial/legendre basis, N = (n^2+3n+2)/2
%    lowerBoundNumIC: the lower bound of NumIC provided by our theory given that the universal constant is 1, and epsilon = 0.5
%    upperBoundNumIC: the upper bound of NumIC 
%    opts: parameters for the optimization algorithm spgl1
%    optPolynomial: 'monomial' or 'legendre' to build the dictionary matrix
%    optEquation: Equation to test the reovery, should be from 1 to n

% Output: recovered coefficients of Equation optEquation

% Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
% Version 1, July 2017 
% Reference: arxiv link
%            Download the optimization package spgl1 from http://www.cs.ubc.ca/~mpf/spgl1/

close all; clear all; clc

%% ODE parameters
% Tuning parameters
n = 50; % number of variables
F = 8.0; % constant of Lorenz 96
optEquation = 10; % Equation to test
NumIC = 100;

% Other parameters

N = (n+1)*(n+2)/2; % number of columns of the dictionary matrix 
dt = 0.001; % time step
SizeOfBurst = 5; % size of each burst

lowerBoundNumIC = round(5*log(N) * log(1/0.5)); % s*log(N)*log(1/varepsilon)
upperBoundNumIC = round(N/SizeOfBurst);

display(['The number of initializations NumIC should be at least ', num2str(lowerBoundNumIC),'c',' and be smaller than ',num2str(upperBoundNumIC)]);

% spgl1 parameters
opts = [];
opts.verbosity = 0; 
opts.iterations = 1000;

 % Option parameters
optPolynomial = 'legendre'; % 'legendre' or 'monomial'

% True Coefficients
c_true_mat = Lorenz96_true_coefficients(n,F);

%% Data generated from K bursts starting from K random initializations 
Xint = 2*rand(n,NumIC)-1; % initialization is a uniform random variable on [-1,1]
[Xfull, Vapproximate,Vexact] =  Lorenz96_XV(F,Xint,dt,SizeOfBurst); 

% Built dictionary
D = dictionary96(Xfull,optPolynomial);


%% Basis Pursuit Denoising Problem
sigma = 2.*norm(Vapproximate(:,optEquation)-Vexact(:,optEquation),2);
soln = basisPursuit_Lorenz96(Vapproximate,D,optEquation,optPolynomial,opts,sigma);

%% Print out
display(['The nonzero terms of the true coefficients from Equation ',num2str(optEquation)])
sparse(c_true_mat(:,optEquation))
display(['The nonzero terms of the recovered coefficients from Equation ',num2str(optEquation)])
sparse(soln)