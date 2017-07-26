function xdot = lorenz96(t,x,n,F)
% Description: RHS of Lorenz 96: 
%     dx_{k}/dt = -x_{k-2} * x_{k-1} + x_{k-1} * x_{k+1} - x_k + F,
%     k=1,...,n

% Input: t: time (default)
%        x: variable in R^n
%        F: the constant parameter of Lorenz 96
% Output: the right hand side of Lorenz 96 (which is used to solve the ODE system
% and to find the exact velocity if needed)

% Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
% Version 1, July 2017 
% Reference: arxiv link

xdot = zeros(n,1);
xdot(1) = - x(n-1) * x(n) + x(n) * x(2) - x(1) + F;
xdot(2) = - x(n) * x(1) + x(1) * x(3) - x(2) + F;
xdot(n) = - x(n-2) * x(n-1) + x(n-1) * x(1) - x(n) + F;
for ind = 3:n-1
    xdot(ind) = - x(ind-2) * x(ind-1) + x(ind-1) * x(ind+1) - x(ind) + F;
end

% x3 = - x1*x2 + x2*x4 - x3 + F -- index 1, 4, N+3