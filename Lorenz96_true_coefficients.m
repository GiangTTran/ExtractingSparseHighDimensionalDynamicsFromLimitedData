function c_true_mat = Lorenz96_true_coefficients(n,F)
% Description: compute the true coefficients of the Lorenz 96 and their indices  
% Input: n: number of variables in Lorenz 96
%        F: the constant of the Lorenz 96
% Output: the true coefficients of the Lorenz 96

% Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
% Version 1, July 2017 
% Reference: arxiv link

N = (n+1)*(n+2)/2;
c_true_mat = zeros(N,n);
for optEquation = 1:n
    % True solution
    switch optEquation
        case 1 % dx1/dt =  - x(N-1) * x(N) + x(N) * x(2) - x(1) + F;
            c_true_index =[1, 2, 3*n, N-1];
       %     c_true_value = [F, -1/sqrt(3), 1/3, - 1/3];
            c_true_value = [F, -1, 1, -1];

        case 2 % dx2/dt =  - x(N) * x(1) + x(1) * x(3) - x(2) + F;
            c_true_index =[1, 3, n+4, 2*n+1];
       %     c_true_value = [F, -1/sqrt(3), 1/3, -1/3 ];
            c_true_value = [F, -1, 1, -1];

        case n % dxN/dt = xdot(N) = - x(N-2) * x(N-1) + x(N-1) * x(1) - x(N) + F;
            c_true_index = [1, n+1, 2*n, N-4];
        %    c_true_value = [F, -1/sqrt(3), 1/3, -1/3];
            c_true_value = [F, -1, 1, -1];

        otherwise
            c_true_index = [1, optEquation+1, (optEquation-2)*(2*n-optEquation+5)/2+2, (optEquation-1)*(2*n-optEquation+4)/2+3];
         %   c_true_value = [F, -1/sqrt(3), -1/3, 1/3 ];
            c_true_value = [F, -1, -1, 1 ];

    end
    c_true_mat(c_true_index(1:4),optEquation) = c_true_value(1:4);
end
