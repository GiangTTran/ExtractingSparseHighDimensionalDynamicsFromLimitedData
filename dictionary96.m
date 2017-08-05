function phiX = dictionary96(U,option)
% Description: Construct the dictionary matrix phiX containing all multivariate monomials up to degree two for the Lorenz 96
% Input: U = [x1(t1) x2(t1) .... xn(t1)
%             x1(t2) x2(t2) .... xn(t2)
%                    ......
%             x1(tm) x2(tm) .... xn(tm)]
%        option = [] (monomial) or 'legendre'
% Output: the dictionary matrix phiX of size m by N, where m= #measurements and N = (n^2+3n+2)/2
% Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
% Version 1, July 2017 
% Reference: arxiv link

if (nargin == 1)
    option = [];
end
m = size(U,1); % number of measurements
n = size(U,2); % dimension of the ODE

phiX = zeros(m,(n+1)*(n+2)/2);% 1 + n + n*(n+1)/2
% 1 - 1 column
phiX(:,1) = ones(m,1);% 1
% X - n columns
phiX(:,2:n+1) = sqrt(3)*U;               % x1 x2 ... xn
% X^2 - n*(n+1)/2 columns

for k = 1:n
        phiX(:, (k*(2*n-k+3)/2+1) : ((k+1)*n -k^2/2 + k/2 +1 )) = 3*repmat(U(:,k),1,n+1-k).*U(:,k:n);
    if (strcmp(option,'legendre'))
        phiX(:, (k*(2*n-k+3)/2+1)) = (sqrt(5)/2.0)*(3*U(:,k).^2-ones(m,1));
    end
    
end
