% Lorenz 96

% Note: Store all snapshots in Xfull
%           Xfull = [x1(t0;1)            x2(t0;1)             ...  xn(t0;1);
%                    x1(t1;1)            x2(t1;1)             ...  xn(t1;1);
%                              ..............
%
%                    x1(t_n1;1)          x2(t_n1;1)           ...  xn(t_n1;1);
%                    x1(t0;2)            x2(t0;2)             ...  xn(t0;2);
%                              ..............
%
%                    x1(t_n2;2)          x2(t_n2;2)           ...  xn(t_n2;2);
%                              ..............
%
%                    x1(t0;NumIC)        x2(t0;NumIC)         ...  xn(t0;NumIC);
%                              ..............
%                    x1(t_nNumIC;NumIC)  x2(t_nNumIC;NumIC)   ...  xn(t_nNumIC;NumIC);]

%
% =========================================================================
function soln = basisPursuit_Lorenz96(Vapproximate,D,optEquation,optPolynomial,opts,sigma)
% Description: Find the coefficient of component optEquation in Lorenz 96
% by solving the basis pursuit problem via spgl1
%            min|c|_1 subject to |D*c  - Vapproximate(:,optEquation)|_2 <= sigma
% Note: Other basis pursuit algorithms: SpaRSA, Douglas-Rachford splitting,
% primal-dual, cvx...

% Input: Vapproximate: approximation of time derivative
%        D: dictionary matrix
%        optPolynomial: should be 'legendre' or 'monomial'
%        optEquation: equation need to recover 
%        opts: parameters for spgl1
%        sigma: |y-Ac|_2<= sigma

% Output: soln(Nx1): the coefficient dxk/dt = Dictionary(Xfull)*soln
% Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward
% Version 1, July 2017 
% Reference: arxiv link
%% Basis Pursuit Denoising Problem
% L1 minimization problem: min |c|_1 subject to D*C = Vfull
n = size(Vapproximate,2);
% Equation to test
Vtest = Vapproximate(:,optEquation) ;

% solve the optimization
soln = zeros(size(D,2),1);

Dnormalized = D./repmat(sqrt(sum(D.^2,1)),size(D,1),1);

[c,~,~,~] = spgl1( Dnormalized, Vtest, 0, sigma, [], opts);

c = c./(sqrt(sum(D.^2,1)))';% only do when using Dnormalized

%  c(abs(c) < (1e-5)) = 0; % optional

% rescale back to the monomial basis
c_recover_index = find(c);
if strcmp(optPolynomial,'legendre')
    for inIter = 1:length(find(c))
        indtmp = c_recover_index(inIter);
        if (indtmp ==1)
            soln(indtmp) = c(indtmp);
        else if (indtmp<=n+1) 
            soln(indtmp) = c(indtmp)*sqrt(3);
        else soln(indtmp) = c(indtmp)*3;
            end
        end
    end
end
return

