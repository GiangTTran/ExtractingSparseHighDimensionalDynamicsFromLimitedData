%==========================================================================
% Description: compute time_derivative using 1st/2nd numerical approximation of time derivative 
%       1st-order approximation:
%                               roc(x(t)) = (x(t+dt) - x(t))/dt
%       2nd-order approximation:
%                               roc(x(t)) = (x(t+dt) - x(t-dt))/(2*dt)
%
% Input: 
%       X(mxn): recorded data, kth row is the measurement value at time k*dt 
%       dt:     time step
%                    where m = number of measurements
%                          n = dimension of the ODE system
%       type: 1 or 2 (1st/2nd order approximation of time derivative)
% Output: 
%       roc( mxn: 1st-order approximation of time derivative
%       roc( mxn: 2nd-order approximation of time derivative
%
% Remark 1: The code can be generalized for non-equal time distribution. For
% example, the 2nd-order approximation will be 
%       roc(x(t_k)) = (x(t_{k+1}) - x(t_{k-1})) / (t_{k+1} - t_{k-1});
%
% Remark 2: For type 2, include derivative estimation at every time
%  dx1/dt = (x2 - x1)/dt  % forward
%  dx2/dt = (x3 - x1)/(2*dt) % central
%   ...
%  dx_{n-1}/dt = (x_n - x_{n-2}) /(2*dt) % central
%  dx_n/dt = (x_n - x_{n-1})/dt % backward

% Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
% Version 1, July 2017 
% Reference: arxiv link

%==========================================================================
function roc = time_derivative(X,dt,type)
[m,n] = size(X);
if (type==1)
    roc = zeros(m,n);
    roc = (X(2:m,:) - X(1:m-1,:))/ dt; % forward Euler
    roc(m,:) = (X(m,:) - X(m-1,:))/dt; % backward
else if (type==2)
        roc = zeros(m,n);
        roc(1,:) = (X(2,:) - X(1,:))/dt; % forward
        
        roc(2:m-1,:) = (X(3:m,1:n) - X(1:m-2,1:n)) /(2*dt); % note roc(k,:) is the time derivative at time t_{k+1}
        
        roc(m,:) = (X(m,:) - X(m-1,:))/dt; % backward
    end
end
return