function [Xfull,Vapproximate,Vexact] = Lorenz96_XV(F,Xint,dt,SizeOfBurst)
% Decription: Construct all-in-one including the data matrix, approximation
% velocity, and exact velocity given F, initialization, timestep, and
% length of the evolution.

% Input: dt: time step
%        SizeOfBurst: size of each burst
%        NumIC: number of initializations
% Output: Xfull of size Mxn
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
%      Vapproximate = approximate derivative from data, using finite difference
%      Vexact = lorenz96(Xfull)

% Copywright: Hayden Schaeffer, Giang Tran, and Rachel Ward.
% Version 1, July 2017 
% Reference: arxiv link

n = size(Xint,1);
NumIC = size(Xint,2);
Tfinal = (SizeOfBurst-1)*dt; % Small Burst

Xfull = []; 
Vapproximate = [];
Vexact = [];

for iter = 1:NumIC
    Xinttmp  = Xint(:,iter);
    % Solve Lorenz 96
    [~,Xtmp] = ode45(@(t,x) lorenz96(t,x,n,F),[0:dt:Tfinal], Xinttmp); % each row is solution at a specific time

    % Store data
    Xfull = [Xfull;Xtmp]; % all snapshots

     % Exact time-derivative by evaluating the RHS at data
    for i=1:SizeOfBurst
        Vtmpexact(i,:) =  lorenz96(Tfinal,Xtmp(i,:),n,F); %exact velocity Vtmp is (SizeOfBurst x N)
    end
    Vexact = [Vexact;Vtmpexact];

    % Approximate time-derivative using finite difference
    Vtmp = time_derivative(Xtmp,dt,2);
    Vapproximate = [Vapproximate;Vtmp];
end

