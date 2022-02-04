
%Load initial stable solution. Should contain the following values:
%   alpha_0, F_0 - values of alpha and F for the initial solution
%   L - mode circumference of the solution
%   N - discretization of solution
%   uout - solution
%   See README for more details



%set bounds on alpha, alim_h should not be extended past alpha = 6
alim_h = 6;
alim_l = -2;

%set step size for alpha and F
da = 0.01; dF = 0.0005;

%set initial direction to move in alpha and F, 1 - positive, -1 - negative
adir = 1; Fdir = -1;

MyPath = pwd;
AddedPath = [MyPath '/functions/'];
addpath(AddedPath, '-end');

[alpha_bound,F_bound,sol,end_flag] = boundary_wrapper(alpha_0,da,adir,F_0,dF,Fdir,L,N,uout,alim_h,alim_l,1);
if end_flag == 2
    [alpha_bound2,F_bound2,~,~] = boundary_wrapper(alpha_bound(1),da,-adir,F_bound(1),dF,Fdir,L,N,...
                                                    sol,alim_h,alim_l,2);
    alpha_bound = [flip(alpha_bound2(2:end)) alpha_bound];
    F_bound = [flip(F_bound2(2:end)) F_bound];
end

%quick plot of boundary
figure
plot(alpha_bound,F_bound,'Linewidth',2)
xlim([alim_l alim_h])
box on
