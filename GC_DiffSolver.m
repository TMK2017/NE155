function [ phi ] = GC_DiffSolver( n,m )
% Green-Corvino Calculator function
% Calculates scalar flux


%% "Fake" Inputs
% Make Grids n by m...
% i = 1,...,n-1
% j = 1,...,m-1
% Problem geometry
% 0 < x < a
% 0 < y < b
a = 0.1; % [cm]
b = 0.1; % [cm]

% Make D, S, Sigma_A, phi
D       = 10*ones(n,m);
Sigma_A = 2*ones(n,m);
S       = 10*ones(n*m,1);
phi     = 10*ones(n,m);






end

