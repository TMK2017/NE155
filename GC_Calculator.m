%% Green-Corvino Calculator Script
% Calculates scalar flux


%% "Fake" Inputs
% Make Grids n by m...
% i = 1,...,n-1
% j = 1,...,m-1
n = 10;
m = 10;

% Problem geometry
% 0 < x < a
% 0 < y < b
a = 0.1; % [cm]
b = 0.1; % [cm]

% Make D, S, Sigma_A, phi
D       = 10*ones(n,m);
Sigma_A = 2*ones(n,m);
S       = 10*ones(n*m,1);
phi     = ones(n*m,1);


%% "Fake" Boundary conditions

% Vacuum bottom
j = 0;
for i = 0:n
    phi(i,j) = 0;
end
% Vacuum Top
j = m;
for i = 0:n
    phi(i,j) = 0;
end
% Vacuum Left
i = 0;
for j = 0:m
    phi(i,j) = 0;
end
% Vacuum Right
i = n;
for j = 0:m
    phi(i,j) = 0;
end


%% Solve for flux

A = zeros(n*m); % square matrix (nm x nm)



