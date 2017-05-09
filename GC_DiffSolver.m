function [phi] = GC_DiffSolver(x,y,D,sigma,S)
% Green-Corvino Calculator function
% Calculates scalar flux


%% Set n,m, epsilon and delta
n = length(x); % in this solver n = n+1 (from math derivation) and m (from code) = m+1 (from math)
m = length(y);
epsil = y(2:end) - y(1:end-1);
del = x(2:end) - x(1:end-1);

%% Pad S, sigma, and V matrices w/ zeros
%Pad S with zeros to account for edge cases (boundaries + 4 corners)
S_expand = zeros(n+1,m+1);
S_expand(2:end-1,2:end-1) = S;

%Pad sigma with zeros to account for edge cases (boundaries + 4 corners)
Sig_expand = zeros(n+1,m+1);
Sig_expand(2:end-1,2:end-1) = sigma;

% Create V (pad with zeros to account for edge cases)
V = zeros(n+1,m+1);
V(2:end-1,2:end-1) = del'*epsil./4;





%% Create S volume average (S_avg)
% S*ij = SijVij + Si+1jVi+1j + Sij+1Vij+1 + Si+1j+1Vi+1j+1
% S*ij = S_A    + S_B        + S_C        + S_D
S_A = S_expand(1:end-1,1:end-1).*V(1:end-1,1:end-1);
S_B = S_expand(2:end,1:end-1).*V(2:end,1:end-1);
S_C = S_expand(1:end-1,2:end).*V(1:end-1,2:end);
S_D = S_expand(2:end,2:end).*V(2:end,2:end);
S_avg = S_A + S_B + S_C + S_D;

% Boundary conditions:
S_avg(1,:) = 0;
S_avg(:,1) = 0;


%% Create sigma volume average (sigma_avg)
% Sig*ij = Sig_ijVij + Sig_i+1jVi+1j + Sig_ij+1Vij+1 + Sig_i+1j+1Vi+1j+1
% Sig*ij = Sig_A    + Sig_B        + Sig_C        + Sig_D
Sig_A = Sig_expand(1:end-1,1:end-1).*V(1:end-1,1:end-1);
Sig_B = Sig_expand(2:end,1:end-1).*V(2:end,1:end-1);
Sig_C = Sig_expand(1:end-1,2:end).*V(1:end-1,2:end);
Sig_D = Sig_expand(2:end,2:end).*V(2:end,2:end);
sigma_avg = Sig_A + Sig_B + Sig_C + Sig_D;



%% Create A matrix (solves A*phi = S_avg)

% Pad del and epsil with a leading & trailing zero
del = [0;del';0];
epsil = [0,epsil,0];

% Pad D with zeros
D_expand = zeros(n+1,m+1);
D_expand(2:end-1,2:end-1) = D;


%% Calculate A_B (bottom)
% Calculate aij_i,j-1
%vector_B = -reshape((D_expand(1:end-1,1:end-1).*del(1:end-1)' + D_expand(2:end,1:end-1).*del(2:end))./(2*epsil(1:end-1)'),[],1)
B_grid = -(D_expand(1:end-1,1:end-1).*repmat(del(1:end-1),1,m) + D_expand(2:end,1:end-1).*repmat(del(2:end),1,m))./(2*repmat(epsil(1:end-1),n,1));



% Boundary Conditions 
B_grid(1,:)  = 0;   % left   (phi=0, i=0,j=0...m)
B_grid(:,1)  = 0; % bottom (phi=0, i=0...n,j=0)


% convert inf / nan values to zero
index = isinf(B_grid) | isnan(B_grid);
B_grid(index) = 0;

% convert to matrix
vector_B= reshape(transpose(B_grid),[],1);
A_B = spdiags(vector_B,0,n*m,n*m);



%% Calculate A_T (top)
% Calculate aij_i,j+1
%vector_T = -reshape((D_expand(1:end-1,2:end).*del(1:end-1)' + D_expand(2:end,2:end).*del(2:end)')./(2*epsil(2:end)'),[],1);
T_grid= -(D_expand(1:end-1,2:end).*repmat(del(1:end-1),1,m) + D_expand(2:end,2:end).*repmat(del(2:end),1,m))./(2*repmat(epsil(2:end),n,1));

% Boundary Conditions 
T_grid(:,1) = 0; % (phi=0, i=0...n,j=0)
T_grid(1,:) = 0; % (phi=0, i=0,j=0...m)

% convert inf / nan values to zero
index = isinf(T_grid) | isnan(T_grid);
T_grid(index) = 0;

% convert to matrix
vector_T= reshape(transpose(T_grid),[],1);
A_T = spdiags(vector_T,0,n*m,n*m);



%% Calculate A_L (left)
% Calculate aij_i-1,j
%vector_L = -reshape((D_expand(1:end-1,1:end-1).*epsil(1:end-1) + D_expand(1:end-1,2:end).*epsil(2:end))./(2*del(1:end-1)),[],1);
L_grid = -(D_expand(1:end-1,1:end-1).*repmat(epsil(1:end-1),n,1) + D_expand(1:end-1,2:end).*repmat(epsil(2:end),n,1))./(2*repmat(del(1:end-1),1,m));


% Boundary Conditions 
L_grid(:,1) = 0; % (phi=0, i=0...n,j=0)
L_grid(1,:) = 0; % (phi=0, i=0,j=0...m)

% convert inf / nan values to zero
index = isinf(L_grid) | isnan(L_grid);
L_grid(index) = 0;

% convert to matrix
vector_L= reshape(transpose(L_grid),[],1);
A_L = spdiags(vector_L,0,n*m,n*m);



%% Calculate A_R (right)
% Calculate aij_i+1,j
%vector_R = -reshape((D_expand(2:end,1:end-1).*epsil(1:end-1) + D_expand(2:end,2:end).*epsil(2:end))./(2*del(2:end)),[],1);
R_grid = -(D_expand(2:end,1:end-1).*repmat(epsil(1:end-1),n,1) + D_expand(2:end,2:end).*repmat(epsil(2:end),n,1))./(2*repmat(del(2:end),1,m));

% Boundary Conditions 
R_grid(:,1) = 0; % (phi=0, i=0...n,j=0)
R_grid(1,:) = 0; % (phi=0, i=0,j=0...m)

% convert inf / nan values to zero
index = isinf(R_grid) | isnan(R_grid);
R_grid(index) = 0;

% convert to matrix
vector_R= reshape(transpose(R_grid),[],1);
A_R = spdiags(vector_R,0,n*m,n*m);


%% Calculate A_C (center)
A_C = spdiags(sigma_avg(:),0,n*m,n*m) - A_B - A_T - A_L - A_R;


%% Calculate A (combine effects from top, bottom, left, right, center)

A = ...
    A_C +...
    spdiags(vector_R(1:end-1),-1,n*m,n*m)' +...
    spdiags(vector_L(2:end),-1,n*m,n*m)  +...
    spdiags(vector_T(1:end-m),-m,n*m,n*m)' +...
    spdiags(vector_B(1+n:end),-n,n*m,n*m);



%% Solve A*phi = S_avg using backslash operator
S_avg = S_avg';
phi = A\S_avg(:);

% Rearrange phi into a (nxm) matrix
phi = reshape(phi,[m,n])';
end

