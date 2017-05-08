function [ phi ] = GC_DiffSolver()
% Green-Corvino Calculator function
% Calculates scalar flux




%% "Fake" Inputs
% Make Grids n by m...
% i = 1,...,n-1
% j = 1,...,m-1
% Problem geometry
% -a < x < a
% -b < y < b
a = 1;
b = 1;
n = 4;
m = 4;


x = linspace(-a,a,n);
y = linspace(-b,b,m);
eps = y(2:end) - y(1:end-1);
del = x(2:end) - x(1:end-1);

% Make D, S, Sigma_A, phi
D       = 10*ones(n-1,m-1);
sigma = 2*ones(n-1,m-1);
S       = 10*ones(n-1,m-1);

%% Pad S, sigma, and V matrices w/ zeros
%Pad S with zeros to account for edge cases (boundaries + 4 corners)
S_expand = zeros(n+1,m+1);
S_expand(2:end-1,2:end-1) = S;

%Pad sigma with zeros to account for edge cases (boundaries + 4 corners)
Sig_expand = zeros(n+1,m+1);
Sig_expand(2:end-1,2:end-1) = sigma;

% Create V (pad with zeros to account for edge cases)
V = zeros(n+1,m+1);
V(2:end-1,2:end-1) = del'*eps./4;





%% Create S volume average (S_avg)
% S*ij = SijVij + Si+1jVi+1j + Sij+1Vij+1 + Si+1j+1Vi+1j+1
% S*ij = S_A    + S_B        + S_C        + S_D
S_A = S_expand(1:end-1,1:end-1).*V(1:end-1,1:end-1);
S_B = S_expand(2:end,1:end-1).*V(2:end,1:end-1);
S_C = S_expand(1:end-1,2:end).*V(1:end-1,2:end);
S_D = S_expand(2:end,2:end).*V(2:end,2:end);
S_avg = S_A + S_B + S_C + S_D;




%% Create sigma volume average (sigma_avg)
% Sig*ij = Sig_ijVij + Sig_i+1jVi+1j + Sig_ij+1Vij+1 + Sig_i+1j+1Vi+1j+1
% Sig*ij = Sig_A    + Sig_B        + Sig_C        + Sig_D
Sig_A = Sig_expand(1:end-1,1:end-1).*V(1:end-1,1:end-1);
Sig_B = Sig_expand(2:end,1:end-1).*V(2:end,1:end-1);
Sig_C = Sig_expand(1:end-1,2:end).*V(1:end-1,2:end);
Sig_D = Sig_expand(2:end,2:end).*V(2:end,2:end);
sigma_avg = Sig_A + Sig_B + Sig_C + Sig_D;



%% Create A matrix (solves A*phi = S_avg)

% Pad del and eps with a leading & trailing zero
del = [0,del,0];
eps = [0,eps,0];

% Pad D with zeros
D_expand = zeros(n+1,m+1);
D_expand(2:end-1,2:end-1) = D;

%% Calculate A_B (bottom)
A_B = -diag(...
    reshape(...
    (D_expand(1:end-1,1:end-1).*del(1:end-1) + D_expand(2:end,2:end).*del(2:end))./(2*eps(1:end-1))...
    ,[],1));
% convert inf / nan values to zero
index = isinf(A_B) | isnan(A_B);
A_B(index) = 0;


%% Calculate A_T (top)
A_T = -diag(...
    reshape(...
    (D_expand(1:end-1,1:end-1).*del(1:end-1) + D_expand(2:end,2:end).*del(2:end))./(2*eps(1:end-1))...
    ,[],1));
% convert inf / nan values to zero
index = isinf(A_T) | isnan(A_T);
A_T(index) = 0;




%% Calculate A_C (center)
A_C = diag(sigma_avg(:)) - A_B - A_T
%A = A_C
A = eye(n*m);





%% Solve A*phi = S_avg using backslash operator
phi = A\S_avg(:);

% Rearrange phi into a (nxm) matrix
phi = reshape(phi,[n,m]);

%phi = zeros(n,m);
%phi(:) = phi_vector;




end

