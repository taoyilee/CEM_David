function [V_tot] = FEM_1D_solver(N_elem,h,omega,L,C,V_in)
% A function to solve the 1D wave equation on a transmission line.
% It is assumed that the line is uniformly discretized, and that 
% the tranmsission line parameters L and C are constant with length.

N = N_elem+1; % number of nodes.
% Construct  M: 
S = zeros(N); % Pre-allocate matrices.
T = zeros(N); 
for ii=1:N;
    S(ii,ii) = 2;
    T(ii,ii) = 4;
    if ii > 1
        S(ii,ii-1) = -1;
        T(ii,ii-1) = 1;
    end
    if ii < N
        S(ii,ii+1) = -1;
        T(ii,ii+1) = 1;
    end
end
S(1,1) = 1; % Set first and last diagonal elements
S(N,N) = 1;
S = 1/(h*L)*S; % Scale entries
T(1,1) = 2; 
T(N,N) = 2;
T = -omega^2*h*C*T/6;
M = zeros(N); % Pre-allocate matrix.
M = S+T;
F = zeros(N-1,1); % Forcing vector
F(N-1) = -M(N-1,N)*V_in; % Voltage at source end.
V = M(1:N-1,1:N-1)\F;

% Post-process
V_tot = zeros(N,1);
V_tot(1:N-1) = V;
V_tot(N) = V_in;
