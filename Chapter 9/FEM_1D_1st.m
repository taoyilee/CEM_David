% A 1D FEM solution of the wave equation using 1st order elements.
% The model problem is a transmission line, with the source at z=0 and
% an open circuit at z=ell.

% To run without interruption, use pause off

% Author: David B Davidson, Stellenbosch University, July 2009. Developed for
% "Computational Electromagnetics for RF and Microwave Engineering",
% Cambridge Univ Press, 2nd edition, to be published.

clear;
close all;

% Setup transmission line parameters. At present, tx line is modelled as of
% uniform inductance and capacitance per unit length (L and C
% respectively). In this script, both are set to 1, resulting in  c=1m/s and Z_0 = 1 Ohm.
L = 1;
C = 1;
c= 1/sqrt(L*C);
Z_0 = sqrt(L/C);
f = 1; % Hz
omega= 2*pi*f;
lambda = c/f;
ell = lambda/2;
V_in = 1; % Voltage at source end of tx line.

% Setup FEM mesh. N_elem is number of elements (hence number of nodes is N_elem+1), h is resulting mesh length.

N_elem=input('Enter initial number of elements (default 2)');
if isempty(N_elem)
    N_elem= 2;
end

num_meshes=input('Enter number of mesh refinement stages (default 1)');
if isempty(num_meshes)
    num_meshes= 1;
end

for nn=1:num_meshes
    % Call solver
    h(nn) = ell/N_elem;
    V = FEM_1D_solver(N_elem,h(nn),omega,L,C,V_in);
    %    plot(linspace(0,ell,N_elem+1),V(:,nn));
    N_int = N_elem*10;
    [V_fem,z] = FEM_pp(V,N_elem,ell,h(nn),N_int);
    beta = 2*pi/lambda;
    V_anl = -cos(2*beta*z*ell);
    plot(z,V_fem,'-',z,V_anl,'-.');
    axis([0,ell,-1,1])
    legend('FEM','Analytical','Location','Best')
    pause
    err = (V_anl-V_fem); % Note - as the analytical solution is zero in places, this cannot be normalized.
    rms_err(nn) = norm(err)/sqrt(length(V_anl));
    N_elem = N_elem*2;
end
if num_meshes >1 % Convergence plot only makes sense with more than one mesh
    figure;
    loglog(h/lambda,rms_err)
    xlabel('h/\lambda')
    ylabel('RMS error')
    p=polyfit(log(h),log(rms_err),1)
    slope=p(1)
end