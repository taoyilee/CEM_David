function [I_vec,phi_c,w,x_c,y_c,cond_num] = MoM_TM_solver(k,N,a,E_0,phi_inc,quad,eta,toeplitz_flag);
% MoM solution for TM incidence on PEC conductor.

% Author: DB Davidson, 22 Feb 2008. 

% k: wavenumber
% N: number of segments
% a: radius of PEC cylinder
% E_0: amplitude of incident E field
% phi_inc: angle incident field incident from.
% quad: input flag to select quadrature or not, see below.
% toeplitz_flag: input flag to exploit Toeplitz symmetry

% Returns vector of complex-valued currents and corresponding real-valued angles at
% segment centres, as well as some geometrical information and the
% condition number of the matrix.

% The MoM solution uses piecewise constant basis functions, and collocation
% (point matching) for the EFIE in E_z. The matrix elements are approximately evaluated, either  
% using a single-point evaluation or a low-order quadrature scheme,
% depending on the input flag QUAD.

% Pre-allocated storage
V_vec = zeros(N,1);
I_vec = zeros(N,1);
Z_mat = zeros(N,N);

Delta_phi = 2*pi/N;  % \delta phi, in radians
% Pre-compute angles of segment mid-points (note that the 1st segment
% is symmetrical about 0 degrees, from -Delta_phi/2 to +Delta_phi/2.
for mm = 1:N,
    phi_c(mm)  = (mm-1)*Delta_phi;      % Angles of mid-points of segments.
    Omega(mm)  = pi/2+(mm-1)*Delta_phi; % Orientation angles of sements (Right-hand rule for current). Only needed for quadrature evaluation.
end
w   = 2*a*tan(Delta_phi/2); % width of strips (assumed constant) (Note that this formula implies that phase centres lie on circle, 
                            % ie polygon circumscribes circular cylinder).
x_c = a*cos(phi_c);         % Coordinates of mid-points of segments
y_c = a*sin(phi_c);

if toeplitz_flag
   M_stop = 1;
else 
   M_stop = N;
end

% Set up matrix equation
for mm = 1:M_stop,
    for nn = 1:N,
      if mm==nn
          gamma = 1.781072418;
          Z_mat(mm,nn)  = k*eta*w/4*(1-j*2/pi*(log(gamma*k*w/4)-1));
      else
          if ~quad % single point evaluation, eq. 2.10
              R_mn = sqrt( (x_c(mm)-x_c(nn))^2 + (y_c(mm)-y_c(nn))^2 );
              Z_mat(mm,nn)  = k*eta/4*w*besselh(0,2,k*R_mn); % Eq.2.30, approximate result.
          else
              % evaluation of 2.8 by quadrature. MATLAB quad routine gives problems, done by hand with zero order rule.
              N_quad = 10; % Number of integration points.
              t_start = 0; % Start value of parameter t along line segment.
              delta_t = w/N_quad;
              %               Z_mat(mm,nn) = 0;
              for ii = 1:N_quad
                  t = t_start+(ii-1/2)*delta_t;
                  Z_mat(mm,nn) = Z_mat(mm,nn) + EFIEkernel(t,a,Omega(nn),phi_c(nn),Delta_phi,w,x_c(mm),y_c(mm),k)*delta_t;
              end
              Z_mat(mm,nn) = k*eta/4*Z_mat(mm,nn);
          end         
      end
    end
end

if toeplitz_flag
    Z_mat = toeplitz(Z_mat(1,1:N),Z_mat(1,1:N));
    % Note that naive use as toeplitz(Z_mat(1,1:N)) produces matrix with Hermitian symmetry, which is incorrect.
else
    %    disp(Z_mat)
end

for mm = 1:N,
    V_vec(mm,1) = E_0*exp(j*k*(x_c(mm)*cos(phi_inc) + y_c(mm)* sin(phi_inc))) ;
end

I_vec = Z_mat\V_vec; % Solve linear system 

cond_num = cond(Z_mat);



function [kernel] = EFIEkernel(t,a,Omega_n,phi_c,Delta_phi,w,x_c_m,y_c_m,k)
tangent = [cos(Omega_n) sin(Omega_n)];
phi_1 = phi_c-Delta_phi/2;
a1 = sqrt(a^2+(w/2)^2);
x = a1*cos(phi_1) +tangent(1)*t;
y = a1*sin(phi_1) +tangent(2)*t;
R_m = sqrt( (x_c_m-x)^2 + (y_c_m-y)^2 );
kernel  = besselh(0,2,k*R_m); % Eq.2.8, exact kernel.
