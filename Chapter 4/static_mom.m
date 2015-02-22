% File to compute charge distribution along a wire

% Author D B Davidson
% March 2003

clear;

% Define geometry:
L = 1 % length in [m]
a=0.001% wire radius in [m]
N = 5 % Number of segments
Delta = L/N
eps_0 = 8.854e-12;
voltage = 1 % [V] Potential of wire 

% Pre-allocated storage
V_vec = zeros(N,1);
I_vec = zeros(N,1);
Z_mat = zeros(N,N);

% Set up matrix equation
for mm = 1:N,
    for nn = 1:N,
      l_m  = sqrt(((mm-nn) * Delta)^2 + a^2);  
      d_plus_mn  =  l_m + Delta/2;
      d_min_mn  =  l_m - Delta/2;
      if mm==nn
          Z_mat(mm,nn)  = 2 *  log ((Delta/2 + sqrt(a^2+(Delta/2)^2))/ a );
      elseif (abs(mm-nn) <= 2) 
          Z_mat(mm,nn)  = log ( ( d_plus_mn + sqrt((d_plus_mn)^2 + a^2) ) / ...
                         (d_min_mn + sqrt((d_min_mn)^2 + a^2)) );
      else 
          Z_mat(mm,nn)  = log ( d_plus_mn  / d_min_mn  ) ;
      end
    end
end

V_vec = 4*pi*eps_0*voltage*ones(N,1); 

I_vec = Z_mat\V_vec;

z_axis = [Delta/2:Delta:(N-1/2)*Delta];

plot(z_axis,I_vec);
xlabel('Distance along wire [m]')
ylabel('Line charge density [C/m]')
pause;
disp('Press any key to continue.')

% Now generate a plot showing the piecewise constant
% approximation: (the Bar plot function works well here)

bar(z_axis,I_vec/1e-12,1.0,'w');
xlabel('Distance along wire [m]')
ylabel('Line charge density [C/m]')

print -deps statMM_5

clear Z_mat
save static_MoM_data