% Matlab m-file which computes the current and input impedance
% for a thin strip dipole on a grounded substrate using the 
% Mixed Potential Integral Equation 
% and entire domain cosinusoidal basis functions. 
% The substrate must be sufficiently thin that only the dominant TM mode can 
% propagate.
% Printed dipole  based on Pozar et al, IEEE T-AP, June 1984, p.602ff - 
% but 0.12 rather than 0.19 lambda thick.
%
% See Chapter 7, Section 7.6, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP 2005 and 2nd ed, CUP 2011.
%
% D.B.Davidson: 10 September 2003

clear all;

c = 2.997925e8           % speed of light in vacuum [m/s] 
eps_0 = 8.854e-12
mu_0 = pi*4e-7
freq_res = 10e9          % resonance (design) frequency
lambda_res = c/freq_res  % wavelength at resonance frequency
d = 0.12*lambda_res      % substrate height
L = 0.39*lambda_res      % strip dipole length (width << L, not required).
freq = 9e9              % operating frequency [Hz]
eps_r_prime = 2.55       % relative permittivity of substrate
tan_delta =   0          % loss tangent - must be zero at present.
eps_r = eps_r_prime *(1-i*tan_delta)
                         % complex relative permittivity

max_mode_index = input('Enter maximum mode index')  % approx half max mode number (=2*max_mode_index-1).
num_int        = input('Enter number of integration points - must be even.')

global k k_0 rho eps_r eps_r_prime d lambda_p 


del_x = L/num_int % size of each interval
x_pot = [-L/2+del_x/2:del_x:L/2-del_x/2];

freq_vec = [0.7:0.025:1.1]*freq_res
%freq_vec = 0.9*freq_res
for ifreq = 1:length(freq_vec)
  freq = freq_vec(ifreq)
  omega = 2*pi*freq        % angular frequency
  lambda_0 = c/freq        % free space wavelength. Do not confuse with lambda,
                           % real part of integration variable!
  k_0 =  2*pi*freq/c       % free space wavenumber
  k =  k_0 *sqrt(eps_r);   % wavenumber in dielectric
  % Find pole of D_TM 
  lambda_p = root_D_TM(k_0,k,eps_r,d,50)

  % Pre-compute interpolation table data.
  % Special treatement to avoid rho=0
  rho = abs(del_x/4); % Value less than smallest distance between offset points 
                      % so that interpolation lies inside table limits. 
  rho_pot(1) = rho;
  Vpot_table(1)   = V_int;
  for jj = 2:num_int % field point index
    rho = abs(x_pot(jj)-x_pot(1));
    rho_pot(jj) = rho;
    Vpot_table(jj)   = V_int;
  end
  % and special "extra point" at end of interpolation table, as above. 
  rho = L;
  rho_pot(num_int+1) = rho;
  Vpot_table(num_int+1)   = V_int;

  % Quadrature points - offset .
  x = [(-L/2+del_x/3):del_x:(L/2-del_x*2/3)];
  xp = [(-L/2+2*del_x/3):del_x:(L/2-del_x/3)];

  % Pre-compute potentials - not functions of modes. 
  for jj = 1:num_int % field point index
    for kk = 1:num_int % source point index
       rho = abs(x(jj)-xp(kk));
       R_0 = rho;
       R_1 = sqrt(rho^2+(2*d)^2);
       Apot(jj,kk) = exp(-j*k_0*R_0)/R_0 - exp(-j*k_0*R_1)/R_1;
       Vpot(jj,kk)   = interp1(rho_pot,Vpot_table,rho);
       if isnan(Vpot(jj,kk))
         jj
         kk
         input('Error in precomputation of potentials, hit break')
       end
    end 
  end

  for m_index = 1:max_mode_index
    mm = 2*m_index-1;
    for n_index = 1:max_mode_index
      nn = 2*n_index-1;
      for jj = 1:num_int % field point index 
        for kk = 1:num_int % source point index
            rho = abs(x(jj)-xp(kk));
            amn_inner(kk) = cos(nn*pi*xp(kk)/L)*Apot(jj,kk);
            vmn_inner(kk) = sin(nn*pi*xp(kk)/L)*Vpot(jj,kk);
        end 
        amn_outer(jj) = cos(mm*pi*x(jj)/L) * del_x*trapz(amn_inner);
        vmn_outer(jj) = sin(mm*pi*x(jj)/L) * del_x*trapz(vmn_inner);
      end 
      A_mat(m_index,n_index) = del_x*trapz(amn_outer);
      V_mat(m_index,n_index) = (mm*nn*pi^2)/L^2*del_x*trapz(vmn_outer);
    end         
  end
  A_mat  = j*omega*(mu_0/(4*pi))*A_mat;
  V_mat  = 1/(j*omega)*(1/(2*pi*eps_0))*V_mat; 
  % For debugging:
 % A_matrix(ifreq) = A_mat
 % V_matrix(ifreq) = V_mat
  
  % End debugging:
  for m_index = 1:max_mode_index
    b_vec(m_index,1) = 1; 
  end 

  Z_mat = A_mat+V_mat;
  I_vec = Z_mat\b_vec;

  % Construct solution at sample points x_s. 
  x_s = [-L/2+del_x/2:del_x:L/2-del_x/2];

  I_s = zeros(1,length(x_s));
  for m_index = 1:max_mode_index
    mm = 2*m_index-1;
    I_s = I_s+I_vec(m_index)*cos(mm*pi*x_s/L);
  end

  Zin(ifreq) = 1/I_s(num_int/2)
%   plot(x_s/L,real(I_s)/1e-3,x_s/L,imag(I_s)/1e-3)
%   xlabel('x/L')'
%   ylabel('I [mA]')
%   legend('Real','Imag',0)
%   print -deps MoM_Som
end
plot(freq_vec/1e9,real(Zin),...
     freq_vec/1e9,imag(Zin),'--')
grid on;
xlabel('Freq[GHz]')'
ylabel('Impedance')
legend('Real','Imag',0)
print -deps MoM_Som_freq

