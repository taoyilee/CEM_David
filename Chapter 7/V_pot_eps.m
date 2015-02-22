% Matlab m-file which computes the scalar potential V
% as a function of distance for various values of epsilon_r
% generating plots similar to Fig. 7.12.
% Usage: V_pot_eps
%
% See Chapter 7, Section 7.5, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP 2005.
%
% Written by DB Davidson, 26 August 2003. Revised 2 Feb 2005.

clear all;

c = 2.997925e8         % speed of light in vacuum [m/s] 
freq = 10e9            % Operating frequency [Hz]
lambda_0 = c/freq      % free space wavelength. Do not confuse with lambda,
                       % real part of integration variable!
%eps_r_prime = 10              % relative permittivity of substrate
tan_delta = 0  % 0.01       % loss tangent
%eps_r = eps_r_prime *(1-i*tan_delta)
                       % complex relative permittivity
k_0 =  2*pi*freq/c     % free space wavenumber
d = 0.05*lambda_0
rho_max = 10/k_0
rho_min = 0.01/k_0
N_rho = 50; % logarithmically spaced
delta_rho = exp((log(rho_max/rho_min))/N_rho);

global k k_0 rho eps_r eps_r_prime d lambda_p 

rel_eps = [1.01 2.2 4.34 9.6];

for jj = 1:4
   
   eps_r_prime = rel_eps(jj);
   eps_r = eps_r_prime *(1-i*tan_delta);
   k =  k_0 *sqrt(eps_r);  % wavenumber in dielectric
   
  % Find pole of D_TM 
  lambda_p = root_D_TM(k_0,k,eps_r,d,50)

  rho = rho_min;  
  for kk = 1:N_rho
      rho*k % for screen feedback
      rho_vec(kk) = rho;
      V_potential(jj,kk) = V_int;
      rho = rho*delta_rho;
  end
end 
norm = abs(V_potential(4,1)); % Normalize to high epsr value.
loglog(rho_vec*k_0,abs(V_potential(1,:))/norm,'k-',...
       rho_vec*k_0,abs(V_potential(2,:))/norm,'k:',...
       rho_vec*k_0,abs(V_potential(3,:))/norm,'k-.',...
       rho_vec*k_0,abs(V_potential(4,:))/norm,'k--')
xlabel('k_0 \rho')
ylabel('|V_n|')
legend(strcat('\epsilon_R=',num2str(rel_eps(1))),...
       strcat('\epsilon_R=',num2str(rel_eps(2))),...
       strcat('\epsilon_R=',num2str(rel_eps(3))),...
       strcat('\epsilon_R=',num2str(rel_eps(4))),0)
print -deps Som_Veps
pause;


