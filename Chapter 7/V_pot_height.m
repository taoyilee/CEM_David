% Matlab m-file which computes the scalar potential V
% as a function of distance for various values of normalized height,
% generating plots similar to Fig. 7.10 & 7.11.
% Usage: V_pot_height
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
eps_r_prime = 10              % relative permittivity of substrate
tan_delta = 0 % 0.01       % loss tangent
eps_r = eps_r_prime *(1-i*tan_delta)
                       % complex relative permittivity
k_0 =  2*pi*freq/c     % free space wavenumber
k =  k_0 *sqrt(eps_r)  % wavenumber in dielectric
rho_max = 10/k
rho_min = 0.01/k
N_rho = 100; % logarithmically spaced
delta_rho = exp((log(rho_max/rho_min))/N_rho);


global k k_0 rho eps_r eps_r_prime d lambda_p 

b = [0.8 0.5 0.2]
  

for jj = 1:3
   d = b(jj)*pi/(2*k_0*sqrt(eps_r_prime-1))  
   
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
norm = abs(V_potential(1,1)) % For small rho, values are asymptotically similar. 
loglog(rho_vec*k,abs(V_potential(1,:))/norm,'k-',...
       rho_vec*k,abs(V_potential(2,:))/norm,'k:',...
       rho_vec*k,abs(V_potential(3,:))/norm,'k-.')
xlabel('k \rho')
ylabel('|V_n|')
legend(strcat('b=',num2str(b(1))),strcat('b=',num2str(b(2))),...
       strcat('b=',num2str(b(3))),0)
print -deps Som_Vh
pause;

semilogx(rho_vec*k,angle(V_potential(1,:))*180/pi,'k-',...
         rho_vec*k,angle(V_potential(2,:))*180/pi,'k:',...
         rho_vec*k,angle(V_potential(3,:))*180/pi,'k-.')
xlabel('k \rho')
ylabel('Phase angle V_n')
axis([k*rho_min k*rho_max -180 +180])
legend(strcat('b=',num2str(b(1))),strcat('b=',num2str(b(2))),...
       strcat('b=',num2str(b(3))),3)
print -deps Som_Vhpz



