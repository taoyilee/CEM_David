% Matlab m-file to plot integrand of scalar potential V for HED, generating
% plots similar to Fig. 7.4, 7.5, 7.6, 7.7 and 7.8 (and some additional
% figures).
%
% See Chapter 7, Section 7.5, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP 2005.

% Based on Mosig, "Integral equation technique", in 
% "Numerical techniques for microwave and millimeter-wave
% passive structures", Wiley 1989. 
% Note error in [Mosig]: In expression for subsrate thickness h (d in this script)
% k_0 h =0.2*pi the pi has been omitted. 

% Written by DB Davidson, August 2003. Revised 2 Feb 2005.

clear all;

c = 2.997925e8         % speed of light in vacuum [m/s] 
freq = 10e9            % Operating frequency [Hz]
lambda_0 = c/freq      % free space wavelength. Do not confuse with lambda,
                       % real part of integration variable!
eps_r_prime =  5       % relative permittivity of substrate
% tan_delta = 0.01       % loss tangent
tan_delta = 0.00       % loss tangent
eps_r = eps_r_prime *(1-i*tan_delta)
                       % complex relative permittivity
k_0 =  2*pi*freq/c     % free space wavenumber
k =  k_0 *sqrt(eps_r) % wavenumber in dielectric
d = 0.2*pi/k_0         % to match [Fig11,p168,Mosig] - h in that reference.
rho =  3/k_0;            % radial distance, ditto

% Find pole of D_TM 
TMroot = root_D_TM(k_0,k,eps_r,d,50)

% Int. variable
lambda = [0:0.0005:5] * k_0;  % k_rho = lambda on the real axis          
k_0_vec = ones(1,length(lambda))*k_0;
k_vec = ones(1,length(lambda))*k;

% Sommerfeld parameters:
u = sqrt(lambda.^2-k_vec.^2);
u_0 = sqrt(lambda.^2-k_0_vec.^2);

plot(lambda/k_0,real(D_TM(eps_r,u_0,u,d)),lambda/k_0,imag(D_TM(eps_r,u_0,u,d)),':');
xlabel('\lambda/k_0')
ylabel('D_{TM}')
legend('Real','Imag',0)
print -deps Som_D_TM
pause;

plot(lambda/k_0,abs(D_TE(u_0,u,d)))
title('D_{TE}')
pause;

Integrand = F(lambda,rho,eps_r,u_0,u,d);

plot(lambda/k_0,real(Integrand),lambda/k_0,imag(Integrand),':');
axis([0 5 -1 2]);
grid;
xlabel('\lambda/k_0')
legend('Real part','Imag part',0)
pause;
print -deps som_int

plot(lambda/k_0,real(Integrand),lambda/k_0,imag(Integrand),'-.');
axis([0.9 1.4 -0.6 0.2]);
grid;
xlabel('\lambda/k_0')
legend('Real part','Imag part',0)
pause;
print -deps som_intz

% -------------Region 1: lambda = [0,k_0] 
%   Change of variable lambda = k_0 cos(t)
clear u u_o k_0_vec k_vec Integrand 
t = [0:0.01:pi/2]; % interval lambda = Re(k_rho) = [0,k_0] on real axis

k_0_vec = ones(1,length(t))*k_0;
k_vec = ones(1,length(t))*k;
u = sqrt(k_0*cos(t).^2-k_vec.^2);
u_0 = sqrt(k_0*cos(t).^2-k_0_vec.^2);

Integrand_t_reg1 = F(k_0*cos(t),rho,eps_r,u_0,u,d)*k_0.*sin(t);

plot(t,real(Integrand_t_reg1),t,imag(Integrand_t_reg1),'-.');
%axis([0.0 2*pi -1 1]);
grid;
xlabel('t')
legend('Real part','Imag part',0)
pause;

% -------------Region 2: lambda = [k_0,k_0*sqrt(eps_r_prime)]
%   Extraction of singularity.
%   Change of variable lambda = k_0 cosh(t)
clear u u_o k_0_vec k_vec 
lambda_p = TMroot % Root as found above. 

% Crude approximation of residue.
% Interpolation through the pole gives a similar value to the correct result
% obtained below. However, the plotted values differ very near the pole,
% presumably due to large argument approximations made in MATLAB.
lambda = [0.9:0.01:1.1] * lambda_p;            
k_0_vec = ones(1,length(lambda))*k_0;
k_vec = ones(1,length(lambda))*k;
u = sqrt(lambda.^2-k_vec.^2);
u_0 = sqrt(lambda.^2-k_0_vec.^2);

Integrand = F(lambda,rho,eps_r,u_0,u,d);
plot(lambda/TMroot,real(Integrand).*(lambda-lambda_p));
grid;
xlabel('\lambda/\lambda_p')
ylabel('Approximate residue')
pause;
print -deps som_res

% Rigorous evaluation of residual.
u = sqrt(lambda_p^2-k^2);
u_0 = sqrt(lambda_p^2-k_0^2);
dD_TE_by_dlambda = lambda_p/u_0 + (lambda_p/u)*coth(u*d) - lambda_p*d*(csch(u*d)^2);
dD_TM_by_dlambda = eps_r*lambda_p/u_0 + (lambda_p/u)*tanh(u*d) + lambda_p*d*(sech(u*d)^2);

Num   = besselj(0,lambda_p*rho) .* lambda_p .* (u_0 + u.*tanh(u*d));
Denom = D_TM(eps_r,u_0,u,d)*dD_TE_by_dlambda + D_TE(u_0,u,d)*dD_TM_by_dlambda;

Residue = Num/Denom

clear u u_o k_0_vec k_vec 
lambda = k_0*[1:0.01:sqrt(eps_r)];
%lambda = k_0*[1:0.01:10];
k_0_vec = ones(1,length(lambda))*k_0;
k_vec = ones(1,length(lambda))*k;
u = sqrt(lambda.^2-k_vec.^2);
u_0 = sqrt(lambda.^2-k_0_vec.^2);
F_full = F(lambda,rho,eps_r,u_0,u,d);
F_sing = Residue*ones(1,length(lambda))./(lambda-lambda_p*ones(1,length(lambda)));

plot(lambda/k_0,real(F_full),...
     lambda/k_0,real(F_full-F_sing),'-.');
%axis([1.0 sqrt(eps_r_prime) -0.6 0.2]);
axis([0 10 -1 +1])
grid;
xlabel('\lambda/k_0')
legend('F','F-F_{sing}',0)
%pause;

clear t F_sing
t = [0:0.01:acosh(sqrt(eps_r_prime))];
k_0_vec = ones(1,length(t))*k_0;
k_vec = ones(1,length(t))*k;
u = sqrt((k_0*cosh(t)).^2-k_vec.^2);
u_0 = sqrt((k_0*cosh(t)).^2-k_0_vec.^2);

F_tx = F(k_0*cosh(t),rho,eps_r,u_0,u,d);
F_sing = Residue*ones(1,length(t))./(k_0*cosh(t) - lambda_p*ones(1,length(t)));
Integrand_t_reg2 = (F_tx - F_sing)*k_0.*sinh(t); % dlambda = k_0*sinh(t) dt

plot(t,abs(F_tx),'-.',t,abs(F_sing),':',t,real(Integrand_t_reg2)/max(abs(Integrand_t_reg2)));
grid;
legend('F','F_Sing','Difference')
axis([0 1.5 -1 3]);
xlabel('t')
%legend('Real part','Imag part',0)
print -deps Som_r2tx
pause;

clear F_full k_0_vec k_vec u u_0
% -------------Region 3: lambda = [k_0*sqrt(eps_r_prime),infty]
lambda = [sqrt(eps_r_prime):0.1:10] * k_0;
k_0_vec = ones(1,length(lambda))*k_0;
k_vec = ones(1,length(lambda))*k;
u = sqrt(lambda.^2-k_vec.^2);
u_0 = sqrt(lambda.^2-k_0_vec.^2);
F_full = F(lambda,rho,eps_r,u_0,u,d);
F_static = besselj(0,lambda*rho)/(1+eps_r);
plot(lambda/k_0,F_full,lambda/k_0,(F_full-F_static),':')
legend('Function','Static term extracted',0)
xlabel('\lambda/k_0')
print -deps Som_r3st
