function y = V_int
% Matlab function which, using numerical integration and various special forms, 
% returns the value of the scalar potential V
% as a function of distance for single-layer microstrip, 
% in the regime where only the TM mode propagates. 
% The 1/(2 pi epsilon_0) scaling factor is NOT included!
% D.B.Davidson, Aug 2003

global k k_0 rho eps_r eps_r_prime d lambda_p Residue

% rho
% Integrate in region 1
V_reg1 = quad(@F_reg1,0,pi/2); % Integrate using lambda = k_0 cos(t) change of variables

% check
%delta_t = 0.001
%t = [0:delta_t:pi];
%Int = F_reg1(t);
%V_reg1_check = sum(Int)*delta_t

% Integrate in region 2
% Rigorous evaluation of residual.
% Note that integral must stop just BELOW eps_r*k_0, otherwise D_TE is
% infinite at this point due to the coth function. This is set by the parameter zeta:
zeta = 0.999;
u = sqrt(lambda_p^2-k^2);
u_0 = sqrt(lambda_p^2-k_0^2);
dD_TE_by_dlambda = lambda_p/u_0 + (lambda_p/u)*coth(u*d) - lambda_p*d*(csch(u*d)^2);
dD_TM_by_dlambda = eps_r*lambda_p/u_0 + (lambda_p/u)*tanh(u*d) + lambda_p*d*(sech(u*d)^2);
Num   = besselj(0,lambda_p*rho) .* lambda_p .* (u_0 + u.*tanh(u*d));
Denom = D_TM(eps_r,u_0,u,d)*dD_TE_by_dlambda + D_TE(u_0,u,d)*dD_TM_by_dlambda;
Residue = Num/Denom;

V_reg2_Xsing = quad(@F_reg2,0,acosh(zeta*sqrt(eps_r_prime))); % Integrate using lambda = k_0 cosh(t) change of variables
V_reg2_sing = Residue*log( (k_0*sqrt(eps_r_prime)-lambda_p)/(lambda_p-k_0)) -i*pi*Residue;
V_reg2 = V_reg2_Xsing + V_reg2_sing;

% Integrate in region 3. Again, integrate in this case from just ABOVE eps_r*k_0, again 
% to avoid problems with D_TE. 
lambda_upper = 10*k_0; % This is an estimate only... but is refined iteratively below.
V_reg3_Xstatic1 = quad(@F_reg3,1/zeta*k_0*sqrt(eps_r_prime),lambda_upper);
V_reg3_Xstatic2 = quad(@F_reg3,1/zeta*k_0*sqrt(eps_r_prime),2*lambda_upper);
err =(V_reg3_Xstatic1 - V_reg3_Xstatic2)/abs(V_reg3_Xstatic2);
eps_rel = abs(V_reg2)*1e-4;
count = 0;
max_count = 10; % to stop infinite loops
while (abs(err) > eps_rel & count < max_count) 
   count  = count +1;
   V_reg3_Xstatic1 = V_reg3_Xstatic2;
   lambda_upper = 2*lambda_upper;
   V_reg3_Xstatic2 = quad(@F_reg3,1/zeta*k_0*sqrt(eps_r_prime),2*lambda_upper);
   err =(V_reg3_Xstatic1 - V_reg3_Xstatic2)/abs(V_reg3_Xstatic2);
end 
V_reg3_Xstatic = V_reg3_Xstatic2;
V_reg3 = V_reg3_Xstatic + 1/(rho*(1+eps_r)) ...
        - quad(@F_static,0,sqrt(eps_r_prime)*k_0);
y = V_reg1 + V_reg2 + V_reg3;   
