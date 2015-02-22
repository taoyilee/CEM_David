function y = F_static(lambda)
% F_reg3_static: returns the value of the integrand of V at lambda.
% This is the static term.
% No change of variables is done. 
global rho eps_r
y = bessel(0,lambda*rho)/(1+eps_r);
