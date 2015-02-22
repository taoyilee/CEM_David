function y = F_reg3(lambda)
% F_reg3: returns the value of the integrand of V at lambda.
% This is region 3, with the static term extracted.
% No change of variables is done in this region. 

global k k_0 rho eps_r d 

k_0_vec = k_0*ones(1,length(lambda));
k_vec   = k*ones(1,length(lambda));
u = sqrt(lambda.^2-k_vec.^2);
u_0 = sqrt(lambda.^2-k_0_vec.^2);
F_static = besselj(0,lambda*rho)/(1+eps_r);
y = (F(lambda,rho,eps_r,u_0,u,d)-F_static);
