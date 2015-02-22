function y = F_reg1(t)
% F_reg1: returns the value of the integrand of V at t.
% This is region 1, using the change of variables lambda = k_0*cos(t), 
% i.e. dt = - k_0 sin(t) dlambda. (The minus sign cancels when the 
% integration limits are interchanged). 

global k k_0 rho eps_r d
k_0_vec = k_0*ones(1,length(t));
k_vec   = k*ones(1,length(t));
lambda = k_0*cos(t);
u = sqrt(lambda.^2-k_vec.^2);
u_0 = sqrt(lambda.^2-k_0_vec.^2);
y = F(lambda,rho,eps_r,u_0,u,d).*k_0.*sin(t);
