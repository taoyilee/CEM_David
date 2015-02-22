function y = F_reg2(t)
% F_reg2: returns the value of the integrand of V at t.
% This is region 2, using the change of variables lambda = k_0*cosh(t)
% and subtracting out the singularity.
% Furthermore, dt = k_0 sinh(t) dlambda
global k k_0 rho eps_r eps_r_prime d lambda_p Residue

k_0_vec = k_0*ones(1,length(t));
k_vec   = k*ones(1,length(t));
lambda = k_0*cosh(t);
u = sqrt(lambda.^2-k_vec.^2);
u_0 = sqrt(lambda.^2-k_0_vec.^2);
F_sing = Residue*ones(1,length(t))./(lambda - lambda_p*ones(1,length(t)));
y  = (F(lambda,rho,eps_r,u_0,u,d)-F_sing).*k_0.*sinh(t);

