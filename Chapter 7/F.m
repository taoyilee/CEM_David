function y = F(lambda,rho,eps_r,u_0,u,d)
% F: returns the value of the integrand of scalar potential V at length(lambda) points.
% eps_r is the (possibly complex) relative dielectric permittivity. 
% u_0 and u are the Sommerfeld parameters in free space and the dielectric
% respectively. (u_0 and u may be vectors). d is substrate thickness [m]. 
% The 1/(2 pi epsilon_0) scaling factor is NOT included. 
y = besselj(0,lambda*rho) .* lambda .* (u_0 + u.*tanh(u*d)) ./ ...
            (D_TE(u_0,u,d).*D_TM(eps_r,u_0,u,d));
