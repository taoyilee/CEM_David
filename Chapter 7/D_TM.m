function y = D_TM(eps_r,u_0,u,d)
% D_TM: returns the value of the TM surface wave.
% eps_r is the (possibly complex) relative dielectric permittivity. 
% u_0 and u are the Sommerfeld parameters in free space and the dielectric
% respectively. (u_0 and u may be vectors). d is substrate thickness [m]. 
y = eps_r*u_0+u.*tanh(u*d);
