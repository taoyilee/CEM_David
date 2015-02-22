function y = D_TE(u_0,u,d)
% D_TE: returns the value of the TE surface wave.
% u_0 and u are the Sommerfeld parameters in free space and the dielectric
% respectively. (u_0 and u may be vectors). d is substrate thickness [m]. 
y = u_0+u.*coth(u*d);
