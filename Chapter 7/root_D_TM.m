function y = root_D_TM(k_0,k,eps_r,d,max_iter)
% Find pole of D_TM using bisection algorithm ("interval chop") 
x2 = k_0*sqrt(real(eps_r)), % Note: approximation of upper limit
u = sqrt(x2^2-k^2);
u_0 = sqrt(x2^2-k_0^2);
fmid = real(D_TM(eps_r,u_0,u,d));
x1 = k_0;
u = sqrt(x1^2-k^2);
u_0 = sqrt(x1^2-k_0^2);
f = real(D_TM(eps_r,u_0,u,d));
if f < 0
   TMroot = x1;
   del_x = x2-x1;
else 
   TMroot = x2;
   del_x = x1-x2;
end 
for kk=1:max_iter,    
   del_x = del_x/2;
   xmid = TMroot+del_x;
   u = sqrt(xmid^2-k^2);
   u_0 = sqrt(xmid^2-k_0^2);
   fmid = real(D_TM(eps_r,u_0,u,d));
   if fmid < 0 
      TMroot = xmid;
   end 
   TMroot/k_0; % Remove semi-colon for screen feedback during root search.
end
y=TMroot;