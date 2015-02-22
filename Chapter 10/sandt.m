function [S,T] = sandt(x1,y1,x2,y2,x3,y3)
% SANDT elemental matrix elements for 2D Whitney (CT/LN) element.
% This function returns the S and T matrices for a triangular element
% with vertex coordinates (x1,y1),(x2,y2) and (x3,y3), using a formulation 
% based on Savage & Peterson, "Higher-order vector finite elements for 
% tetrahedral cells," IEEE T-MTT, June 1996, pp.874-879.  

% Author: D B Davidson, October 2003.

global LOCALEDGENODES

area = 0.5*abs(det([1 x1 y1; 1 x2 y2; 1 x3 y3]));
temp = inv([x1 x2 x3; y1 y2 y3; 1 1 1]);
b = temp(:,1);
c = temp(:,2);
a = temp(:,3);

for ii=1:3,
   for jj = 1:3,
       phi(ii,jj) = b(ii)*b(jj)+c(ii)*c(jj) ;
       v_z(ii,jj) = b(ii)*c(jj)-b(jj)*c(ii) ;
   end 
end

M = 1/12*[2 1 1; 1 2 1; 1 1 2];

% Compute S and T
for ii=1:3,
   for jj = 1:3,
     i1 = LOCALEDGENODES(ii,1);
     i2 = LOCALEDGENODES(ii,2);
     j1 = LOCALEDGENODES(jj,1);
     j2 = LOCALEDGENODES(jj,2);
     S(ii,jj) = 4*area*v_z(i1,i2)*v_z(j1,j2);
     T(ii,jj) = area*( phi(i2,j2)*M(i1,j1) - phi(i2,j1)*M(i1,j2) - phi(i1,j2)*M(i2,j1) + phi(i1,j1)*M(i2,j2));
   end 
end
