function [S] = s_nodal(x1,y1,x2,y2,x3,y3)
% S_nodal elemental matrix elements for a first order nodal element.
% This function returns the S matrix for a triangular element
% with vertex coordinates (x1,y1),(x2,y2) and (x3,y3), using a formulation 
% described in D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", Chapter 10 of 2nd edn, presently in  preparation. 

% Author: D B Davidson, Aug 20009

area = 0.5*abs(det([1 x1 y1; 1 x2 y2; 1 x3 y3]));
temp = inv([x1 x2 x3; y1 y2 y3; 1 1 1]);
b = temp(:,1);
c = temp(:,2);
a = temp(:,3);

for ii=1:3,
   for jj = 1:3,
       phi(ii,jj) = b(ii)*b(jj)+c(ii)*c(jj) ;
   end 
end


% Compute S 
for ii=1:3,
   for jj = 1:3,
     S(ii,jj) = area*phi(ii,jj);
   end 
end
