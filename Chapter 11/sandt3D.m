function [S,T] = sandt3D(vertex1,vertex2,vertex3,vertex4)
% SANDT elemental matrix elements for 3D Whitney (CT/LN) element.
%
% This function returns the S and T matrices for a tetrahedral element
% with four vertex coordinates as input, using the formulation
% of Savage & Peterson, "Higher-order vector finite elements for
% tetrahedral cells," IEEE T-MTT, June 1996, pp.874-879.
% See also D.B.Davidson, Chapter 10, "Computational Electromagnetics for RF and Microwave
% Engineering", 2nd edn, presently in preparation.

% Author: D B Davidson, Sept 2009.

x1=vertex1(1);
y1=vertex1(2);
z1=vertex1(3);
x2=vertex2(1);
y2=vertex2(2);
z2=vertex2(3);
x3=vertex3(1);
y3=vertex3(2);
z3=vertex3(3);
x4=vertex4(1);
y4=vertex4(2);
z4=vertex4(3);

global LOCALEDGENODES

volume = 1/6*abs(det([1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]));
temp = inv([x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4; 1 1 1 1]);
b = temp(:,1);
c = temp(:,2);
d = temp(:,3);
a = temp(:,4);

nabla_lambda=zeros(4,3);
for ii=1:4 % node loop
    nabla_lambda(ii,:) = [b(ii) c(ii) d(ii)];
    for jj = 1:4 % node loop
        phi(ii,jj) = b(ii)*b(jj)+c(ii)*c(jj)+d(ii)*d(jj) ;
        nabla_lambda(jj,:) = [b(jj) c(jj) d(jj)];
        % keyboard
        v(ii,jj,:) = cross(nabla_lambda(ii,:),nabla_lambda(jj,:));
    end
end

M = 1/20*[2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2];

% Compute S and T
for ii=1:6 % edge loop
    for jj = 1:6 % edge loop
        i1 = LOCALEDGENODES(ii,1);
        i2 = LOCALEDGENODES(ii,2);
        j1 = LOCALEDGENODES(jj,1);
        j2 = LOCALEDGENODES(jj,2);
        S(ii,jj) = 4*volume*dot(v(i1,i2,:),v(j1,j2,:));
        T(ii,jj) = volume*( phi(i2,j2)*M(i1,j1) - phi(i2,j1)*M(i1,j2) - phi(i1,j2)*M(i2,j1) + phi(i1,j1)*M(i2,j2));
    end
end
