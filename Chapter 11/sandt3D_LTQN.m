function [S,T] = sandt3D_LTQN(vertex1,vertex2,vertex3,vertex4)
% SANDT elemental matrix elements for 3D LT/QN element.
%
% This function returns the S and T matrices for a tetrahedral element
% using LT/QN (mixed 2nd order elements)
% with four vertex coordinates as input, using the formulation
% of Savage & Peterson, "Higher-order vector finite elements for
% tetrahedral cells," IEEE T-MTT, June 1996, pp.874-879.
% See also D.B.Davidson, Chapter 10, "Computational Electromagnetics for RF and Microwave
% Engineering", 2nd edn, presently in preparation.

% The matrix elements are packed as follows (e1 are the CTLN functions, e2 the LTLN functions, 
% and f1 and f2 the face-based LTQN functions)

% Se1e1 Se1e2 Se1f1 Se1f2
% Se2e1 Se2e2 Se2f1 Se2f2
% Sf1e1 Sf1e2 Sf1f1 Sf1f2
% Sf2e1 Sf2e2 Sf2f1 Sf2f2

% Author: D B Davidson, Oct 2009.

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

volume = 1/6*abs(det([1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]));
temp = inv([x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4; 1 1 1 1]);
b = temp(:,1);
c = temp(:,2);
d = temp(:,3);
a = temp(:,4);

nabla_lambda=zeros(4,3);
for ii=1:4 % node loop
    nabla_lambda(ii,:) = [b(ii) c(ii) d(ii)];
%     for jj = 1:4 % node loop
%         phi(ii,jj) = b(ii)*b(jj)+c(ii)*c(jj)+d(ii)*d(jj) ;
%         nabla_lambda(jj,:) = [b(jj) c(jj) d(jj)];
%         % keyboard
%         v(ii,jj,:) = cross(nabla_lambda(ii,:),nabla_lambda(jj,:));
%     end
end

S=zeros(20,20);
T=zeros(20,20);
quad_pts=11;
[w,lambda] = tet_quad(11);
for nn=1:quad_pts
    %lambda(nn,:)
    %nabla_lambda
    funcs = LTQN3D(lambda(nn,:),nabla_lambda);
    curl_funcs = curl_LTQN3D(lambda(nn,:),nabla_lambda);
    for ii=1:20
        for jj=1:20
            S(ii,jj) = S(ii,jj)+w(nn)*dot(curl_funcs(ii,:),curl_funcs(jj,:));
            T(ii,jj) = T(ii,jj)+w(nn)*dot(funcs(ii,:),funcs(jj,:));
        end
    end
end
S=S*volume;
T=T*volume;
