function [S,T] = sandt_LTQN(x1,y1,x2,y2,x3,y3)
% SANDT elemental matrix elements for 2D LT/QN element.
% This function returns the S and T matrices for a triangular element
% with vertex coordinates (x1,y1),(x2,y2) and (x3,y3), using a formulation
% based on Savage & Peterson, "Higher-order vector finite elements for
% tetrahedral cells," IEEE T-MTT, June 1996, pp.874-879, and extensions for
% higher order, see D.B.Davidson, Chapter 10, "Computational Electromagnetics for RF and Microwave
% Engineering", 2nd edn, presently in preparation.

% The matrix elements are packed as follows (e1 are the CTLN functions, e2 the LTLN functions, and f the LTQN functions)

% Se1e1 Se1e2 Se1f
% Se2e1 Se2e2 Se2f
% Sfe1  Sfe2  Sff

% Author: D B Davidson, Aug 2009.

global LOCALEDGENODES

area = 0.5*abs(det([1 x1 y1; 1 x2 y2; 1 x3 y3]));
temp = inv([x1 x2 x3; y1 y2 y3; 1 1 1]);
b = temp(:,1);
c = temp(:,2);
a = temp(:,3);

nabla_lambda(1,:)=[b(1) c(1)];
nabla_lambda(2,:)=[b(2) c(2)];
nabla_lambda(3,:)=[b(3) c(3)];

S=zeros(8,8);
T=zeros(8,8);
quad_pts=6;
[w,lambda] = tri_quad(6);
for nn=1:quad_pts
    %lambda(nn,:)
    %nabla_lambda
    funcs = LTQN(lambda(nn,:),nabla_lambda);
    curl_funcs = curl_LTQN(lambda(nn,:),nabla_lambda);
    for ii=1:8
        for jj=1:8
            S(ii,jj) = S(ii,jj)+w(nn)*dot(curl_funcs(ii,:),curl_funcs(jj,:));
            T(ii,jj) = T(ii,jj)+w(nn)*dot(funcs(ii,:),funcs(jj,:));
        end
    end
end
S=S*area;
T=T*area;