function [w,lambda] = tri_quad (n)
% This function returns the weights and evaluation points (in simplex
% coordinates lambda) from D.A.Dunavant, "Gaussian quadrature formulas for
% triangles", IJNME, vol 21,  1985, pp. 1129-1148.

switch n
    case 6 % symmetric 6 point rule, degree of precision 4
        w = zeros(6,1);
        alpha = zeros(6,1);
        beta = zeros(6,1);
        gamma = zeros(6,1);
        w(1) =     0.223381589678011;
        w(2) = w(1);
        w(3) = w(1);
        w(4) =     0.109951743655322;
        w(5) = w(4);
        w(6) = w(4);
        alpha(1) = 0.108103018168070;
        beta(1)  = 0.445948490915965;
        gamma(1) = beta(1);
        alpha(2) = beta(1);
        alpha(3) = beta(1);
        beta(2)  = alpha(1);
        beta(3)  = beta(1);
        gamma(2) = beta(1);
        gamma(3) = alpha(1);
        alpha(4) = 0.816847572980459;
        beta(4)  = 0.091576213509771;
        gamma(4) = beta(4);
        alpha(5) = beta(4);
        alpha(6) = beta(4);
        beta(5)  = alpha(4);
        beta(6)  = beta(4);
        gamma(5) = beta(4);
        gamma(6) = alpha(4);
    otherwise
        error 'Unimplemented rule'
end
lambda(1:n,1) = alpha;
lambda(1:n,2) = beta;
lambda(1:n,3) = gamma;
