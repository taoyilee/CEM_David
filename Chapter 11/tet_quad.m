function [w,lambda] = tet_quad (n)
% This function returns the weights and evaluation points (in simplex
% coordinates lambda) from P.Keast, "Moderate degree tetrahedral quadrature formulaes",
% Comput Metho Appl Mech Eng, Vol 55, 1986, pp. 339-348.
% Entered by DB Davidson, Oct 2009
switch n
    case 11 % symmetric 11 point rule, degree of precision 4
        zz = 1.0/4.0;
        aa = 0.714285714285714285e-01;
        bb = 0.785714285714285714;
        cc = 0.399403576166799219;
        dd = 0.100596423833200785;
        w1 = -0.1315555555555555550e-01;
        w2 = 0.7622222222222222222e-02;
        w3 = 0.2488888888888888880e-01;
        lambda = [zz zz zz zz;
            aa aa aa bb;
            aa aa bb aa;
            aa bb aa aa;
            bb aa aa aa;
            cc cc dd dd;
            cc dd cc dd;
            cc dd dd cc;
            dd dd cc cc;
            dd cc dd cc;
            dd cc cc dd];
        w = 6*[w1 w2 w2 w2 w2 w3 w3 w3 w3 w3 w3]'; 
        % Keast's rules include a factor of 1/6 for unitary tet, which is
        % handled elsewhere in this code (by multiplying the integrand by
        % the volume).

    otherwise
        error 'Unimplemented rule'
end

