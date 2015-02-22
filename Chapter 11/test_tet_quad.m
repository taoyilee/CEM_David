function test_tet_quad(n)
% TEST_TET_QUAD tests the n-pt quadrature rules for a variety of 
% integrals up to 4th order.
% For analytical results, see Appendix, Useful formulas for simplex
% coordinates.

[w,lambda]=tet_quad(n)
% Test the integral of lambda1^2*lambda1^2:
quad_int=0;
for ii=1:n
    quad_int = quad_int+w(ii)*lambda(ii,1)^2*lambda(ii,2)^2; % times V
end
int=factorial(3)*factorial(2)^2/factorial(7); % times V
err=abs(int-quad_int)/quad_int % V's cancel
% Test the integral of lambda1*lambda2*lambda3*lambda4:
quad_int=0;
for ii=1:n
    quad_int = quad_int+w(ii)*lambda(ii,1)*lambda(ii,2)*lambda(ii,3)*lambda(ii,4);
end
int=factorial(3)/factorial(7);
err=abs(int-quad_int)/quad_int

% Test the integral of a constant 
quad_int=0;
for ii=1:n
    quad_int = quad_int+w(ii);
end
int=factorial(3)/factorial(3);
err=abs(int-quad_int)/quad_int

    