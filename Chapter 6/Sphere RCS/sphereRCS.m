function y = sphereRCS(a,k,N)
% Function to compute monostatic RCS of a PEC sphere.
% Inputs: a; radius of sphere [m]
%         k; wavenumber [rad/m] (May be a vector)
% Returns the RCS [m^2] (NB! This is the absolute RCS). 
% Method: Implements eqn. (11-247), p. 657, CA Balanis, Advanced EM Engineering, Wiley 89, .
% See Chapter 6, Section 6.3.1, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP 2005.
% Written by DB Davidson, 10 Nov 2003. Revised 2 Feb 2005.

lambda = 2*pi./k;
ka = k*a;
RCS = zeros(size(k));
for n=1:N,
  [tempH,ierr] =  besselh(n+1/2,2,ka);
  % Do some error checking on the function evaluation
  for kk = 1:length(k)
      switch ierr(kk) 
          case(1)
              disp('Illegal arguments.');pause
          case(2)
              disp('Overflow.  Return Inf.');pause
          case(3)
              disp('Some loss of accuracy in argument reduction.');pause
          case(4)
              disp('Complete loss of accuracy, z or nu too large.');pause
          case(5)
              disp('No convergence.  Return NaN.');pause
          otherwise
              % No problems with function evaluation encountered, continue.
      end
  end
  Hcaret2n = sqrt(pi*ka/2) .* besselh(n+1/2,2,ka);
  Hcaret2nprime = 1/2*sqrt(pi./(2*ka)).*besselh(n+1/2,2,ka)+...
                  sqrt(pi*ka/2) .* (-besselh(n+3/2,2,ka) + ...
                  (n+1/2)*besselh(n+1/2,2,ka)./ka);
  RCS = RCS + (-1)^n * (2*n+1)./(Hcaret2n.*Hcaret2nprime); 
end
y= lambda.^2/(4*pi).*abs(RCS).^2;

