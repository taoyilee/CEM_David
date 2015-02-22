function y = cyl_TM_echo_width(a,k,N,phi)
% Function to compute echo width of a PEC circular cylinder, TM_z
% incidence, using analytical result
% Inputs: a; radius of cylinder [m]
%         k; wavenumber [rad/m] (May be a vector)
%         phi; angle of scattering in radians, = pi for monostatic
% Returns the echo width [m] (NB! This is the absolute RCS).
% Method: Implements eqn. (11-102), p. 607, CA Balanis, Advanced EM Engineering, Wiley 89, .
% Written by DB Davidson, 25 Feb 2008.


lambda = 2*pi./k;
ka = k*a;
echo = zeros(size(k));
for n=0:N,
    [tempH,ierr] =  besselj(n,ka);
    for kk = 1:length(k)
        check_err(ierr(kk))
    end
    [tempH,ierr] =  besselh(n,2,ka);
    % repeat on Hankel function
    for kk = 1:length(k)
        check_err(ierr(kk))
    end
    if n==0
        epsilon_n=1;
    else
        epsilon_n=2;
    end
    echo = echo + epsilon_n*besselj(n,ka)./besselh(n,2,ka)*cos(n*phi);
end
y= 2*lambda/pi.*abs(echo).^2;

function check_err(ierr) % Do some error checking on the function evaluation
switch ierr
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
