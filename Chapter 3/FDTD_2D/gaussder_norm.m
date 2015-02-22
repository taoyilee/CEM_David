% A properly normalized Gaussian derivative pulse
function y=gaussder_norm(t,m,sigma)
  %y= -1/sqrt(2*pi)*(t-m)/sigma^3*exp(-(t-m)^2/(2*sigma^2));
  y= -exp(0.5)*(t-m)/sigma*exp(-(t-m)^2/(2*sigma^2));  % Better scaled version
end
