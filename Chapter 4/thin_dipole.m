% MoM code for thin wires, using collocation and
% piecewise sinusoisdal basis functions. 
% Originally written by Isak Theron as a task. Extended by 
% by David Davidson. 
% Geometry is read in; uses default b for air-filled co-axial
% cable in magnetic frill model.
% Note: best results obtained with 0.48 lambda wire; 0.005 lambda radius, 
% and 50-100 segments. Magnetic frill converges slowly, especially for thinner
% wires.  

clear;
N_seg=input('Enter number of segments (must be even): ')
L=input('Enter wire length (in wavelengths): ')
a=input('Enter wire radius (in wavelengths): ')
b = 3.5*a % for a 75 ohm system.

B = 2*pi; % This is beta, for \lambda=1m.
dz=L/(N_seg); % This is \Delta_z
Bdz=B*dz;
N=N_seg -1 ; % Number of degrees of freedom is one less than segments. 
for n=(1:N)
  Rnp1 = sqrt(a^2+((n-2)*dz)^2);        % Rn+1 
  Rn = sqrt(a^2+((n-1)*dz)^2);          % Rn
  Rnm1 = sqrt(a^2+(n*dz)^2);            % Rn-1
  A = exp(-i*B*Rnp1)/Rnp1;
  BB = exp (-i*B*Rn)*sin(2*Bdz)/(Rn*sin(Bdz));
  C = exp(-i*B*Rnm1)/Rnm1;
  c(1,n) = -i*30*(A-BB+C)/sin(Bdz);
end;
Z = toeplitz(c,c); % This uses the symmetry property of the matrix for 
                   % efficient fill. Alternately, the above loop should be 
                   % replaced by a nested loop in m and n. 
  
% Magnetic frill feed model for dipole
CentreSeg=ceil(N/2);
for m=(CentreSeg:N)
  zz=(m-CentreSeg)*dz;
  R1=sqrt(a^2+zz^2);
  R2=sqrt(b^2+zz^2);
  E(m) = (exp(-i*B*R1)/R1-exp(-i*B*R2)/R2)/(2*log(b/a));
  E(N+1-m)=E(m);
end
V = -conj(E');

CentreSeg = ceil(N/2);
I_mag_frill = Z \ V; 
Zin_mag_frill = 1/(I_mag_frill(CentreSeg))


% Delta gap feed model
deltaZ=dz;
for n = (1:N);
  V(n) = 0;
end;
V(CentreSeg) = -1/deltaZ;
I_delta_gap = Z \ V;
Zin_delta_gap=1/I_delta_gap(CentreSeg) %  this is the dipole impedance
input('Press (almost) any key to continue');

Length=linspace(-L/2+dz/2,L/2-dz/2,N);
plot(Length,abs(I_delta_gap),Length,abs(I_mag_frill),'-.');
%title(['Current on ',num2str(L),'\lambda dipole : ',num2str(N_seg), ' segments']);
xlabel('z/\lambda');
ylabel('|I| [A/m]');
legend('Delta gap', 'Magnetic frill')
pause; 

print -deps thin_dipI

% compute reflection coefficients in 75 Ohm system (magnitude, dB)
gamma_frill = 20*log10(abs((Zin_mag_frill-75)/(Zin_mag_frill+75)))
gamma_delta = 20*log10(abs((Zin_delta_gap-75)/(Zin_delta_gap+75)))



% Obtain the 2 and 1 norm condition numbers (if required)
%TwoNorm = cond(Z)
%DigitsLost = log10(TwoNorm)
%OneNorm = 1/rcond(Z)
