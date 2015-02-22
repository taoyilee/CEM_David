% Determine dominant eigenfrequency for a rectangular cavity using the FDTD method. 
% Usage: fdtd_3d_demo
% The code assumes a uniform free-space filling, and that the dominant mode
% is TE101, for a 1x0.5x0.75 PEC cavity (although other cavities could be
% analyzed with minor changes).
% See D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", Chapter 3 of 2nd edn, presently in preparation.


% Written by D B Davidson, 9 Nov 2009.

clear all;
close all;

eps_0 = 8.854e-12;
mu_0=4*pi*1e-7;
eps_r = 1.0;
mu_r = 1.0;

L_x=1; % x-dimensions of resonant cavity in m
L_y=0.5; % ditto y
L_z=0.75; % ditto z

N_x = 8; % Number of Yee cells in x, y and z directions.
N_y = 4;
N_z = 6;
refine=input('Factor to refine mesh by? (default 2): ');
if isempty(refine)
    refine = 2
end
N_x = N_x*refine;
N_y = N_y*refine;
N_z = N_z*refine;


Delta_x = L_x/N_x; % Spatial cell size in x, y and z directions
Delta_y = L_y/N_y; 
Delta_z = L_z/N_z; 

Delta_s = min([Delta_x Delta_y Delta_z]) % Courant limit given by smallest spatial step
c= 1/sqrt(eps_0*mu_0); 
CFL = Delta_s/(sqrt(3)*c); % Courant limit

N_t = 2^10*refine; % power of 2
Delta_t = CFL*1; % 


% Dimension the field vectors - will be overwritten in place
E_x = zeros(N_x,  N_y+1,N_z+1);
E_y = zeros(N_x+1,N_y,  N_z+1);
E_z = zeros(N_x+1,N_y+1,N_z);
H_x = zeros(N_x+1,N_y,  N_z);
H_y = zeros(N_x,  N_y+1,N_z);
H_z = zeros(N_x,  N_y,  N_z+1);

mu = mu_0*mu_r; % In a general code, this could vary from cell to cell
epsilon = eps_0*eps_r; % Ditto.

% Inject source - random noise with zero mean
% Only TE_z modes sought - only H_z excited.

H_z(:,:,:) = rand(N_x,N_y,N_z+1)-0.5;

% Allocate storage for field samples
t = zeros(N_t,1);
t(1) = Delta_t;
H_z_t = zeros(N_t,1);

% Time step the fields 
for n=1:N_t
    n
    % Save Hz field in centre of box
    t(n) = (n)*Delta_t;
    H_z_t(n) = H_z(round(N_x/2),round(N_y/2),round(N_z/2));
    % Update H fields for t^n+1/2. Dimensions are conforming.
    H_x = H_x + Delta_t/mu * (diff(E_y,1,3)/Delta_z - diff(E_z,1,2)/Delta_y);        
    H_y = H_y + Delta_t/mu * (diff(E_z,1,1)/Delta_x - diff(E_x,1,3)/Delta_z);        
    H_z = H_z + Delta_t/mu * (diff(E_x,1,2)/Delta_y - diff(E_y,1,1)/Delta_x);        
    % Update E fields for t^n+1. Dimensions require explicit specifiation to make conforming. 
    E_x(:,2:N_y,2:N_z) = E_x(:,2:N_y,2:N_z) + Delta_t/epsilon * (diff(H_z(:,:,2:N_z),1,2)/Delta_y - diff(H_y(:,2:N_y,:),1,3)/Delta_z);
    E_y(2:N_x,:,2:N_z) = E_y(2:N_x,:,2:N_z) + Delta_t/epsilon * (diff(H_x(2:N_x,:,:),1,3)/Delta_z - diff(H_z(:,:,2:N_z),1,1)/Delta_x);
    E_z(2:N_x,2:N_y,:) = E_z(2:N_x,2:N_y,:) + Delta_t/epsilon * (diff(H_y(:,2:N_y,:),1,1)/Delta_x - diff(H_x(2:N_x,:,:),1,2)/Delta_y);   
end

Delta_f = 1/(N_t*Delta_t);
for kk = 1:N_t/2
    freq(kk) = (kk-1)*Delta_f;
end
% Find first few analytical results:
k_exact = [5.2359878 7.55154489 8.1788743 8.9472598]; % rad/m
f_exact = k_exact * c/(2*pi); % Hz
H_z_f = fft(H_z_t);
%scale = max(abs(H_z_f(1:N_t/2));
%dummy = scale*[1 1 1];
dummy = [0 0 0 0];
plot(freq,abs(H_z_f(1:N_t/2)),f_exact,dummy,'x');
pause;

% Improve the plot - only plot out to 500 MHz. Also, don't include DC term;
% this can be quite significant, depending on initial conditions.
figure
end_pt = floor(5e8/Delta_f)
norm = max(abs(H_z_f(2:end_pt)))
plot(freq(2:end_pt)/1e6,abs(H_z_f(2:end_pt))/norm) % ,'k-',f_exact,dummy,'x');
f_exact = f_exact/1e6;
line([f_exact(1) f_exact(1)],[0 1],'LineStyle','--','Color','k')
line([f_exact(2) f_exact(2)],[0 1],'LineStyle','--','Color','k')
line([f_exact(3) f_exact(3)],[0 1],'LineStyle','--','Color','k')
line([f_exact(4) f_exact(4)],[0 1],'LineStyle','--','Color','k')
xlabel('f (MHz)')
ylabel('H_z (normalized)')




