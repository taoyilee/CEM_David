% A MATLAB implementation of the 2D FDTD theory for scattering off a PEC cylinder,
% as described in Chapter 3, Section 3.2, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave Engineering", CUP 2005.
% Usage: fdtd_2D_demo_v2

% Options are provided at run-time to include the  PEC cylinder or not; and
% also to refine the mesh in a simple fashion. Enter 1 for reasonable
% results, with reasonable execution time. 

% The first plot shows the position of the cylinder, and the
% scattered/total field boundary. Subsequent plots show the pulse
% propagating through the mesh, and scattering off the cylinder. This is
% repeated as a movie.

% Also shown is a plot in the time domain of the scattered field at coordinates (point1_x,point1_y).
% Data thus saved can be used to produce results similar to Fig. 3.13, using the
% "cyl_plot_v2.m" script.

% Commented out sections of code indicate how systematic testing may be
% undertaken, as discussed in Section 3.2.6.

% Author: D.B.Davidson
% Originally written 22 Feb 2003, revised 2 Feb 2005.
% Revised again 7 June 2005, to include the scattered/total field
% formulation on ALL sides (not just the on the left side of Fig 3.1, p.
% 72). See also additional text material.
% Revised again 09 Oct 2007, to correct dimensions of E and H fields, DBD.
% Correction also in time offset terms. 
% Minor revisions 08 Mar 2010. DBD.

% Known issues: 
% 1. A small signal is leaking into the simulation region from the 
% the corners of the front total/scattered field interface. It should be
% gated out.
% 2. In the R2009a release of MATLAB, the movie function is not working
% correctly.


clear;
clf;
format compact;
% Selectively includes PEC cylinder, centered at (N_centre_x,N_centre_y)
cyl_present=input('Include PEC cylinder? (1 or true)')

% Physical constants

c = 2.997925e8;    % [m/s] Speed of light in vacuum
eps_0 = 8.854e-12; % [F/m] epsilon_0
mu_0 = 4*pi*1e-7;  % [H/m] mu_0
eta_0 = sqrt(mu_0/eps_0);
                   % [Ohm] Wave impedance of free space
% Set parameters for simulation

refine = input('Factor to refine mesh? 1 standard')    % A quick method to refine the mesh 
pulse_compress = 1 %2 % A quick way to shorten the pulse. For Figs 3.11-14, use 2. For 3.15, use 1.

N_x = refine*200 % refine*400     % number of cells in x-direction
N_y = refine*100 % refine*200            % ditto y. 

M = refine*512   % refine*1024      % Number of time steps
L = round(N_x/2)    % scat/tot field boundary on left side of Fig. 3.1 as Section 3.2
L1 = 20  % 10              % New scat/tot field boundaries on upper, lower and right side of Fig. 3.1
%delta_s = 0.005/refine % [m] spatial step
delta_s = 0.005*2/refine % [m] spatial step. This generates a physically larger computational volume.
R = 1        % fraction of Courant limit. Must be <= 1
delta_t = R* delta_s/(c * sqrt(2)) 
              % [s] Time step size
sigma = 1.0e-10/pulse_compress
              % Controls spectral content of Gaussian derivative pulse - 
              % equals 1/omega_max
f_max = (1/sigma)/(2*pi) % Freq of largest significant spectral component, in Hz. Information only
m_offset = 4*sigma;  % Controls switch-on time
Peak = 1;      % Peak amplitude of E field
% Set parameters for PEC cylinder
radius = 0.03 % [m] radius of cylinder
N_centre_x = round(0.75*N_x)
N_centre_y = round(0.5*N_y)

% Check that the simulation specification is valid:

if ((N_centre_x-L)*delta_s <= radius)
   disp('Error in simulation data. Scattered/total field not entirely to the left of target')
   return
end


%------------------------------------------

% Set up material grid (free space to start)
C_Ex = ones(N_x+1,N_y+1)*delta_t/(eps_0*delta_s); 
C_Ey = ones(N_x+1,N_y+1)*delta_t/(eps_0*delta_s); 
D_Hz = ones(N_x,N_y)*delta_t/(mu_0*delta_s); 
% Now force the electric fields to zero inside (and on the surface of) the PEC
% Note that the indices of the centre are treated as per usual FDTD
% indices, i.e. the actual location is:
%x_c=(N_centre_x-1))*delta_s ; y_c=(N_centre_x-1))*delta_s
if cyl_present % Otherwise just leave it as free space
  for ii = 1:N_x+1
    for jj = 1:N_y+1
        if ( sqrt( ((ii-1/2-(N_centre_x-1))*delta_s)^2 +  ((jj-1-(N_centre_y-1))*delta_s)^2 ) <=  radius ) 
           C_Ex(ii,jj) = 0;
        end
        if ( sqrt( ((ii-1-(N_centre_x-1))*delta_s)^2 +  ((jj-1/2-(N_centre_y-1))*delta_s)^2 ) <=  radius ) 
           C_Ey(ii,jj) = 0;
        end
     end
  end
end 
%------------------------------------------

% Set up storage for time histories.
H_z_point1 = zeros(1,M);
E_y_point1 = zeros(1,M);
point1_x = N_x/4;
% This is another hack! also remove!!
%point1_x = 100;
point1_y = N_y/2;
H_z_point2 = zeros(1,M);
E_y_point2 = zeros(1,M);
point2_x = round((N_x+L)/2);
point2_y = N_y/2;


% Produce a simple graphical output, showing the cylinder, scat/tot zone
% interface and the point at which the scattered field will be computed. 
mesh_pic=zeros(N_x+1,N_y+1);
for ii = 1:N_x+1
    for jj = 1:N_y+1
      if C_Ex(ii,jj) == 0
        mesh_pic(ii,jj) = N_x/2; % To get vertical scale in plot OK when plotted
      elseif ii==L
        mesh_pic(ii,jj) = N_x/4;
      elseif (ii==point1_x & jj==point1_y)
        mesh_pic(ii,jj) = N_x/2;
      end
    end 
end
mesh((1:1:N_y+1),(1:1:N_x+1),mesh_pic);
axis image;
title('Simulation region')
disp('Press any key to continue')
pause;


% First time step - Initialize values for H_z, E_x and E_y
H_z_nmin1 = zeros(N_x,N_y); %
E_x_nmin1 = zeros(N_x+1,N_y+1); %
E_y_nmin1 = zeros(N_x+1,N_y+1); %
% Pre-allocation
H_z_n = zeros(N_x,N_y); %
E_x_n = zeros(N_x+1,N_y+1); %
E_y_n = zeros(N_x+1,N_y+1); %
% Get CPU time
start_time=cputime;
%------------------------------------------

movie_count = 1;
movie_interval = refine*50 % ,refine*10 % 

report_time_interval = refine*50; 

% Time loop
for m = 2:M,

  % A note on dimensions. Due to the nature of the Yee cell, the E and H fields 
  % have slightly different dimensions. 
  % E_x: N_x+1 by N_y+1
  % E_y: N_x+1 by N_y+1
  % H_z: N_x   by N_y
  % Essentially, it is because the H field are in the interior of the cells, whereas the E fields are
  % on the exterior of the cells. 
  % This must be explicitly taken into account in the matrix sizes and the update equations.
    
  % ---------------------------- H field update -----------------------------------------------------------------                  
  H_z_n(1:N_x,1:N_y) = H_z_nmin1(1:N_x,1:N_y) ...
      + D_Hz(1:N_x,1:N_y).*(  E_x_nmin1(1:N_x,2:N_y+1)   - E_x_nmin1(1:N_x,1:N_y) ...
                           + E_y_nmin1(1:N_x,1:N_y) - E_y_nmin1(2:N_x+1,1:N_y) ) ;
  % Drive a test line source - used to check basic operation
  % H_z_n(N_x/2,N_y/2) = gaussder((m-1)*delta_t,m_offset,sigma);

  % Special update on scat/tot field boundary
  E_y_nmin1_inc_front = ones(1,N_y+1)*Peak*gaussder_norm((m-1)*delta_t - (L-1)*delta_s/c,m_offset,sigma) ;
  % The H_z field is the total field.
  H_z_n(L,L1:N_y-L1) = H_z_nmin1(L,L1:N_y-L1) ...
      + D_Hz(L,L1:N_y-L1).*(  E_x_nmin1(L,L1+1:N_y-L1+1)     - E_x_nmin1(L,L1:N_y-L1) ...
                        + E_y_nmin1(L,L1:N_y-L1) + E_y_nmin1_inc_front(L1:N_y-L1) - E_y_nmin1(L+1,L1:N_y-L1)) ;
  % Special update on additional new scat/tot field boundary (only needed for Ey) 
  % on the right hand side of Fig. 3.1, at N_x - L1
  % Note now the the "far" side of the ABC is now the SCATTERED, not TOTAL,
  % field, so the role of the E_y fields swops around.
  E_y_nmin1_inc_back = ones(1,N_y+1)*Peak*gaussder_norm((m-1)*delta_t - (N_x-L1+1-1)*delta_s/c,m_offset,sigma) ;
  % E_y_nmin1_inc can be overwritten since it is not used again.   Again,
  % the H_z field is the total field.
   H_z_n(N_x-L1,L1:N_y-L1) = H_z_nmin1(N_x-L1,L1:N_y-L1) ...
       + D_Hz(N_x-L1,L1:N_y-L1).*(  E_x_nmin1(N_x-L1,L1+1:N_y-L1+1)     - E_x_nmin1(N_x-L1,L1:N_y-L1) ...
                         + E_y_nmin1(N_x-L1,L1:N_y-L1) - E_y_nmin1_inc_back(L1:N_y-L1) - E_y_nmin1(N_x-L1+1,L1:N_y-L1)) ;
  % Special update on additional new scat/tot field boundary (only needed for Ey) 
  % on the upper side of Fig. 3.1, at N_y - L1. This plays the same role as
  % the ABC at N_x-L1. 
  for ii=1:N_x+1,
    E_x_nmin1_inc_top(ii,1) = 0 ; % For this specific polarization (Ey,Hz) the incident x-field is zero
  end;
  % The H_z field that follows is the total field. The E_x field above the interface is the scattered field, that 
  % below, the total field.
  % Note: since the incident field is zero, the following code stub has not
  % thus been properly tested. 
  

  H_z_n(L:N_x-L1,N_y-L1) = H_z_nmin1(L:N_x-L1,N_y-L1) ...
     + D_Hz(L:N_x-L1,N_y-L1).*(  E_x_nmin1(L:N_x-L1,N_y-L1+1) + E_x_nmin1_inc_top(L:N_x-L1,1) - E_x_nmin1(L:N_x-L1,N_y-L1) ...
                       + E_y_nmin1(L:N_x-L1,N_y-L1) - E_y_nmin1(L+1:N_x-L1+1,N_y-L1)) ;                  
  % Ditto at scat/tot field boundary (only needed for Ey) 
  % on the lower side of Fig. 3.1, at L1
  for ii=1:N_x+1,
    E_x_nmin1_inc_bottom(ii,1) = 0 ; % Again, for this specific polarization (Ey,Hz) the incident x-field is zero
  end;
  % Note: same caution as above. Also note that role of E_x fields swops
  % around - E_x field above interface is now total field.  

  H_z_n(L:N_x-L1,L1) = H_z_nmin1(L:N_x-L1,L1) ...
     + D_Hz(L:N_x-L1,L1).*(  E_x_nmin1(L:N_x-L1,L1+1)  - E_x_nmin1(L:N_x-L1,L1) - E_x_nmin1_inc_bottom(L:N_x-L1,1) ...
                       + E_y_nmin1(L:N_x-L1,L1) - E_y_nmin1(L+1:N_x-L1+1,L1)) ;                  

                                   
  % Very special update - 4 single points on corners of scat/tot field boundary
  % The H_z field is the total field.
  % bottom front corner
  H_z_n(L,L1) = H_z_nmin1(L,L1) ...
      + D_Hz(L,L1).*(  E_x_nmin1(L,L1+1) -  E_x_nmin1_inc_bottom(L)  - E_x_nmin1(L,L1) ...
                        + E_y_nmin1(L,L1) + E_y_nmin1_inc_front(L1) - E_y_nmin1(L+1,L1)) ;
  % top front corner
  H_z_n(L,N_y-L1) = H_z_nmin1(L,N_y-L1) ...
      + D_Hz(L,N_y-L1).*(  E_x_nmin1(L,N_y-L1+1) + E_x_nmin1_inc_top(L)  - E_x_nmin1(L,N_y-L1) ...
                        + E_y_nmin1(L,N_y-L1) + E_y_nmin1_inc_front(N_y-L1) - E_y_nmin1(L+1,N_y-L1)) ;

  % At present, following two stubs are commented out - they inject a spurious signal at these corners. It is not clear why.  
%                     
%   % back bottom corner
%   H_z_n(N_x-L1,L1) = H_z_nmin1(N_x-L1,L1) ...
%       + D_Hz(N_x-L1,L1).*(  E_x_nmin1(N_x-L1,L1+1) - E_x_nmin1_inc_bottom(N_x-L1)  - E_x_nmin1(N_x-L1,L1) ...
%                         + E_y_nmin1(N_x-L1,L1) - E_y_nmin1_inc_back(L1) - E_y_nmin1(N_x-L+1,L1)) ;
%   % back top corner
%   H_z_n(N_x-L1,N_y-L1) = H_z_nmin1(N_x-L1,N_y-L1) ...
%       + D_Hz(N_x-L1,N_y-L1).*(  E_x_nmin1(N_x-L1,N_y-L1+1) + E_x_nmin1_inc_top(N_x-L1)  - E_x_nmin1(N_x-L1,N_y-L1) ...
%                         + E_y_nmin1(N_x-L1,N_y-L1) - E_y_nmin1_inc_back(N_y-L1) - E_y_nmin1(N_x-L+1,N_y-L1)) ;



  % ---------------------------- End H field update -----------------------------------------------------------------                  

  % ---------------------------- E field update -----------------------------------------------------------------                  
  % Update E fields: (note that latest H fields must be used!)
  E_x_n(1:N_x,2:N_y) = E_x_nmin1(1:N_x,2:N_y) ...
      + C_Ex(1:N_x,2:N_y).*(  H_z_n(1:N_x,2:N_y)  - H_z_n(1:N_x,1:N_y-1) ) ;
  E_y_n(2:N_x,1:N_y) = E_y_nmin1(2:N_x,1:N_y) ...
      - C_Ey(2:N_x,1:N_y).*(  H_z_n(2:N_x,1:N_y)  - H_z_n(1:N_x-1,1:N_y) ) ;
  
  % Special update on scat/tot field boundary (only needed for Ey) 
  % as in Fig. 3.1. The E_y field is the scattered field.
  H_z_n_inc = ones(1,N_y)*(Peak/eta_0)*gaussder_norm((m-1/2)*delta_t - (L-1/2)*delta_s/c,m_offset,sigma) ;
  E_y_n(L,L1:N_y-L1) = E_y_nmin1(L,L1:N_y-L1) ...
      - C_Ey(L,L1:N_y-L1).*(  H_z_n(L,L1:N_y-L1)  - H_z_n_inc(L1:N_y-L1) - H_z_n(L-1,L1:N_y-L1) ) ;
  
  % Special update on additional new scat/tot field boundary (only needed for Ey) 
  % on the right hand side of Fig. 3.1, at N_x - L1. Again, the E_y field is the
  % scattered field.
  H_z_n_inc = ones(1,N_y)*(Peak/eta_0)*gaussder_norm((m-1/2)*delta_t - (N_x-L1-1/2)*delta_s/c,m_offset,sigma) ;
  E_y_n(N_x-L1+1,L1:N_y-L1) = E_y_nmin1(N_x-L1+1,L1:N_y-L1) ...
      - C_Ey(N_x-L1+1,L1:N_y-L1).*(  H_z_n(N_x-L1+1,L1:N_y-L1)  + H_z_n_inc(L1:N_y-L1) - H_z_n(N_x-L1,L1:N_y-L1) ) ;
  % Special update on additional new scat/tot field boundary (only needed for Ex) 
  % on the upper side of Fig. 3.1, at N_y - L1. This plays the same role as
  % the ABC at N_x-L1. 
  for ii=1:N_x,
    % Re-writing gaussder_norm to handle a vector call would permit this to
    % be recoded in vector notation. Note that since N_x may not be equal to N_y, we defined a new variable  H_z_n_inc2
    H_z_n_inc2(ii,1) = (Peak/eta_0)*gaussder_norm((m-1/2)*delta_t - (ii-1/2)*delta_s/c,m_offset,sigma) ;
  end;
  % Upper interface
  % The H_x field above the interface is scattered, that below, total. The
  % E_x field is the scattered field.
  E_x_n(L:N_x-L1,N_y-L1+1) = E_x_nmin1(L:N_x-L1,N_y-L1+1) ... 
      + C_Ex(L:N_x-L1,N_y-L1+1).*(  H_z_n(L:N_x-L1,N_y-L1+1)  -  H_z_n(L:N_x-L1,N_y-L1) +  H_z_n_inc2(L:N_x-L1,1) ) ;
  % Lower interface
  % Now, the H_x field above the interface is total, that below, scattered.
  % Again, the E_x field is the scattered field.
  E_x_n(L:N_x-L1,L1)       = E_x_nmin1(L:N_x-L1,L1) ... 
      + C_Ex(L:N_x-L1,L1)      .*(  H_z_n(L:N_x-L1,L1)        -  H_z_n_inc2(L:N_x-L1,1) -  H_z_n(L:N_x-L1,L1-1)     ) ; 
  
  % ---------------------------- End E field update -----------------------------------------------------------------                  
  
  % ---------------------------- ABC exterior treatment -----------------------------------------------------------------                  
  % Impose ABC on sides - assumes free space cell on boundary.
  % Left/right boundaries:
  E_y_n(1,1:N_y+1)   = E_y_nmin1(1,1:N_y+1)*(1-c*delta_t/delta_s)   + c*delta_t/delta_s*E_y_nmin1(2,1:N_y+1); 
  E_y_n(N_x+1,1:N_y+1) = E_y_nmin1(N_x+1,1:N_y+1)*(1-c*delta_t/delta_s) + c*delta_t/delta_s*E_y_nmin1(N_x,1:N_y+1); 
  % Top/bottom boundaries:
  E_x_n(1:N_x+1,1)   = E_x_nmin1(1:N_x+1,1)*(1-c*delta_t/delta_s)   + c*delta_t/delta_s*E_x_nmin1(1:N_x+1,2); 
  E_x_n(1:N_x+1,N_y+1) = E_x_nmin1(1:N_x+1,N_y+1)*(1-c*delta_t/delta_s) + c*delta_t/delta_s*E_x_nmin1(1:N_x+1,N_y); 
  % ---------------------------- End ABC exterior treatment -----------------------------------------------------------------                  
  
  % Fix outer values of E_tangential as PEC:
  %E_y_n(1,:) = 0;
  %E_y_n(N_x,:) = 0;
  %E_x_n(:,1) = 0;
  %E_x_n(:,N_y) = 0;
  
  % Store data
  % Movie
  if mod(m,movie_interval) == 0 
    mesh(eta_0*H_z_n) % Normalize
    title(strcat('\eta_o H_z field at timestep ',num2str(m)))
    H_z_Movie(movie_count) = getframe;
    mesh(E_x_n) 
    title(strcat('E_x field at timestep ',num2str(m)))
    E_x_Movie(movie_count) = getframe;
    mesh(E_y_n) 
    title(strcat('E_y field at timestep ',num2str(m)))
    E_y_Movie(movie_count) = getframe;      
    movie_count = movie_count +1;
    % pause
  end
  % Time history
  H_z_point1(m) = H_z_n(point1_x,point1_y);
  H_z_point2(m) = H_z_n(point2_x,point2_y);
  E_y_point1(m) = E_y_n(point1_x,point1_y);
  E_y_point2(m) = E_y_n(point2_x,point2_y);

  % Update for next iteration
  H_z_nmin1 = H_z_n;
  E_x_nmin1 = E_x_n;
  E_y_nmin1 = E_y_n;
  
  % Output some indication of how far the code is:
  % disp('.')
  if(rem(m,report_time_interval)==0)
     m
  end 
  
  %  pause;
end
% End of time stepping
%------------------------------------------
disp('CPU time in seconds for time stepping')
run_time = cputime-start_time % seconds
floprate = 11*N_x*N_y*(M-1)/run_time

% movie(H_z_Movie,1,4);

time=[0:M-1]*delta_t;
plot(time/1e-9,E_y_point1)
xlabel('Time [ns]')
ylabel('E_y [V/m]')


print -deps RFGaussDerSimABC
disp('Press any key to continue')
pause;

% Compute and store time histories of source pulse

for mm = 1:M
  t(mm) = delta_t*(mm-1); 
  GaussDerSource(mm) = gaussder_norm(t(mm),m_offset,sigma);
end

delta_f = 1/((M-1)*delta_t);
freq=delta_f*(0:M-1);

clf;
run_movie = input('Run movie? (y/n)','s')    % Movie not always wanted, and difficult to stop. 
if run_movie =='y' | run_movie =='Y'
   movie(H_z_Movie,1,1)
end

% Remove large variables from workspace before saving
clear C_Ex C_Ey D_Hz 
clear E_x_n E_x_nmin1 E_y_n E_y_nmin1 H_z_n H_z_nmin1 mesh_pic

save cyl_fdtd_data 

