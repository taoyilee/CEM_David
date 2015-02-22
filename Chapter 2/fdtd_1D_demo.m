% A MATLAB implementation of the 1D FDTD theory described in Chapter 2,
% Section 2.4, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP 2005, and 2edn, in preparation.
% Usage: fdtd_1D_demo
%
% The plot shows the time-domain evolution of the voltage at the load, 
% and then the Fourier transform results are plotted, to give results similar
% (but not identical) to Fig. 2.6. The code terminates on error criteria eps, or 
% when the number of time steps exceeds the value specified by the
% user via variable Nk.
%
% Note that the algorithm becomes unstable (due to the BC's) when 
% delta_t exceeds about 70% of the Courant limit. 
% Small values of Rl can also make the algorithm unstable. 
% Author: D.B.Davidson
% v1:   Originally written 23 Feb 2005.
% v2_2: Revised 8 March 2005, DBD and corrected 28 Feb 2006, DBD.
% v2_3: Corrected again 15 August 2007, DBD.

clear;
figure; % Permits comparative runs, by keeping previous results open in a figure window.

% h = 0.25 % length [m]
h = 0.25 % [m]
C = 1
L = 1
c = 1/sqrt(L*C)
Z_0 = sqrt(L/C)
Rs = 1 % Ohm
Rl = input('Load resistance? (Z_0 = 1 Ohm) Default: 2') 
if isempty(Rl)
    Rl = 2;
end
V0 = 1 % Amplitude of source voltage
%Rl = 1 % Ohm

freq = 4 % freq. of applied sinusoid [Hz]

Nz = 11;

delta_z = h/(Nz-1)  % Number of points. 

% Nk = 400 % Number of time steps for 0.25
T= 1/freq; % Period of signal
delta_f = 1/T;   % See p.58 for discussion of FFT parameters.
M = 64;    % Samples per period.
%delta_t = T/(M-1)
delta_t = T/M % Note slight change. The reason is that M in this case is the 
              % number of points PER PERIOD, and we need to ensure that the last time sample in one period
              % is not also the same time sample in the next period.
Nk=16*M;      % A maximum number of periods to run if the convergence criteria eps is not achieved.

growth = 1.5 % An abritary growth factor indicating instability
             % Don't make too close to 1, since the early time 
             % behaviour can indeed grow quite quickly.


beta1 = 2*delta_t/(Rs*C*delta_z)
beta2 = 2*delta_t/(Rl*C*delta_z)
r = (delta_t)^2/(L*C*delta_z^2)


% First time step - Initialize.
V_nmin1 = zeros(Nz,1); %
I_nmin1 = zeros(Nz,1); %

% Pre-allocation
V_n = zeros(Nz,1); %
I_n = zeros(Nz,1); %

V_period = zeros(M,Nz);
V_prev_period_freq = zeros(M,Nz);

movie_count = 1;
movie_interval = 10,


% Time loop
for nn = 2:Nk,
  Vo_nmin1 = V0*cos(2*pi*freq*(nn-1-1)*delta_t); % Source.
  V_n(1) = (1-beta1)*V_nmin1(1) - 2*I_nmin1(1) + (2/Rs)*Vo_nmin1;
% Loop code
%  for kk = 2:(Nz-1),
%    V_n(kk) = V_nmin1(kk) - (I_nmin1(kk) - I_nmin1(kk-1));
%  end
% Vector equivalent:
  V_n(2:(Nz-1)) = V_nmin1(2:(Nz-1)) - (I_nmin1(2:(Nz-1)) - I_nmin1(1:(Nz-2)));
  V_n(Nz) = (1-beta2)*V_nmin1(Nz) + 2*I_nmin1(Nz-1);

% Loop code
%  for kk = 1:Nz-1,
%    I_n(kk) = I_nmin1(kk) - r*(V_n(kk+1) - V_n(kk));
%  end
% Vector equivalent:
  I_n(1:Nz-1) = I_nmin1(1:Nz-1) - r*(V_n(2:Nz) - V_n(1:Nz-1));

  if mod(nn,movie_interval) == 0 
    plot(delta_t/(C*delta_z)*V_n) 
    axis([1 Nz -1 +1]);
    title(strcat('Voltage at timestep ',num2str(nn)))
    Voltage_Movie(movie_count) = getframe;
    movie_count = movie_count +1;
  end
  plot(I_n)
  % Check for possible instability (as a result of boundary conditions)
  if norm(V_n) > growth*norm(V_nmin1)
    if nn > 5 % Don't run test on first few passes.
      'Algorithm unstable. Time step',nn 
      break
    end 
  end
  % Update for next iteration
  V_nmin1 = V_n;
  I_nmin1 = I_n;
  index = mod(nn,M);
  if index == 0 
    index = M;
  end
  % Stores the data, overwriting with period M, for subsqeuent Fourier
  % transform.
  V_period(index,:) = V_n(:)'; % Transpose needed for this real-valued operation.
  if index == M % end of period 
    V_period_freq = fft(V_period); % Note: fft operates by column on multi-dimensional 
                                   % arrays and this must correspond to
                                   % time. NB! do NOT use fft2 here!
                                   % delta_t factor required to
                                   % get scaling correct
                                   % neglected by fft. See
                                   % discussion, p.58).                                          

    k=2; % See discussion toward end of file regarding the index k.
    eps = norm(V_period_freq(k,:) - V_prev_period_freq(k,:))/norm(V_period_freq(k,:))
    % Note that RMS norm includes inverse of root of length of vector, but it cancels above.  
    % The FFTs in the numerator and denominator of the above expression are
    % both unscaled - the scale factors cancel here. See later comments regarding correct scaling of the FFT.
    
    % Exit loop, or overwrite for next period:
    if eps < 0.002 
        break;
    end
    V_prev_period_freq = V_period_freq;
  end
end


% Get constants correct: 

V_n = delta_t/(C*delta_z)*V_n;
V_period_freq = delta_t/(C*delta_z)*V_period_freq; 

% Now compute exact reults (p.36) and compare to the simulated ones. 

% For the phasor results, what is needed is the first harmonic of the
% Fourier series expansion. The discussion in the textbook on p.58 refers to the general use of the
% FFT to approximate the Fourier transform. In this case, we are actually
% doing a Fourier series analysis. By sampling over exactly one period of the 
% (eventually) periodic signal, the DFT (and hence FFT) exactly represent
% the (sampled) Fourier series. See, for instance, the discussion in E.O.Brigham, "The
% Fast Fourier Transform and its applications", Prentice-Hall 1988,
% p.98-100 where this is discussed in detail.
% In this case, there is an additional factor of the period 1/T (present in the Fourier
% series coefficients, but absent in the Fourier integral) which must be
% taken in to account. 
% Furthermore, the factor of 2 comes from the negative and positive frequency
% components of the Fourier integral.

V_period_freq = (2*delta_t/T)*V_period_freq; % Scale as per above discussion

% Note also that in this case, the first (fundamental) harmonic is what is
% sought, at k=2. (k=1 is the DC term).

k = 2;
z= delta_z*(0:1:Nz-1);
lambda = c/freq;
beta = 2*pi/lambda;
Gamma = (Rl - Z_0)/ (Rl + Z_0); % Eq. (2.16)
V_plus = 1/2*V0; % for matched source, Eq. (2.15)
z_exact= (0:0.01:h);
V_exact = V_plus*( exp(-j*beta*(z_exact-h)) + Gamma * exp(j*beta*(z_exact-h))); % Eq. (2.14)

% Finally, plot results.
plot(z_exact,real(V_exact),'-', z_exact,imag(V_exact),'--',...
     z,real(V_period_freq(k,:)),'o',z,imag(V_period_freq(k,:)),'+');
legend('Real, exact','Imag, exact','Real, FDTD','Imag, FDTD',0);
xlabel('z (m)');
ylabel('Steady-state voltage (V)');
print -deps yeevex
