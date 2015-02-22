% MATLAB implementation of FDTD 1D wideband application.
% DBD 13 May 2003
% Usage: fdtd_1D_WB_demo
% Note that the algorithm becomes unstable (due to the BC's) when
% delta_t exceeds about 70% of the Courant limit.
% Small values of Rl can also make the algorithm unstable.

% This code solves directly for the voltages and does not
% use the normalized voltage.

% _v1_0 Minor correction 15 Aug 2007, DBD.
% _v1_1 Minor improvements in user input, 26 Nov 2009, DBD.

clear;

C = 1;
L = 1;
Z_0 = sqrt(L/C);
c = 1/sqrt(L*C);
flag=input('Single pulse analysis ? (1 or true for TD, 0 or false or return for tx function)')
if isempty(flag)
    flat = 0;
end


%For time-domain demo, use Rs=1, Rl=3 for nice results. Start out with both
%=1
if flag
    %h = 0.25; % length [m]
    h=0.5
    Rs = 1; % Ohm
    Rl = input('Load resistance? Default 3 Ohm') % Ohm
    if isempty(Rl)
        Rl = 3;
    end
    Courant = 0.5 % 0.7;
    sigma = 1/(8*pi);
    Nk = 512 % Number of time steps. Power of two best (for FFT).
else % Transfer function. Use Rs=1/2, Rl=2, reduce Courant limit and run over wider band.
    h = 0.25; % length [m]
    Rs = 0.5; % Ohm
    Rl = 2; % Ohm
    Courant = 0.4
    sigma = 1/(16*pi); % more bandwidth
    Nk = 1024 % Number of time steps. Power of two best (for FFT).
end
offset = 4*sigma;

Nz = 11 % 41; % 11, 21 and 41 demonstrate effect of dispersion nicely.

delta_z = h/(Nz-1)  % Number of points.


delta_t = Courant*delta_z/c %

growth = 20 % An abritary growth factor indicating instability
% Don't make too close to 1, since the early time
% behaviour can indeed grow quite quickly.

beta1 = 2*delta_t/(Rs*C*delta_z);
beta2 = 2*delta_t/(Rl*C*delta_z);
const=delta_t/(C*delta_z);
% r = (delta_t)^2/(L*C*delta_z^2)


% First time step - Initialize.
V_nmin1 = zeros(Nz,1); %
I_nmin1 = zeros(Nz,1); %

% Pre-allocation
V_n = zeros(Nz,1); %
I_n = zeros(Nz,1); %
Vs = zeros(Nk,1); %
Vl = zeros(Nk,1); %
time = zeros(Nk,1); %


movie_count = 1;
movie_interval = 10,

% Time loop
for nn = 2:Nk,
    time(nn) = (nn-1)*delta_t;
    Vo_nmin1 = gaussder_norm((nn-1)*delta_t,offset,sigma);
    Vs(nn) = Vo_nmin1;
    V_n(1) = (1-beta1)*V_nmin1(1) - 2*const*I_nmin1(1) + (2*const/Rs)*Vo_nmin1;
    % Loop code
    %  for kk = 2:(Nz-1),
    %    V_n(kk) = V_nmin1(kk) - (I_nmin1(kk) - I_nmin1(kk-1));
    %  end
    % Vector equivalent:
    V_n(2:(Nz-1)) = V_nmin1(2:(Nz-1)) - const*(I_nmin1(2:(Nz-1)) - I_nmin1(1:(Nz-2)));
    V_n(Nz) = (1-beta2)*V_nmin1(Nz) + 2*const*I_nmin1(Nz-1);
    Vl(nn) = V_n(Nz);

    % Loop code
    %  for kk = 1:Nz-1,
    %    I_n(kk) = I_nmin1(kk) - r*(V_n(kk+1) - V_n(kk));
    %  end
    % Vector equivalent:
    I_n(1:Nz-1) = I_nmin1(1:Nz-1) - delta_t/(L*delta_z)*(V_n(2:Nz) - V_n(1:Nz-1));

    if mod(nn,movie_interval) == 0
        plot(delta_t/(C*delta_z)*V_n)
        axis([1 Nz -1 +1]);
        title(strcat('Voltage at timestep ',num2str(nn)))
        Voltage_Movie(movie_count) = getframe;
        movie_count = movie_count +1;
    end
    % plot(I_n)
    % Check for possible instability (as a result of boundary conditions)
    norm(V_n);
    norm(V_nmin1);
    if norm(V_n) > growth*norm(V_nmin1)
        if nn > 10 % Don't run test on first few passes.
            'Algorithm unstable. Time step',nn
            break
        end
    end
    % Update for next iteration
    V_nmin1 = V_n;
    I_nmin1 = I_n;
end

%axis([1 Nz -1 +1]);
%plot(V_n);
%pause;

save fdtd_1D_WB_demo;

if flag
    tau_s = Z_0/(Rs+Z_0);
    figure
    plot(time,Vs*tau_s,time,Vl,'--');
    legend('V_+','V_L');
    title('Positive-going voltage and load voltage.')
    xlabel('t [s]')
    ylabel('Voltage')
    pause;
else
    plot(time,Vl);
    title('Load voltage .')
    xlabel('t [s]')
    ylabel('Voltage')
    pause;
    TxFunc = (fft(Vl)./fft(Vs));
    delta_f = 1/(Nk*delta_t);
    freq = zeros(Nk/2,1);
    for kk=1:Nk/2
        freq(kk) = (kk-1)*delta_f;
    end
    plot(freq(1:Nk/32),abs(TxFunc(1:Nk/32)));
    axis([0 6 0 2])
    title('Transfer function V_l/V_o.')
    xlabel('Frequency [Hz]')
    ylabel('|V_l/V_o|')
end

% movie(Voltage_Movie,1,4); % play 1 times at 4 frames per second.


