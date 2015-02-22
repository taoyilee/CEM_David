% M-file to compute TM scattering from PEC cylinder using MoM MFIE.

% Author D B Davidson
% 3 March 2008.

% The code uses pulse basis functions, with either a simple mid-point
% evaluation for off-diagonal elements, or a simple quadrature scheme.
% Details of the MoM solution are in the function MoM_TM_solver.
% Reference: AF Peterson, S.L.Ray, R.Mittra, "Computational methods for
% electromagnetics", IEEE Press, 1998, Section 2.1.

% The code can either compute and display the surface current at one
% frequency, or can compute radar width as a function of frequency.

% Note that very small numbers of segments render the approximation of the
% geometry highly inaccurate. 

clear;
clf;

% Define geometry:

r2d    = 180/pi;
eps_0 = 8.854e-12;  % [F/m]
mu_0  = 4*pi*1e-7;  % [H/m]
eta_0 = sqrt(mu_0/eps_0); % Wave impedance of free space
E_0   = 1;
H_0   = E_0/eta_0;          % H field corresponding to incident E field.

analysis_type =input('Enter analysis type, 1 for current, 2 for RCS over frequency range')
switch (analysis_type)
    case {1}
        N = input('Enter number of segments, or return for default (10)') % Number of segments
        if (isempty(N))
            N=10;
        end
        lambda = 1          % free-space wavelength [m]
        k=2*pi/lambda       % free-space wavenumber [rad/m]
        a = 1*lambda/(2*pi) % cylinder radius in [m]
        phi_inc = pi;       % Angle incident field incident from.

        [J_s,phi_c]      = MoM_TM_solver(k,N,a,E_0,phi_inc,0,eta_0,1); % call MoM solver.
        [J_s_quad,phi_c] = MoM_TM_solver(k,N,a,E_0,phi_inc,1,eta_0,1); % call MoM solver.

%        ref_J_exact =       [0.8882 0.7050 0.8914 1.410 1.656 1.707]; % Data from Peterson et al, Table 2.5, p.49.
%        ref_J_MFIE_eq2_26 = [0.9742 0.7405 0.9363 1.517 1.778 1.829]; % Ditto

%        ref_J_ang = [0 36 72 108 144 180];

        plot(phi_c*r2d,abs(J_s)/H_0,'+',...
            phi_c*r2d,abs(J_s_quad)/H_0,'^')
        axis([0 180 0 3])
        xlabel('\Phi')
        %ylabel('Surface current density normalized by |H_{inc}|')
        ylabel('|J_s/H^{inc}|')
        legend('MoM EFIE (this code, single point)','MoM MFIE (this code, quadrature)','Location','Best')
        print -deps TM_Js_1
    case {2}
        N_per_lambda = input('Enter discretization per wavelength, or return for default (20)') % Number of segments
        if (isempty(N_per_lambda))
            N_per_lambda=10;
        end
        %lambda = [0.03:0.01:0.1]
        %k=2*pi./lambda;       % free-space wavenumber [rad/m]
        a = 0.03;            % cylinder radius in [m]
        k= [0.1:0.05:10]/a;      % free-space wavenumber [rad/m]
%        k= [0.3:1:10]/a;      % free-space wavenumber [rad/m]
        lambda = 2*pi./k;    % free-space wavelength [m]
        N_freq = length(k);
        phi_inc = pi;       % Angle of incident field incident from
        for ff = 1:N_freq
            disp(['Frequency ',num2str(ff),' of ',num2str(N_freq)]);
            N(ff)=max(ceil(2*pi*a/lambda(ff)*N_per_lambda),20);
            % Ensure that there are at least 20 segments even at low
            % frequencies to adequately resolve geometry.
            [J_s,     phi_c,w,x_c,y_c,cond_num] =MoM_TM_solver(k(ff),N(ff),a,E_0,phi_inc,0,eta_0,1);  % call MoM solver.
            cond_num_TM_sp(ff) = cond_num;
            [J_s_quad,phi_c,w,x_c,y_c,cond_num] =MoM_TM_solver(k(ff),N(ff),a,E_0,phi_inc,1,eta_0,1);  % call MoM solver with quadrature
            cond_num_TM_quad(ff) = cond_num;
            % RCS computation
            rcs_TM(ff) = 0;
            rcs_TM_quad(ff) = 0;
            phi_bistat = phi_inc; % For monostatic RCS
            for nn= 1:N(ff)
               rcs_TM(ff)      = rcs_TM(ff)      +      J_s(nn)*w*exp(j*k(ff)*(x_c(nn)*cos(phi_bistat)+y_c(nn)*sin(phi_bistat))) ;
               rcs_TM_quad(ff) = rcs_TM_quad(ff) + J_s_quad(nn)*w*exp(j*k(ff)*(x_c(nn)*cos(phi_bistat)+y_c(nn)*sin(phi_bistat))) ;
            end 
            rcs_TM(ff)      = ((k(ff)*eta_0^2)/4)*abs(rcs_TM(ff))^2     /abs(E_0)^2 ;
            rcs_TM_quad(ff) = ((k(ff)*eta_0^2)/4)*abs(rcs_TM_quad(ff))^2/abs(E_0)^2 ;
            clear J_s J_s_quad phi_c w x_c y_c Omega; % clear before next call.
        end
        loglog(k*a,rcs_TM_quad/(pi*a),'-',k*a,rcs_TM/(pi*a),'-.');
        xlabel('ka')
        ylabel('\sigma_{TM}/\pi a')
        legend('This code, eq.2.10','This code, eq.2.8','Location','Best')
        print -deps TM_echo
        pause; 
        N_Bessel=25; % number of terms in Bessel functions
        echo_width_exact_TM =cyl_TM_echo_width(a,k,N_Bessel,pi);
        plot(a./lambda,echo_width_exact_TM./lambda,'k-',a./lambda,rcs_TM_quad./lambda,'k:',a./lambda,rcs_TM./lambda,'k-.');
        xlabel('a/\lambda')
        ylabel('\sigma_{TM}/\lambda')
        legend('Eigenfunction','This code, quadrature','This code, single point','Location','Best')
        print -deps TM_echo2
        save MoM_TM_data % save for future use
     
    otherwise
        disp('Unknown method.')
end