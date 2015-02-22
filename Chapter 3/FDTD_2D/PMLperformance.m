% m-file to compute ABC performance from 200x200 and 400x400 grid.

clear;

load cyl_fdtd_data_pml200thin
E_200_thin = E_y_point1;
time_200 =time;
d_thin = d_cell;
poly_thin = poly_m;

load cyl_fdtd_data_pml200thick
E_200_thick = E_y_point1;
d_thick = d_cell;
poly_thick = poly_m;

load cyl_fdtd_data_pml400
E_400 = E_y_point1;
time_400 =time;
plot(1e9*time_200,E_200_thin,time_400,E_400,':')
T = length(time_200);
refl_thin = E_200_thin-E_400(1:T);
refl_thick = E_200_thick-E_400(1:T);
plot(time_200*1e9,refl_thin,time_200*1e9,refl_thick);
pause;
gate=round(0.75*T);
refl_thin_gated(1:gate) = refl_thin(1:gate);
refl_thin_gated_norm = refl_thin_gated./max(abs(E_400(1:T)));
refl_thick_gated(1:gate) = refl_thick(1:gate);
refl_thick_gated_norm = refl_thick_gated./max(abs(E_400(1:T)));
plot(1e9*time_200(1:gate),refl_thin_gated_norm,'k-',...
     1e9*time_200(1:gate),refl_thick_gated_norm,'k-.');
xlabel('t [ns]')
ylabel('\Gamma_{norm}')
legend(['d=',num2str(d_thin),'; m=',num2str(poly_thin)],...
       ['d=',num2str(d_thick),'; m=',num2str(poly_thick)])
axis([0 5 -10e-4  10e-4])
pause;
print -deps PMLgamma

plot(1e9*time_200(1:gate),20*log10(refl_thin_gated_norm),'k-',...
     1e9*time_200(1:gate),20*log10(refl_thick_gated_norm),'k-.');
axis([0 5 -100 -50])
xlabel('t [ns]')
ylabel('|\Gamma_{norm}| [dB]')
legend(['d=',num2str(d_thin),'; m=',num2str(poly_thin)],...
       ['d=',num2str(d_thick),'; m=',num2str(poly_thick)])
pause;
print -deps PMLgamma_dB


