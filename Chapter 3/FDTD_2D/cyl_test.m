
clear;
load cyl_fdtd_data_coarse_old;
E_y_point1_old = E_y_point1;
load cyl_fdtd_data_coarse_new;
E_y_point1_new = E_y_point1;
load cyl_fdtd_data_coarse_new_longer;
E_y_point1_new_longer = E_y_point1;

plot(time/1e-9,E_y_point1_old,time/1e-9,E_y_point1_new,':',time/1e-9,E_y_point1_new_longer,'--')
ylabel('E_y [V/m]')
xlabel('t [ns]')
legend('Old formulation N_x = 400','New formulation N_x = 400','New formulation N_x = 500')
disp('Press any key to continue')
pause;

