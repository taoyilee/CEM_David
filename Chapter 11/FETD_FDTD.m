% Analysis of CTLN functions for FETD-FDTD equivalence

y=[0:0.01:1];
we1 = 1-y;
we2=y;
we1_dot_we2 = we1.*we2;
plot(y,we1,'--',y,we2,'-.',y,we1_dot_we2,'-',quad_points,values,'ko');
xlabel('y')
ylabel('f(y)')
quad_points=[0,1];
values=[0,0]; % For plotting only
legend('w_{e1}','w_{e2}','w_{e1} \cdot w_{e2}')
pause;
print -dpng FETD_CTLN_funcs
