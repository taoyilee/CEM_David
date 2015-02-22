% A script to plot the convergence of the RCS of a PEC sphere,
% and illustrate the use of function "sphereRCS". The resulting plot
% is similar to Fig. 6.5
%
% See Chapter 6, Section 6.3, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP 2005.
%
% Written by DB Davidson, 10 Nov 2003. Revised 2 Feb 2005.

clear;
clf;
a = 1;
k = [0.1:0.05:10];
ka=k*a;
RCS_1 = sphereRCS(a,k,1)/(pi*a^2);
RCS_5 = sphereRCS(a,k,5)/(pi*a^2);
RCS_10 = sphereRCS(a,k,10)/(pi*a^2);
RCS_50 = sphereRCS(a,k,50)/(pi*a^2);
lambda=2*pi./k;
semilogy(ka,RCS_1,'k:',ka,RCS_5,'k--',...
         ka,RCS_10,'k-.',ka,RCS_50,'k-');
%semilogy(a./lambda,RCS_1,'k:',a./lambda,RCS_5,'k--',...
%         a./lambda,RCS_10,'k-.',a./lambda,RCS_50,'k-');
%xlabel('a/\lambda')
xlabel('ka')
ylabel('\sigma/\pi a^2')
legend('N=1','N=5','N=10','N=50',4)
print -deps cnvg_sph
