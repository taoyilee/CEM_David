function rel_err = eig_err3D(a,b,d,eigvalues,N,num_zero_eig)
% EIG_ERR3D computes the error in the first N non-zero eigenvalues for the
% specific problem of a=1,b=0.5, d=0.75
% Note the modes may be either TEmnp or TMmnp, and some care is needed to
% establish the exact sequence of modes. See
% J.S.Savage, A.F.Peterson, "Higher-order vector finite elements for
% tetrahedral cells", IEEE MTT, June 1996.
eig_exact = [5.23599 7.02481 7.55145 7.55145 8.17887 8.17887 8.88577 8.94726];
rel_err = abs((eig_exact(1:N)'-eigvalues(num_zero_eig+1:num_zero_eig+N))./eig_exact(1:N)');
% keyboard
