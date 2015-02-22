function rel_err = TEeig_err(a,b,TEeigvalues,N,num_zero_eig)
% TEEIG_ERR computes the error in the first N "non-zero" eigenvalues.
% Compute the first 100 TE eigenvalues (for safety) 
ii = 1;
for mm= 0:10, % 
  for nn = 0:9, % 
    if (mm == 0 & nn==0) 
        % do nothing
    else    
        k_c(ii) = sqrt((mm*pi/a)^2+(nn*pi/b)^2);
        ii = ii + 1;
    end
  end
end
k_c_sort = sort(k_c);
TEexact = k_c_sort(1:N);
rel_err = abs((TEexact'-TEeigvalues(num_zero_eig+1:num_zero_eig+N))./TEexact');
% keyboard
