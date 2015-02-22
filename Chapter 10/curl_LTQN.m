function curl_LTQN_functions=curl_LTQN(lambda,nabla_lambda)
% LTQN returns the curl of the eight LTQN basis functions 
% given lambda and nabla_lamba.
% The first three entries are the Whitney functions, 
% the next three the LTLN functions, and the last two
% the "face" functions.
% Note that the association
% edge 1: 12
% edge 2: 13
% edge 3: 23 
% is hard-wired into this function.
curl_LTQN_functions(1,:) = 2*cross([nabla_lambda(1,:) 0],[nabla_lambda(2,:) 0]);  
curl_LTQN_functions(2,:) = 2*cross([nabla_lambda(1,:) 0],[nabla_lambda(3,:) 0]); 
curl_LTQN_functions(3,:) = 2*cross([nabla_lambda(2,:) 0],[nabla_lambda(3,:) 0]);  
curl_LTQN_functions(4,:) = 0;  
curl_LTQN_functions(5,:) = 0; 
curl_LTQN_functions(6,:) = 0;  
curl_LTQN_functions(7,:) = -3*lambda(2)*cross([nabla_lambda(1,:) 0],[nabla_lambda(3,:) 0])...
                           -3*lambda(1)*cross([nabla_lambda(2,:) 0],[nabla_lambda(3,:) 0]);
curl_LTQN_functions(8,:) = -3*lambda(3)*cross([nabla_lambda(2,:) 0],[nabla_lambda(1,:) 0])...
                           -3*lambda(2)*cross([nabla_lambda(3,:) 0],[nabla_lambda(1,:) 0]);
