function LTQN_functions=LTQN(lambda,nabla_lambda)
% LTQN returns the eight LTQN basis functions 
% given lambda and nabla_lamba.
% The first three entries are the Whitney functions, 
% the next three the LTLN functions, and the last two
% the "face" functions.
% Note that the association
% edge 1: 12
% edge 2: 13
% edge 3: 23 
% is hard-wired into this function.
LTQN_functions(1,:) = lambda(1)*nabla_lambda(2,:) - lambda(2)*nabla_lambda(1,:);  
LTQN_functions(2,:) = lambda(1)*nabla_lambda(3,:) - lambda(3)*nabla_lambda(1,:); 
LTQN_functions(3,:) = lambda(2)*nabla_lambda(3,:) - lambda(3)*nabla_lambda(2,:);  
LTQN_functions(4,:) = lambda(1)*nabla_lambda(2,:) + lambda(2)*nabla_lambda(1,:);  
LTQN_functions(5,:) = lambda(1)*nabla_lambda(3,:) + lambda(3)*nabla_lambda(1,:); 
LTQN_functions(6,:) = lambda(2)*nabla_lambda(3,:) + lambda(3)*nabla_lambda(2,:);  
LTQN_functions(7,:) = lambda(2)*lambda(3)*nabla_lambda(1,:)+ lambda(1)*lambda(3)*nabla_lambda(2,:)...
    -2*lambda(1)*lambda(2)*nabla_lambda(3,:);
LTQN_functions(8,:) = lambda(3)*lambda(1)*nabla_lambda(2,:)+ lambda(2)*lambda(1)*nabla_lambda(3,:)...
    -2*lambda(2)*lambda(3)*nabla_lambda(1,:);
