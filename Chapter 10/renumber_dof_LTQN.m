function [dof_e1,dof_e2,dof_f1,dof_f2] = renumber_dof_LTQN(dof_free_flag)
% RENUMBER_DOF renumbers degrees of freedom on the free edges 
% Author: D B Davidson, August 2009.
global NUM_EDGES NUM_DOFS NUM_ELEMENTS
dof_e1 = zeros(1,NUM_EDGES); % 0 indicates prescribed dof.
counter = 0; 
for i_edge = 1:NUM_EDGES
   if (dof_free_flag(i_edge))
     counter = counter + 1;
     dof_e1(i_edge) = counter;
   end 
end
dof_e2 = zeros(1,NUM_EDGES); % 0 indicates prescribed dof.
for i_edge = 1:NUM_EDGES
   if (dof_free_flag(i_edge))
     counter = counter + 1;
     dof_e2(i_edge) = counter;
   end 
end
for i_elem = 1:NUM_ELEMENTS
    counter = counter +1;
    dof_f1(i_elem)= counter;
end 
for i_elem = 1:NUM_ELEMENTS
    counter = counter +1;
    dof_f2(i_elem)= counter;
end 
NUM_DOFS= counter;