function [dof_e1,num_free_edges] = renumber_dof(edge_free_flag)
% RENUMBER_DOF renumbers the free edges 
% (which are the degrees of freedom when using Whitney elements).  
% Author: D B Davidson, October 2003.
global NUM_EDGES 
dof_e1 = zeros(1,NUM_EDGES); % 0 indicates prescribed dof.
counter = 0; 
for i_edge = 1:NUM_EDGES
   if (edge_free_flag(i_edge))
     counter = counter + 1;
     dof_e1(i_edge) = counter;
   end 
end
num_free_edges = counter;
    