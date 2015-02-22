function [dof_e1,dof_e2,dof_f1,dof_f2,num_free_edges,num_free_faces] = renumber_dof_LTQN(edge_free_flag,face_free_flag)
% RENUMBER_DOF renumbers the degrees of freedom for mixed 2nd order
% elements.
% Author: D B Davidson, October 2009.
global NUM_EDGES NUM_FACES
dof_e1 = zeros(1,NUM_EDGES); % 0 indicates prescribed dof in all cases.
dof_e2 = zeros(1,NUM_EDGES); 
dof_f1 = zeros(1,NUM_FACES); 
dof_f2 = zeros(1,NUM_FACES); 
counter = 0; 
num_free_edges = 0;
num_free_faces = 0;
for i_edge = 1:NUM_EDGES
   if (edge_free_flag(i_edge))
     num_free_edges = num_free_edges +1;
     counter = counter + 1;
     dof_e1(i_edge) = counter;
     counter = counter + 1;
     dof_e2(i_edge) = counter;
   end 
end
for i_face = 1:NUM_FACES
   if (face_free_flag(i_face))
     num_free_faces = num_free_faces +1;
     counter = counter + 1;
     dof_f1(i_face) = counter;
     counter = counter + 1;
     dof_f2(i_face) = counter;
   end 
end
