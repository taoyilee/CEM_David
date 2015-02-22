function dof_flag = free_dof(a,b)
% FREE_DOF flags the degrees of freedom as free (1) or prescribed (0). 
% Author: D B Davidson, October 2003.
global NUM_EDGES EDGES  NODE_COORD
dof_flag(1:NUM_EDGES) = 1; % Flag all as free to start.

for i_edge = 1:NUM_EDGES
   node1 = EDGES(i_edge,1);
   node2 = EDGES(i_edge,2);
   if abs(NODE_COORD(node1,2))<eps  & abs(NODE_COORD(node2,2)) <eps % ie y=0
       dof_flag(i_edge) = 0;
   end
   if abs(NODE_COORD(node1,2)-b)<eps  & abs(NODE_COORD(node2,2)-b) <eps % ie y=b
       dof_flag(i_edge) = 0;
   end
   if abs(NODE_COORD(node1,1))<eps  & abs(NODE_COORD(node2,1)) <eps % ie x=0
       dof_flag(i_edge) = 0;
   end
   if abs(NODE_COORD(node1,1)-a)<eps  & abs(NODE_COORD(node2,1)-a) <eps % ie x=a
       dof_flag(i_edge) = 0;
   end
end

       
   
    