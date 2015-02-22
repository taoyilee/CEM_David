function edge_free_flag = free_dof3D(a,b,c,NODE_COORD,symmetry_flag,toler)
% FREE_DOF3D flags the edges  as free (1) or prescribed (0).
% PEC walls are assumed at x=0,y=0,z=0 and z=c.
% Magnetic symmetry is assumed at x=a, y=b.
% Author: D B Davidson, Sept 2009.
global NUM_EDGES EDGES
edge_free_flag(1:NUM_EDGES) = 1; % Flag all as free to start.

for i_edge = 1:NUM_EDGES
    node1 = EDGES(i_edge,1);
    node2 = EDGES(i_edge,2);
    if abs(NODE_COORD(node1,1))<toler  & abs(NODE_COORD(node2,1)) <toler % ie x=0
        edge_free_flag(i_edge) = 0;
    end
    if ~symmetry_flag
        if abs(NODE_COORD(node1,1)-a)<toler  & abs(NODE_COORD(node2,1)-a) <toler % ie x=a
            edge_free_flag(i_edge) = 0;
        end
    end
    if abs(NODE_COORD(node1,2))<toler  & abs(NODE_COORD(node2,2)) <toler % ie y=0
        edge_free_flag(i_edge) = 0;
    end
    if ~symmetry_flag
        if abs(NODE_COORD(node1,2)-b)<toler  & abs(NODE_COORD(node2,2)-b) <toler % ie y=b
            edge_free_flag(i_edge) = 0;
        end
    end
    if abs(NODE_COORD(node1,3))<toler  & abs(NODE_COORD(node2,3)) <toler % ie z=0
        edge_free_flag(i_edge) = 0;
    end
    if abs(NODE_COORD(node1,3)-c)<toler  & abs(NODE_COORD(node2,3)-c) <toler % ie z=c
        edge_free_flag(i_edge) = 0;
    end
end
