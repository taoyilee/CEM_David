function [node_flag,num_free_nodes] = free_nodes3D(a,b,c,NODE_COORD,symmetry_flag,toler)
% FREE_NODES flags nodes as free (1) or prescribed (0) for a PEC cavity.
% PEC walls are assumed at x=0,y=0,z=0 and z=c.
% Magnetic symmetry is assumed at x=a, y=b if the symmetry flag is set,
% otherwise, PEC walls are assumed there too.
% Author: D B Davidson, Sept 2009.
global NUM_NODES
node_flag(1:NUM_NODES) = 1; % Flag all as free to start.
for inode = 1:NUM_NODES,
    % PEC: x=0
    if abs(NODE_COORD(inode,1)) < toler
        node_flag(inode) = 0;
    end
    if ~symmetry_flag
        if abs(NODE_COORD(inode,1)-a) < toler % x=a
            node_flag(inode) = 0;
        end
    end
    % PEC: y=0
    if abs(NODE_COORD(inode,2)) < toler
        node_flag(inode) = 0;
    end
    if ~symmetry_flag
        if abs(NODE_COORD(inode,2)-b) < toler % y=b
            node_flag(inode) = 0;
        end
    end
    % PEC: z=0
    if abs(NODE_COORD(inode,3)) < toler
        node_flag(inode) = 0;
    end
    % PEC: z=c
    if abs(NODE_COORD(inode,3)-c) < toler
        node_flag(inode) = 0;
    end

end
num_free_nodes=length(find(node_flag));
