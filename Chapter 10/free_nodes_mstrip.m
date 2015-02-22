function [node_flag,num_free_nodes] = free_nodes_mstrip(a,b,h,w)
% FREE_NODES_MSTRIP flags nodes as free (1) or prescribed (0)
% for boxed microstrip, assuming symmetry.
% Author: D B Davidson, August 2009.
% Note that it is implicitly assumed that a first order element, with nodes
% at the vertices only, is used.
global NUM_NODES NODE_COORD
node_flag(1:NUM_NODES) = 1; % Flag all as free to start.
for inode = 1:NUM_NODES,
    % PEC: y=0
    if abs(NODE_COORD(inode,2)) < eps % y=0
        node_flag(inode) = 0;
    end
    % PEC: x=a
    if abs(NODE_COORD(inode,1)-a/2) < eps % x=a/2
        node_flag(inode) = 0;
    end
    % PEC: y=b
    if abs(NODE_COORD(inode,2)-b) < eps % y=b
        node_flag(inode) = 0;
    end
    if abs(NODE_COORD(inode,2)-h) < eps && abs(NODE_COORD(inode,1)) <= w/2 % y=h and |x|<= w/2, ie on the centre conductor
        % (will later be flagged as inhomogenous prescribed)
        node_flag(inode) = 0;
    end
end
num_free_nodes=length(find(node_flag));
