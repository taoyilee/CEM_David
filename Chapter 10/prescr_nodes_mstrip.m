function [node_flag,num_pre_nodes] = prescr_nodes_mstrip(a,b,h,w)
% PRESCRI_NODES_MSTRIP flags nodes as inhomogenous prescribed (1) for boxed microstrip,
% assuming symmetry.
% Author: D B Davidson, August 2009.
% Note that it is implicitly assumed that a first order element, with nodes
% at the vertices only, is used.
global NUM_NODES NODE_COORD
node_flag(1:NUM_NODES) = 0; % Flag all as not inhomogenous prescribed to start.

for inode = 1:NUM_NODES,
    if abs(NODE_COORD(inode,2)-h) < eps && abs(NODE_COORD(inode,1)) <= w/2 % y=h and |x|<= w/2, ie on the centre conductor
        node_flag(inode) = 1;
    end
end
num_pre_nodes=length(find(node_flag));
