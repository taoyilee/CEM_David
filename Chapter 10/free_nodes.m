function [node_flag,num_free_nodes] = free_nodes(a,b)
% FREE_NODES flags nodes as free (1) or prescribed (0) for waveguide. 
% Author: D B Davidson, October 2003.
global NUM_NODES NODE_COORD
node_flag(1:NUM_NODES) = 1; % Flag all as free to start.
for inode = 1:NUM_NODES,
  % PEC: x=0
  if abs(NODE_COORD(inode,1)) < eps
    node_flag(inode) = 0;
  end  
  % PEC: y=0
  if abs(NODE_COORD(inode,2)) < eps
    node_flag(inode) = 0;
  end
  % PEC: x=a
  if abs(NODE_COORD(inode,1)-a) < eps
    node_flag(inode) = 0;
  end  
  % PEC: y=b
  if abs(NODE_COORD(inode,2)-b) < eps
    node_flag(inode) = 0;
  end
end
num_free_nodes=length(find(node_flag));
