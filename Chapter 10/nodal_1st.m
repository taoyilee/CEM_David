function alpha=nodal_1st(elem_num,x,y)
% NODAL_1st returns the nodal basis functions alpha_i
% in element i_elem at point x,y which is assumed to be within 
% element i_elem

global ELEMENTS NODE_COORD
lambda  = simplex2D(elem_num,x,y);
alpha =lambda;
