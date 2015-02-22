function [lambda,tri_area]=simplex2D(elem_num,xp,yp)
% SIMPLEX2D returns the simplex coordinates of the point
% with coordinates (xp,yp) within element elem_num
% as well as the element area

global ELEMENTS NODE_COORD 

trinodes = ELEMENTS(elem_num,:); 
x1 = NODE_COORD(trinodes(1),1);
y1 = NODE_COORD(trinodes(1),2);
x2 = NODE_COORD(trinodes(2),1);
y2 = NODE_COORD(trinodes(2),2);
x3 = NODE_COORD(trinodes(3),1);
y3 = NODE_COORD(trinodes(3),2);

sigma  = det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
sigma1 = det([1 xp yp; 1 x2 y2; 1 x3 y3]);
sigma2 = det([1 x1 y1; 1 xp yp; 1 x3 y3]);
sigma3 = det([1 x1 y1; 1 x2 y2; 1 xp yp]);
lambda(1) = sigma1/sigma;
lambda(2) = sigma2/sigma;
lambda(3) = sigma3/sigma;
tri_area = 0.5*abs(sigma);
