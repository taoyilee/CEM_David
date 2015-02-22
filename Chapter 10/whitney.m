function whitney_functions=whitney(elem_num,x,y)
% WHITNEY returns the Whitney basis functions w_ij
% in element i_elem at point x,y which is assumed to be within 
% element i_elem
% Note that the association
% edge 1: 12
% edge 2: 13
% edge 3: 23 
% is hard-wired into this function.
global ELEMENTS NODE_COORD
lambda  = simplex2D(elem_num,x,y);

trinodes = ELEMENTS(elem_num,:); 
x1 = NODE_COORD(trinodes(1),1);
y1 = NODE_COORD(trinodes(1),2);
x2 = NODE_COORD(trinodes(2),1);
y2 = NODE_COORD(trinodes(2),2);
x3 = NODE_COORD(trinodes(3),1);
y3 = NODE_COORD(trinodes(3),2);

%area = 0.5*abs(det([1 x1 y1; 1 x2 y2; 1 x3 y3]));
temp = inv([x1 x2 x3; y1 y2 y3; 1 1 1]);
b = temp(:,1);
c = temp(:,2);

nabla_lambda(1,:) = [b(1) c(1)];
nabla_lambda(2,:) = [b(2) c(2)];
nabla_lambda(3,:) = [b(3) c(3)];

whitney_functions(1,:) = lambda(1)*nabla_lambda(2,:) - lambda(2)*nabla_lambda(1,:);  
whitney_functions(2,:) = lambda(1)*nabla_lambda(3,:) - lambda(3)*nabla_lambda(1,:); 
whitney_functions(3,:) = lambda(2)*nabla_lambda(3,:) - lambda(3)*nabla_lambda(2,:);  
%keyboard