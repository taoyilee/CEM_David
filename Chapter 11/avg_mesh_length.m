function [h] = avg_mesh_length(vertices)
% Function to find the average mesh length of a tetrahedral mesh, given the vertices.
global ELEMENTS LOCALEDGENODES NUM_ELEMENTS
for ielem = 1:NUM_ELEMENTS
    for jedge = 1:6
        vertex_1 = vertices(ELEMENTS(ielem,LOCALEDGENODES(jedge,1)),:);
        vertex_2 = vertices(ELEMENTS(ielem,LOCALEDGENODES(jedge,2)),:);        
        length(jedge) = norm(vertex_2-vertex_1);
    end
h_elem(ielem)=mean(length); % average mesh length per element
end
h=mean(h_elem); % average overall mesh length 
       