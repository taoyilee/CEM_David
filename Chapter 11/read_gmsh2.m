function[num_nodes,x_nodes,y_nodes,z_nodes,line_counter,lines,tri_counter,triangles,tet_counter,tets]=read_gmsh2(flname)
% Function to read mesh in from gmsh file, version 2 format.

% Original version DB Davidson, Univ Stellenbosch, March
% 2008. 1D, 2D and 3D simplicial elements read (ie lines, triangles, tets)
% and stored. (Note - tet code not tested).
% Corrections to tet code 08 Sept 2009, DBD.

%clear x_nodes y_nodes lines triangles tets

%flname='xband.msh'
fid=fopen(flname);
temp=textscan(fid,'%u',1,'Headerlines',4);
num_nodes=temp{1,1};
for inode=1:num_nodes
    node_data=textscan(fid,'%u %f %f %f',1);
    x_nodes(inode) = node_data{1,2};
    y_nodes(inode) = node_data{1,3};
    z_nodes(inode) = node_data{1,4};    
end;
temp=textscan(fid,'%u',1,'Headerlines',3);
num_elements=temp{1,1};
point_counter= 0;
line_counter = 0;
tri_counter =  0;
tet_counter =  0;
lines =     0; % Mesh files may contain only certain elements.
triangles = 0;
tets =      0;
for ii=1:num_elements
    elem_num=textscan(fid,'%u',1);
    elem_type=textscan(fid,'%u',1);
    if elem_type{1,1} == 15 % 1-node point
        elem_data=textscan(fid,'%u',4); % Read rest of line
        if elem_data{1,1} ~= 2;
            error('Only 2 tags supported in mesh file for 1-node points')
        end
        point_counter = point_counter +1;
        temp_point = elem_data{1,1};
        points(point_counter,1) = temp_point(end);
    elseif elem_type{1,1} == 1 % 1D line element
        elem_data=textscan(fid,'%u',6); % Read rest of line
        if elem_data{1,1} ~= 3
            error('Only default of 3 tags supported in mesh file')
        end
        line_counter = line_counter +1;
        temp_line = elem_data{1,1};
        lines(line_counter,1:2) = temp_line(end-1:end);
    elseif elem_type{1,1} == 2 % 2D triangular element
        elem_data=textscan(fid,'%u',7); % Read rest of line
        if elem_data{1,1} ~= 3;
            error('Only default of 3 tags supported in mesh file')
        end
        tri_counter = tri_counter +1;
        temp_tri= elem_data{1,1};
        triangles(tri_counter,1:3) = temp_tri(end-2:end);
    elseif elem_type{1,1} == 4 % 3D tetrahedral element  - note, code not tested.
        elem_data=textscan(fid,'%u',8); % Read rest of line
        if elem_data{1,1} ~= 3
            error('Only default of 3 tags supported in mesh file')
        end
        tet_counter = tet_counter +1;
        temp_tet = elem_data{1,1};
        tets(tet_counter,1:4) = temp_tet(end-3:end);
    end
end;

fclose(fid)
