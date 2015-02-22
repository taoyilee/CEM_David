function[vertices] = brick_mesh(x_len,y_len,z_len,N_x,N_y,N_z)
% This function returns a mesh of grid points, forming a brick (regular
% hexahedral) mesh. It is filled by column major order.
Delta_x = x_len/(N_x);
Delta_y = y_len/(N_y);
Delta_z = z_len/(N_z);
for kk = 1:N_z+1
    for jj = 1:N_y+1
        for ii = 1:N_x+1
            node_counter = ii+(jj-1)*(N_x+1)+(kk-1)*(N_y+1)*(N_x+1);
            vertices(node_counter,1) = (ii-1)*Delta_x;           
            vertices(node_counter,2) = (jj-1)*Delta_y;
            vertices(node_counter,3) = (kk-1)*Delta_z;
        end
    end
end
