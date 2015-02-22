function [V_fem,z] = FEM_pp(V,N_elem,ell,h,N_int)
% FEM_pp performs some post-processing using the FEM solution.
% The FEM basis functions are used to interpolate the solution at N_int points.
% Note that because solution is continuous at nodes, values at nodes can be 
% computed using either element.
z = linspace(0,ell,N_int);
for ii = 1:N_int
    for jj=1:N_elem
        z_l = (jj-1) * h;
        z_r = jj *h;
        if z_l <= z(ii) && z(ii) <= z_r
            V_fem(ii) = V(jj) * (z_r-z(ii))/h + V(jj+1) * (z(ii)-z_l)/h;
        end
    end
end
