% Two-dimensional TEM mode computation using a triangular mesh
% Usage: Static2D
% See % D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", Chapter 10 of 2nd edn, presently in
% preparation.
%
% Mesh density is controlled by the 3rd and 4th arguments passed to
% function trimesh. The mesher is very specialized, and has only
% been tested with even-numbered values (4,2); (8,4) etc.
% Execution time and memory requirements grow very rapidly with mesh size, since the code
% has not been optimised at all w.r.t. speed or memory usage.
%
% Shares some functions with Eigen2D, in particular the mesher. 
% Written by D B Davidson, 14 August 2009.

clear all;
close all;

global ELEMENTS NODE_COORD NUM_NODES NUM_ELEMENTS NUM_DOFS
eps_0 = 8.854e-12;
mu_0=4*pi*1e-7;
mu_r = 1.0;

%a=5; % box width in units (x direction)
a=2.5;
b=a/2; % box height in units (y direction)
h=0.5; %substrate height in units
w=h;    % centre conductor thickness in units
%eps_r_sub = 1 % Air dielectric
eps_r_sub = 2.55 % Substrate dielectric constant as in static eg
%eps_r_sub =  9; % Substrate dielectric constant as in dynamic eg
V = 1; % voltage on centre conductor

x_mesh = 20; % Mesh in x
y_mesh = 20; % Ditto y
[x_nodes,y_nodes] = trimesh(a/2,b,x_mesh,y_mesh); % generate triangular mesh

% Plot mesh, with node and element numbers.
triplot(ELEMENTS,NODE_COORD(:,1),NODE_COORD(:,2));
axis([0 a/2 0 b]);
for inode = 1:NUM_NODES
    text(NODE_COORD(inode,1),NODE_COORD(inode,2),num2str(inode))
end
for ielem = 1:NUM_ELEMENTS
    x_c(ielem) = mean(NODE_COORD(ELEMENTS(ielem,1:3),1));
    y_c(ielem) = mean(NODE_COORD(ELEMENTS(ielem,1:3),2));
    text(x_c(ielem),y_c(ielem),num2str(ielem));
end

% Set up element material properties.
eps_r = ones(NUM_ELEMENTS,1); % Default material is free space
for ielem = 1:NUM_ELEMENTS
    if y_c(ielem) < h % height of substrate
        eps_r(ielem) = eps_r_sub;
    end
end
print -deps FEM2D_msh

% Plot mesh, with material properties
figure
triplot(ELEMENTS,NODE_COORD(:,1),NODE_COORD(:,2));
axis([0 a/2 0 b]);
for ielem = 1:NUM_ELEMENTS
    text(x_c(ielem),y_c(ielem),num2str(eps_r(ielem)));
end
%axis equal
%axis image
print -deps FEM2D_msh_ms


[node_free_flag,num_free_nodes] = free_nodes_mstrip(a,b,h,w);
node_pre_flag = not(node_free_flag); % Prescribed nodes
%dof_n1 = renumber_dof_nodes(node_free_flag);
[node_prenz_flag,num_prenz_nodes] = prescr_nodes_mstrip(a,b,h,w);
%pre_n1 = renumber_pre_nodes(node_pre_flag);
phi_pre = zeros(NUM_NODES,1); % Zero all (possibly) prescribed nodes
for inode =1:NUM_NODES
    if node_prenz_flag(inode)
        phi_pre(inode) = V; % Set non-zero prescribed nodes to their value
    end
end

S_mat = zeros(NUM_NODES,NUM_NODES);
b_vec = zeros(NUM_NODES,1);

for ielem=1:NUM_ELEMENTS % Assemble by elements, using alternate approach to handling prescribed boundaries.
    %for ielem=1:1
    trinodes = ELEMENTS(ielem,:);
    [S_elem] = s_nodal( NODE_COORD(trinodes(1),1),NODE_COORD(trinodes(1),2),...
        NODE_COORD(trinodes(2),1),NODE_COORD(trinodes(2),2),...
        NODE_COORD(trinodes(3),1),NODE_COORD(trinodes(3),2) );
    S_elem = S_elem*eps_r(ielem); % Scale by material properties
    % Assemble into global matrix.
    for jnode = 1:3
        jj = ELEMENTS(ielem,jnode);
        for knode = 1:3
            kk = ELEMENTS(ielem,knode);
            if node_free_flag(jj) && node_free_flag(kk)  % i.e. both free
                S_mat(jj,kk) = S_mat(jj,kk)+S_elem(jnode,knode);
            elseif node_free_flag(jj) && node_pre_flag(kk) % i.e. one free, one prescribed
                b_vec(jj) = b_vec(jj) - S_elem(jnode,knode)*phi_pre(kk);
                S_mat(kk,kk) =1;
                b_vec(kk) = phi_pre(kk);
                % Note - since the loop cycles through all values of jnode
                % and knode, the case where the jj-th node is prescribed
                % and the kk-th node is free will be taken care of 
                % subsequently, or has already been dealt with.
            end
        end
    end

end
phi_tot = sparse(S_mat)\b_vec;

% Create a 2D matrix of potential values
phi_mat = zeros(x_mesh,y_mesh);
for yy = 1:y_mesh+1
    for xx = 1:x_mesh+1
        phi_mat(xx,yy) = phi_tot(xx+(x_mesh+1)*(yy-1));
    end
end
XX = linspace(0,a/2,x_mesh+1);
YY = linspace(0,b,y_mesh+1);

figure
contour(XX,YY,phi_mat')
hold
[Ex,Ey]=gradient(phi_mat'); % Note that gradient operates across columns, then rows.
quiver(XX',YY',-Ex,-Ey)
xlabel('x (units)')
ylabel('y (units)')
print -deps FEM2D_ms_stat_1st

% Compute capacitance per unit length from E=1/2 C V^2 and save results for
% subsequent extrapolation.

E_tot = 0;
for ielem=1:NUM_ELEMENTS % Assemble by elements, using direct approach to handling prescribed boundaries.
    %for ielem=1:1
    trinodes = ELEMENTS(ielem,:);
    [S_elem] = s_nodal( NODE_COORD(trinodes(1),1),NODE_COORD(trinodes(1),2),...
        NODE_COORD(trinodes(2),1),NODE_COORD(trinodes(2),2),...
        NODE_COORD(trinodes(3),1),NODE_COORD(trinodes(3),2) );
    S_elem = S_elem*eps_r(ielem); % Scale by material properties
    phi_elem =[phi_tot(trinodes(1)) phi_tot(trinodes(2)) phi_tot(trinodes(3))]';
    E_elem = 0.5*phi_elem'*S_elem*phi_elem;
    E_tot = E_tot+E_elem;
end
E_tot = 2*eps_0*E_tot; % factor of 2 due to symmetry
if eps_r_sub == 1
    C0 = 2*E_tot/V^2;
    savefile=['Data_C0_',num2str(NUM_ELEMENTS),'_aeq',num2str(round(a)),'_beq',num2str(round(b))]
    save(savefile,'C0','eps_r','a','b','w','h','NUM_ELEMENTS','NUM_NODES','x_mesh','y_mesh')
else
    C = 2*E_tot/V^2;
    savefile=['Data_C_',num2str(NUM_ELEMENTS),'_aeq',num2str(round(a)),'_beq',num2str(round(b))]
    save(savefile,'C','eps_r','a','b','w','h','NUM_ELEMENTS','NUM_NODES','x_mesh','y_mesh')
end
