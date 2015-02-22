% Two-dimensional eigenmode computation using 
% Whitney elements and a triangular mesh
% for hollow waveguide with PEC walls. 
% Generates plots similar to Figs. 9.9 and 9.10 of 1st edn.
% Usage: Eigen2D
% See Chapter 9, Section 9.7, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP 2005, and Chapter 10 of 2nd edn, presently in
% preparation. 
%
% Mesh density is controlled by the 3rd and 4th arguments passed to
% function trimesh. The mesher is very specialized, and has only 
% been tested with even-numbered values (4,2); (8,4) etc. 
% Execution time and memory requirements grow very rapidly with mesh size, since the code
% has not been optimised at all w.r.t. speed or memory usage.
%
% Written by D B Davidson, October 2003. Revised 03 February 2005, and again for
% 2nd edn, 12 August 2009.

clear all;
close all;

global ELEMENTS NODE_COORD EDGES ELEMENT_EDGES NUM_NODES NUM_ELEMENTS NUM_EDGES LOCALEDGENODES EDGECONXELEM EDGECONXELEMEDGES ...
       NUM_DOFS
eps_0 = 8.854e-12;
mu_0=4*pi*1e-7;
eps_r = 1.0;
mu_r = 1.0;
% Caution - do not change following three lines, hard-wired into other 
% parts of code.
LOCALEDGENODES(1,:) = [1 2]; 
LOCALEDGENODES(2,:) = [1 3];
LOCALEDGENODES(3,:) = [2 3];
   
%a=20;
%b=10;

a=2.286e-2; % X band guide in cm
b=1.016e-2;
scaling=1; % to convert to SI units (unity in this case)

[x_nodes,y_nodes] = trimesh(a,b,8,4); % generate triangular mesh

triplot(ELEMENTS,NODE_COORD(:,1),NODE_COORD(:,2));
axis([0 a 0 b]);
for inode = 1:NUM_NODES
  text(NODE_COORD(inode,1),NODE_COORD(inode,2),num2str(inode))
end
for ielem = 1:NUM_ELEMENTS
  x_c = mean(NODE_COORD(ELEMENTS(ielem,1:3),1));
  y_c = mean(NODE_COORD(ELEMENTS(ielem,1:3),2));
  text(x_c,y_c,num2str(ielem));
end
print -deps tri_mesh
pause;
edgemake;
dof_e1_free_flag=free_dof(a,b);
dof_e1 = renumber_dof(dof_e1_free_flag);
S = zeros(NUM_DOFS,NUM_DOFS);
T = zeros(NUM_DOFS,NUM_DOFS);
for ielem=1:NUM_ELEMENTS % Assemble by elements
%for ielem=1:1
  trinodes = ELEMENTS(ielem,:); 
  [S_elem,T_elem] = sandt( NODE_COORD(trinodes(1),1),NODE_COORD(trinodes(1),2),...
                           NODE_COORD(trinodes(2),1),NODE_COORD(trinodes(2),2),...
                           NODE_COORD(trinodes(3),1),NODE_COORD(trinodes(3),2) ) ;
  % Assemble into global matrix.
  for jedge = 1:3
    jj = dof_e1(ELEMENT_EDGES(ielem,jedge));
    for kedge = 1:3      kk = dof_e1(ELEMENT_EDGES(ielem,kedge));
        if jj & kk  % i.e. both free 
          S(jj,kk) = S(jj,kk)+S_elem(jedge,kedge);
          T(jj,kk) = T(jj,kk)+T_elem(jedge,kedge);
        end
    end
  end
end 

% Find eigensolution.
[eigvectors,eigvaluessq] = (eig(S/mu_r,T*eps_r));
[TEeigvalues,sort_vec] = sort(sqrt(diag((eigvaluessq))));
[node_flag,num_free_nodes] = free_nodes(a,b);
rel_err = TEeig_err(a,b,TEeigvalues,min(NUM_DOFS-num_free_nodes,8),num_free_nodes); 

filenm=['eigdata_',num2str(NUM_ELEMENTS)]
save  (filenm)

XX = [0:0.05:1]*a;
YY = [0:0.1:1]*b;
for ii = 1:6 % Plot first six non-trivial eigenmodes
    eigmode = eigvectors(:,sort_vec(num_free_nodes+ii))
    dofs=eigmode;
    plot_field(dofs,dof_e1,XX,YY,3,2,ii,TEeigvalues(num_free_nodes+ii)/scaling);
    print ('-deps',filenm)
end
pause;
figure;
filenm2=['spur_eigdata_',num2str(NUM_ELEMENTS)]
for ii = 1:6 % Plot first six spurious eigenmodes
    eigmode = eigvectors(:,ii)
    dofs=eigmode;
    plot_field(dofs,dof_e1,XX,YY,3,2,ii,TEeigvalues(ii)/scaling);
    print ('-deps',filenm2)
end

