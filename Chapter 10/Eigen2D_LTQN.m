% Two-dimensional eigenmode computation using
% LTQN elements and a triangular mesh
% for hollow waveguide with PEC walls.
% Generates plots similar to Figs. 9.9 and 9.10 of 1st edn.
% Usage: Eigen2D
% See Chapter 10, D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", CUP, 2nd edn, presently in
% preparation.
%
% Mesh density is controlled by the 3rd and 4th arguments passed to
% function trimesh. The mesher is very specialized, and has only
% been tested with even-numbered values (4,2); (8,4) etc.
% Execution time and memory requirements grow very rapidly with mesh size, since the code
% has not been optimised at all w.r.t. speed or memory usage.
%
% Written by D B Davidson, 18 August 2009

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

[x_nodes,y_nodes] = trimesh(a,b,4,2); % generate triangular mesh

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
edge_free_flag=free_dof(a,b);
[dof_e1,dof_e2,dof_f1,dof_f2] = renumber_dof_LTQN(edge_free_flag);
%keyboard;
S = zeros(NUM_DOFS,NUM_DOFS);
T = zeros(NUM_DOFS,NUM_DOFS);
for ielem=1:NUM_ELEMENTS % Assemble by elements
    %for ielem=1:1
    trinodes = ELEMENTS(ielem,:);
    [S_elem,T_elem] = sandt_LTQN( NODE_COORD(trinodes(1),1),NODE_COORD(trinodes(1),2),...
        NODE_COORD(trinodes(2),1),NODE_COORD(trinodes(2),2),...
        NODE_COORD(trinodes(3),1),NODE_COORD(trinodes(3),2) ) ;
    %    keyboard
    % Assemble into global matrix.
    ll_f1=dof_f1(ielem);
    ll_f2=dof_f2(ielem);
    for jedge = 1:3
        jj_e1 = dof_e1(ELEMENT_EDGES(ielem,jedge));
        jj_e2 = dof_e2(ELEMENT_EDGES(ielem,jedge));
        % Edge-edge contributions
        for kedge = 1:3
            kk_e1 = dof_e1(ELEMENT_EDGES(ielem,kedge));
            kk_e2 = dof_e2(ELEMENT_EDGES(ielem,kedge));
            if jj_e1 && kk_e1  % i.e. both free  (same holds for e2)
                % Edge-edge and edge-face contributions (upper triangle)
                S(jj_e1,kk_e1) = S(jj_e1,kk_e1)+S_elem(jedge,kedge);
                S(jj_e1,kk_e2) = S(jj_e1,kk_e2)+S_elem(jedge,kedge+3);
                S(jj_e2,kk_e2) = S(jj_e2,kk_e2)+S_elem(jedge+3,kedge+3);

                T(jj_e1,kk_e1) = T(jj_e1,kk_e1)+T_elem(jedge,kedge);
                T(jj_e1,kk_e2) = T(jj_e1,kk_e2)+T_elem(jedge,kedge+3);
                T(jj_e2,kk_e2) = T(jj_e2,kk_e2)+T_elem(jedge+3,kedge+3);
            end
        end
        % Edge-face contributions
        if jj_e1 % i.e. edge free (same holds for e2, and ll faces are free
            S(jj_e1,ll_f1) = S(jj_e1,ll_f1)+S_elem(jedge,7);
            S(jj_e1,ll_f2) = S(jj_e1,ll_f2)+S_elem(jedge,8);
            S(jj_e2,ll_f1) = S(jj_e2,ll_f1)+S_elem(jedge+3,7);
            S(jj_e2,ll_f2) = S(jj_e2,ll_f2)+S_elem(jedge+3,8);

            T(jj_e1,ll_f1) = T(jj_e1,ll_f1)+T_elem(jedge,7);
            T(jj_e1,ll_f2) = T(jj_e1,ll_f2)+T_elem(jedge,8);
            T(jj_e2,ll_f1) = T(jj_e2,ll_f1)+T_elem(jedge+3,7);
            T(jj_e2,ll_f2) = T(jj_e2,ll_f2)+T_elem(jedge+3,8);

        end
    end
    % Face-face contributions (upper triangle)
    S(ll_f1,ll_f1) = S(ll_f1,ll_f1) + S_elem(7,7);
    S(ll_f1,ll_f2) = S(ll_f1,ll_f2) + S_elem(7,8);
    S(ll_f2,ll_f2) = S(ll_f2,ll_f2) + S_elem(8,8);
    T(ll_f1,ll_f1) = T(ll_f1,ll_f1) + T_elem(7,7);
    T(ll_f1,ll_f2) = T(ll_f1,ll_f2) + T_elem(7,8);
    T(ll_f2,ll_f2) = T(ll_f2,ll_f2) + T_elem(8,8);
%    keyboard
end
% Fill in other entries by symmetry
for jj=2:NUM_DOFS
    for kk=1:jj-1
        S(jj,kk) = S(kk,jj);
        T(jj,kk) = T(kk,jj);
    end
end

% Find eigensolution.
[eigvectors,eigvaluessq] = (eig(S/mu_r,T*eps_r));
[TEeigvalues,sort_vec] = sort(sqrt(diag((eigvaluessq))));

[node_flag,num_free_nodes] = free_nodes(a,b);
num_free_edges = nnz(edge_free_flag);
num_zero_eigvals = num_free_nodes+num_free_edges;
rel_err = TEeig_err(a,b,TEeigvalues,min(NUM_DOFS-num_zero_eigvals,8),num_zero_eigvals);

filenm=['eigdata_LTQN_',num2str(NUM_ELEMENTS)]
save  (filenm)
% Note - visualization of LTQN functions not implemented at present.
