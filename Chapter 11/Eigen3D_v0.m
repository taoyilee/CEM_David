% Three-dimensional eigenmode computation using
% Whitney elements and a tetrahedral mesh
% for a hollow cavity with PEC walls, with two planes of symmetry (magnetic wall) assumed
% at the central values of x and y.
% Usage: Eigen3D
% See D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", Chapter 11 of 2nd edn, presently in preparation.
% The mesh is generated in a two-step process. Firstly, a regular
% hexahedralmesh is generated. The vertices of this are then passed to a 3D Delaunay
% mesher, which returns a tetrahedral mesh (each hex element is split into
% six tets). Mesh density is controlled by N_x, N_y and N_z, which are the
% mesh

% Written by D B Davidson, 1 Sept 2009.

clear all;
close all;



global ELEMENTS EDGES ELEMENT_EDGES NUM_NODES NUM_ELEMENTS NUM_EDGES LOCALEDGENODES  NUM_DOFS



eps_0 = 8.854e-12;
mu_0=4*pi*1e-7;
eps_r = 1.0;
mu_r = 1.0;
% Caution - do not change following three lines, hard-wired into other
% parts of code.
LOCALEDGENODES(1,:) = [1 2];
LOCALEDGENODES(2,:) = [1 3];
LOCALEDGENODES(3,:) = [1 4];
LOCALEDGENODES(4,:) = [2 3];
LOCALEDGENODES(5,:) = [2 4];
LOCALEDGENODES(6,:) = [3 4];

symmetry_flag = input('Use symmetry? 1 yes, 0 no (default no)');
if isempty(symmetry_flag)
    symmetry_flag = 0;
end

L_x=1; % x-dimensions of resonant cavity in m
L_y=0.5; % ditto y
L_z=0.75; % ditto z

N_x=8; % Number of quadrilateral elements in x
N_y=4; % ditto y
N_z=6; % ditto z

NUM_NODES=(N_x+1)*(N_y+1)*(N_z+1);

if symmetry_flag
    vertices = brick_mesh(L_x/2,L_y/2,L_z,N_x,N_y,N_z);
else
    vertices = brick_mesh(L_x,L_y,L_z,N_x,N_y,N_z);
end
% generate mesh points on a regular grid.

if N_x==1 && N_y==1 && N_z==1  % Special treatment needed to avoid a problem in delaunay3
    TetMesh=delaunay3(vertices(:,1),vertices(:,2),vertices(:,3),{'QJ'});
else
    TetMesh=delaunay3(vertices(:,1),vertices(:,2),vertices(:,3));
end
ELEMENTS=sort(TetMesh,2); % Sort nodes into ascending order for each element
NUM_ELEMENTS = length(ELEMENTS)
if symmetry_flag
    filenm=['eigdata3D_sym',num2str(NUM_ELEMENTS)]
else
    filenm=['eigdata3D_',num2str(NUM_ELEMENTS)]
end
tetramesh(TetMesh,vertices)
print('-depsc',filenm)
print('-dpng',filenm)

edgemake3D;
dof_e1_free_flag=free_dof3D(L_x,L_y,L_z,vertices,symmetry_flag);
dof_e1 = renumber_dof(dof_e1_free_flag);
NUM_DOFS
S = zeros(NUM_DOFS,NUM_DOFS);
T = zeros(NUM_DOFS,NUM_DOFS);
for ielem=1:NUM_ELEMENTS % Assemble by elements
    %for ielem=1:1
    tetnodes = ELEMENTS(ielem,:);
    [S_elem,T_elem] = sandt3D( vertices(tetnodes(1),:),vertices(tetnodes(2),:),...
        vertices(tetnodes(3),:),vertices(tetnodes(4),:) );
    % Assemble into global matrix.
    for jedge = 1:6
        jj = dof_e1(ELEMENT_EDGES(ielem,jedge));
        for kedge = 1:6
            kk = dof_e1(ELEMENT_EDGES(ielem,kedge));
            if jj & kk  % i.e. both free
                S(jj,kk) = S(jj,kk)+S_elem(jedge,kedge);
                T(jj,kk) = T(jj,kk)+T_elem(jedge,kedge);
            end
        end
    end
end
% Find eigensolution.
[eigvectors,eigvaluessq] = (eig(S/mu_r,T*eps_r));
% S=S/mu_r;
% T=T*eps_r;
% [eigvectors,eigvaluessq] = (eigs(sparse(S),sparse(T)));
%clear S T % clear matrics before saving system, or doing further work.
[TEeigvalues,sort_vec] = sort(sqrt(diag((eigvaluessq))));
[node_flag,num_free_nodes] = free_nodes3D(L_x,L_y,L_z,vertices,symmetry_flag);

rel_err = TEeig_err3D(L_x,L_y,L_z,TEeigvalues,min(NUM_DOFS-num_free_nodes,8),num_free_nodes)

save  (filenm)

% XX = [0:0.05:1]*a;
% YY = [0:0.1:1]*b;
% for ii = 1:6 % Plot first six non-trivial eigenmodes
%     eigmode = eigvectors(:,sort_vec(num_free_nodes+ii))
%     dofs=eigmode;
%     plot_field(dofs,dof_e1,XX,YY,3,2,ii,TEeigvalues(num_free_nodes+ii)/scaling);
%     print ('-deps',filenm)
% end
% pause;
% figure;
% filenm2=['spur_eigdata_',num2str(NUM_ELEMENTS)]
% for ii = 1:6 % Plot first six spurious eigenmodes
%     eigmode = eigvectors(:,ii)
%     dofs=eigmode;
%     plot_field(dofs,dof_e1,XX,YY,3,2,ii,TEeigvalues(ii)/scaling);
%     print ('-deps',filenm2)
% end
%
