% Three-dimensional eigenmode computation using
% Whitney elements and a tetrahedral mesh
% for a hollow cavity with PEC walls, with two planes of symmetry (magnetic wall) assumed
% at the central values of x and y.
% Usage: Eigen3D
% See D.B.Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", Chapter 11 of 2nd edn, presently in preparation.
% The mesh is can either be generated using an internal, or by reading in
% an externally-generated mesh .
% The internal mesh is generated using a two-step process. Firstly, a regular
% hexahedral (ie brick) mesh is generated. The vertices of this are then passed to a 3D Delaunay
% mesher, which returns a tetrahedral mesh (each hex element is split into
% six tets). Mesh density is controlled by N_x, N_y and N_z, which are the
% dimensions of the hex mesh.
% The externally-generated mesh is read in gmsh format. The x,y and z
% dimensions must be L_x, L_y and L_z respectively.

% Written by D B Davidson, 1 Sept 2009.
% Extensions DBD 9 Sept 2009, to read gmsh files.

clear all;
close all;

global ELEMENTS EDGES ELEMENT_EDGES NUM_NODES NUM_ELEMENTS NUM_EDGES LOCALEDGENODES  NUM_DOFS



eps_0 = 8.854e-12;
mu_0=4*pi*1e-7;
eps_r = 1.0;
mu_r = 1.0;
% Caution - do not change following SIX lines, hard-wired into other
% parts of code.
LOCALEDGENODES(1,:) = [1 2];
LOCALEDGENODES(2,:) = [1 3];
LOCALEDGENODES(3,:) = [1 4];
LOCALEDGENODES(4,:) = [2 3];
LOCALEDGENODES(5,:) = [2 4];
LOCALEDGENODES(6,:) = [3 4];

int_mesh_flag=input('Generate mesh file internally? 1 yes, 0 no (default yes)')
if isempty(int_mesh_flag)
    int_mesh_flag = 1
end
if ~int_mesh_flag
    %    mesh_filenm = 'cube.msh' % Presently hard-coded.
    mesh_filenm=input('Enter mesh file name: ','s')
    if isempty(mesh_filenm)
        mesh_filenm = 'box.msh' % box.msh is a 1x0.5x0.75 rectangular box.
        % cube.mesh cube is a 1x1x1m box.
    end
end

symmetry_flag = input('Use symmetry? 1 yes, 0 no (default no)');
if isempty(symmetry_flag)
    symmetry_flag = 0;
end

L_x=1; % x-dimensions of resonant cavity in m
L_y=0.5; % ditto y
L_z=0.75; % ditto z
% L_x = 1
% L_y = 1
% L_z = 1


if int_mesh_flag
    N_x=2; % Number of quadrilateral elements in x
    N_y=1; % ditto y
    N_z=2; % ditto z
    NUM_NODES=(N_x+1)*(N_y+1)*(N_z+1);
    % generate mesh points on a regular grid.
    if symmetry_flag
        vertices = brick_mesh(L_x/2,L_y/2,L_z,N_x,N_y,N_z);
    else
        vertices = brick_mesh(L_x,L_y,L_z,N_x,N_y,N_z);
    end

    if N_x==1 && N_y==1 && N_z==1  % Special treatment needed to avoid a problem in delaunay3
        TetMesh=delaunay3(vertices(:,1),vertices(:,2),vertices(:,3),{'QJ'});
    else
        TetMesh=delaunay3(vertices(:,1),vertices(:,2),vertices(:,3));
    end
else
    [NUM_NODES,x_nodes,y_nodes,z_nodes,line_counter,lines,tri_counter,triangles,tet_counter,TetMesh]=read_gmsh2(mesh_filenm);
    vertices(:,1) = x_nodes;
    vertices(:,2) = y_nodes;
    vertices(:,3) = z_nodes;
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
toler=min([L_x L_y L_z])*1e-5;
edge_free_flag=free_dof3D(L_x,L_y,L_z,vertices,symmetry_flag,toler);
[dof_e1,num_free_edges] = renumber_dof(edge_free_flag);
NUM_DOFS = num_free_edges;
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
[node_flag,num_free_nodes] = free_nodes3D(L_x,L_y,L_z,vertices,symmetry_flag,toler );

rel_err = eig_err3D(L_x,L_y,L_z,TEeigvalues,min(NUM_DOFS-num_free_nodes,8),num_free_nodes)
h=avg_mesh_length(vertices)
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
