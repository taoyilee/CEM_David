% Three-dimensional eigenmode computation using
% Whitney elements and a tetrahedral mesh
% for a hollow cavity with PEC walls, with two planes of symmetry (magnetic wall) assumed
% at the central values of x and y.
% Usage: Eigen3D_LTQN
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
% Extensions DBD 8 Oct 2009, adding option of using LT/QN elements.

clear all;
close all;

global ELEMENTS EDGES ELEMENT_EDGES FACES ELEMENT_FACES NUM_NODES NUM_ELEMENTS NUM_EDGES NUM_FACES ...
    LOCALEDGENODES LOCALFACENODES NUM_DOFS

eps_0 = 8.854e-12;
mu_0=4*pi*1e-7;
eps_r = 1.0;
mu_r = 1.0;
% Caution - do not change following TEN lines, hard-wired into other
% parts of code.
LOCALEDGENODES(1,:) = [1 2];
LOCALEDGENODES(2,:) = [1 3];
LOCALEDGENODES(3,:) = [1 4];
LOCALEDGENODES(4,:) = [2 3];
LOCALEDGENODES(5,:) = [2 4];
LOCALEDGENODES(6,:) = [3 4];
LOCALFACENODES(1,:) = [1 2 3];
LOCALFACENODES(2,:) = [1 2 4];
LOCALFACENODES(3,:) = [1 3 4];
LOCALFACENODES(4,:) = [2 3 4];

elem_order=input('Order of elements? 1 mixed 1st order, 2 mixed 2nd order (default mixed 1st): ');
if isempty(elem_order)
    elem_order = 1
end

int_mesh_flag=input('Generate mesh file internally? 1 yes, 0 no (default yes): ');
if isempty(int_mesh_flag)
    int_mesh_flag = 1
end
if ~int_mesh_flag
    mesh_filenm=input('Enter mesh file name: ','s');
    if isempty(mesh_filenm)
        mesh_filenm = 'box.msh' % box.msh is a 1x0.5x0.75 rectangular box.
        % cube.mesh cube is a 1x1x1m box.
    end
end

symmetry_flag = input('Use symmetry? 1 yes, 0 no (default no): ');
if isempty(symmetry_flag)
    symmetry_flag = 0;
end

eig_solve_flag = input('Compute only eigenvalues (to save memory)? 1 yes, 0 no (default no) ');
if isempty(eig_solve_flag)
    eig_solve_flag = 1;
end


L_x=1; % x-dimensions of resonant cavity in m
L_y=0.5; % ditto y
L_z=0.75; % ditto z
% L_x = 1
% L_y = 1
% L_z = 1


if int_mesh_flag
    N_x=4 %8; % Number of quadrilateral elements in x
    N_y=2 %4; % ditto y
    N_z=3 %6; % ditto z
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
switch elem_order
    case 1
        if symmetry_flag
            filenm=['eigdata3D_sym_CTLN_',num2str(NUM_ELEMENTS)]
        else
            filenm=['eigdata3D_CTLN_',num2str(NUM_ELEMENTS)]
        end
    case 2
        if symmetry_flag
            filenm=['eigdata3D_sym_LTQN_',num2str(NUM_ELEMENTS)]
        else
            filenm=['eigdata3D_LTQN_',num2str(NUM_ELEMENTS)]
        end
    otherwise
        error('Unimplemented element order')
end

tetramesh(TetMesh,vertices)
print('-depsc',filenm) % colour print
%print('-deps',filenm) % B&W print
print('-dpng',filenm)
pause;

edgemake3D;
facemake3D;
toler=min([L_x L_y L_z])*1e-5;
[edge_free_flag,face_free_flag]=free_dof3D_LTQN(L_x,L_y,L_z,vertices,symmetry_flag,toler);
% Renumber dofs and compute NUM_DOFS
switch elem_order
    case 1
        [dof_e1,num_free_edges] = renumber_dof(edge_free_flag);
        NUM_DOFS = num_free_edges;
    case 2
        [dof_e1,dof_e2,dof_f1,dof_f2,num_free_edges,num_free_faces] = renumber_dof_LTQN(edge_free_flag,face_free_flag);
        NUM_DOFS = 2*num_free_edges+2*num_free_faces;
    otherwise
        error('Unimplemented element order')
end
S = zeros(NUM_DOFS,NUM_DOFS);
T = zeros(NUM_DOFS,NUM_DOFS);
for ielem=1:NUM_ELEMENTS % Assemble by elements
    %for ielem=1:1
    tetnodes = ELEMENTS(ielem,:);
    switch elem_order
        case 1
            [S_elem,T_elem] = sandt3D( vertices(tetnodes(1),:),vertices(tetnodes(2),:),...
                vertices(tetnodes(3),:),vertices(tetnodes(4),:) );
%             [Stemp,Ttemp] = sandt3D_LTQN( vertices(tetnodes(1),:),vertices(tetnodes(2),:),...
%                 vertices(tetnodes(3),:),vertices(tetnodes(4),:) );
%             S_elem(1:6,1:6)=Stemp(1:6,1:6);
%             T_elem(1:6,1:6)=Ttemp(1:6,1:6);
%             % Assemble into global matrix.
            for jedge = 1:6
                jj = dof_e1(ELEMENT_EDGES(ielem,jedge));
                for kedge = 1:6
                    kk = dof_e1(ELEMENT_EDGES(ielem,kedge));
                    if jj && kk  % i.e. both free
                        S(jj,kk) = S(jj,kk)+S_elem(jedge,kedge);
                        T(jj,kk) = T(jj,kk)+T_elem(jedge,kedge);
                    end
                end
            end
        case 2
                        [S_elem,T_elem] = sandt3D_LTQN( vertices(tetnodes(1),:),vertices(tetnodes(2),:),...
                            vertices(tetnodes(3),:),vertices(tetnodes(4),:) );
%             [Stemp,Ttemp] = sandt3D_LTQN( vertices(tetnodes(1),:),vertices(tetnodes(2),:),...
%                 vertices(tetnodes(3),:),vertices(tetnodes(4),:) );
%             S_elem=zeros(20);
%             T_elem=zeros(20);
%             S_elem(1:6,1:6)=Stemp(1:6,1:6);
%             T_elem(1:6,1:6)=Ttemp(1:6,1:6);
            %             S_elem(1:12,1:12)=Stemp(1:12,1:12);
            %             T_elem(1:12,1:12)=Ttemp(1:12,1:12);
            % Assemble into global matrix
            for jedge = 1:6 % edge-edge interactions
                jj_e1 = dof_e1(ELEMENT_EDGES(ielem,jedge));
                jj_e2 = dof_e2(ELEMENT_EDGES(ielem,jedge));
                for kedge = 1:6
                    kk_e1 = dof_e1(ELEMENT_EDGES(ielem,kedge));
                    kk_e2 = dof_e2(ELEMENT_EDGES(ielem,kedge));
                    if jj_e1 && kk_e1  % i.e. both edges free (note that all edge dofs have same free/prescribed properties)
                        S(jj_e1,kk_e1) = S(jj_e1,kk_e1)+S_elem(jedge,kedge);
                        T(jj_e1,kk_e1) = T(jj_e1,kk_e1)+T_elem(jedge,kedge);
                        S(jj_e2,kk_e1) = S(jj_e2,kk_e1)+S_elem(jedge+6,kedge);
                        T(jj_e2,kk_e1) = T(jj_e2,kk_e1)+T_elem(jedge+6,kedge);
                        S(jj_e2,kk_e2) = S(jj_e2,kk_e2)+S_elem(jedge+6,kedge+6);
                        T(jj_e2,kk_e2) = T(jj_e2,kk_e2)+T_elem(jedge+6,kedge+6);
                        % Use symmetry for e1e2 interaction
                        S(kk_e1,jj_e2) = S(jj_e2,kk_e1);
                        T(kk_e1,jj_e2) = T(jj_e2,kk_e1);
                    end
                end
            end

            % Face-edge and vice-versa interactions.
            for jface = 1:4 
                jj_f1 = dof_f1(ELEMENT_FACES(ielem,jface));
                jj_f2 = dof_f2(ELEMENT_FACES(ielem,jface));
                for kedge = 1:6
                    kk_e1 = dof_e1(ELEMENT_EDGES(ielem,kedge));
                    kk_e2 = dof_e2(ELEMENT_EDGES(ielem,kedge));
                    if jj_f1 && kk_e1  % i.e. face and edge free
                        S(jj_f1,kk_e1) = S(jj_f1,kk_e1)+S_elem(jface+12,kedge);
                        T(jj_f1,kk_e1) = T(jj_f1,kk_e1)+T_elem(jface+12,kedge);
                        S(jj_f1,kk_e2) = S(jj_f1,kk_e2)+S_elem(jface+12,kedge+6);
                        T(jj_f1,kk_e2) = T(jj_f1,kk_e2)+T_elem(jface+12,kedge+6);
                        S(jj_f2,kk_e1) = S(jj_f2,kk_e1)+S_elem(jface+16,kedge);
                        T(jj_f2,kk_e1) = T(jj_f2,kk_e1)+T_elem(jface+16,kedge);
                        S(jj_f2,kk_e2) = S(jj_f2,kk_e2)+S_elem(jface+16,kedge+6);
                        T(jj_f2,kk_e2) = T(jj_f2,kk_e2)+T_elem(jface+16,kedge+6);
                        % Use symmetry for face-edge entries:
                        S(kk_e1,jj_f1) = S(jj_f1,kk_e1);
                        T(kk_e1,jj_f1) = T(jj_f1,kk_e1);
                        S(kk_e2,jj_f1) = S(jj_f1,kk_e2);
                        T(kk_e2,jj_f1) = T(jj_f1,kk_e2);
                        S(kk_e1,jj_f2) = S(jj_f2,kk_e1);
                        T(kk_e1,jj_f2) = T(jj_f2,kk_e1);
                        S(kk_e2,jj_f2) = S(jj_f2,kk_e2);
                        T(kk_e2,jj_f2) = T(jj_f2,kk_e2);
                    end
                end
            end

            % face-face interactions.
            for jface = 1:4 %
                jj_f1 = dof_f1(ELEMENT_FACES(ielem,jface));
                jj_f2 = dof_f2(ELEMENT_FACES(ielem,jface));
                for kface = 1:4
                    kk_f1 = dof_f1(ELEMENT_FACES(ielem,kface));
                    kk_f2 = dof_f2(ELEMENT_FACES(ielem,kface));
                    if jj_f1 && kk_f1  % i.e. all free
                        S(jj_f1,kk_f1) = S(jj_f1,kk_f1)+S_elem(jface+12,kface+12);
                        T(jj_f1,kk_f1) = T(jj_f1,kk_f1)+T_elem(jface+12,kface+12);
                        S(jj_f2,kk_f1) = S(jj_f2,kk_f1)+S_elem(jface+16,kface+12);
                        T(jj_f2,kk_f1) = T(jj_f2,kk_f1)+T_elem(jface+16,kface+12);
                        S(jj_f2,kk_f2) = S(jj_f2,kk_f2)+S_elem(jface+16,kface+16);
                        T(jj_f2,kk_f2) = T(jj_f2,kk_f2)+T_elem(jface+16,kface+16);
                        S(kk_f1,jj_f2) = S(jj_f2,kk_f1);
                        T(kk_f1,jj_f2) = T(jj_f2,kk_f1);
                    end
                end
            end

        otherwise
            error('Unimplemented element order')
    end
end
% Find eigensolution.

% Scale S and T before eigensolver call, to optimize memory use.
S=S/mu_r;
T=T*eps_r;
if eig_solve_flag
    eigvaluessq = eig(S,T);
    eigvalues = sort(sqrt((eigvaluessq)));
else
    [eigvectors,eigvaluessq] = eig(S,T);
    [eigvalues,sort_vec] = sort(sqrt(diag((eigvaluessq))));
end

% [eigvectors,eigvaluessq] = (eigs(sparse(S),sparse(T)));
clear S T % clear matrics before saving system, or doing further work.
[node_flag,num_free_nodes] = free_nodes3D(L_x,L_y,L_z,vertices,symmetry_flag,toler );

%Find dimension of null space (plus trivial solution)
switch elem_order
    case 1
        num_zero_eig = num_free_nodes;
    case 2
        % For 2nd order solution, dimension is number of free vertices and free mid-side nodes (identical to free edges).
        num_zero_eig = num_free_nodes+num_free_edges;
    otherwise
        error('Unimplemented element order')
end

rel_err = eig_err3D(L_x,L_y,L_z,eigvalues,min(NUM_DOFS-num_zero_eig,8),num_zero_eig)
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
