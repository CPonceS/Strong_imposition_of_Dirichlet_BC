clc, close all, clear all, format shortG
tic

%% Physical properties
rho   = 7800;         % Density of the material
E     = 2.1e11;       % Young's Modulus
v     = 0.3;          % Poisson Ratio
G     = E/(2*(1+v));  % Shear modulus
g     = 9.806;        % Acceleration of gravity
kappa = 5/6;          % Correction factor
L     = 50/100;       % Length of the beam
b     = 3/100;        % Width of the beam
h     = 1/1000;       % Thickness of the beam
A0    = b*h;          % Cross section area
I0    = b*(h^3)/12;   % Second moment of inertia of the cross section

prop.dom = [0 L]; 

prop.M  = [rho*I0 0;0 rho*A0];
prop.K  = [E*I0 0;0 kappa*G*A0];
prop.b  = [0;-rho*g*A0]*0;
prop.F0 = [0 0;-1 0]; 
prop.F1 = [1 0;0 1];

% Boundary conditions
prop.Dirichlet = prop.dom(1);  % Coordinate where Dirichlet BC are applied
prop.Neumann   = prop.dom(2);  % Coordinate where Neumann BC are applied
prop.Didx = [0 0];             % Logical index relative to which generalized velocities are applied as control inputs.
prop.Nidx = [0 1];             % Logical index relative to which generalized tractions are applied as control inputs.


Nelements = 3;

%% Discretization (Theorem 1, P2)

n_PDE = size(prop.F0,2);  % number of generalized displacements 
m_PDE = size(prop.F0,1);  % number of generalized strains

NE = Nelements;
NN = 2*NE + 1;          

coor  = linspace(prop.dom(1),prop.dom(2),NN); 
conec = [(1:2:2*NE-1)' (2:2:2*NE)' (3:2:2*NE+1)'];
Le    = (coor(conec(:,3)) - coor(conec(:,1)))'; 

Nodes_D = [isequal(prop.Dirichlet,prop.dom(1))*1 isequal(prop.Dirichlet,prop.dom(2))*NN isequal(prop.Dirichlet,[prop.dom(1) prop.dom(2)])*[1 NN] isequal(prop.Dirichlet,[prop.dom(1); prop.dom(2)])*[1 NN]];
Nodes_D = Nodes_D(Nodes_D>0);
Nodes_N = [isequal(prop.Neumann,prop.dom(1))*1 isequal(prop.Neumann,prop.dom(2))*NN isequal(prop.Neumann,[prop.dom(1) prop.dom(2)])*[1 NN] isequal(prop.Neumann,[prop.dom(1); prop.dom(2)])*[1 NN]];
Nodes_N = Nodes_N(Nodes_N>0);
Nodes_F = setdiff(1:NN,Nodes_D);

dof_D = reshape(Nodes_D(:)+(0:n_PDE-1)*NN,[],1);
dof_N = reshape(Nodes_N(:)+(0:n_PDE-1)*NN,[],1);
dof_F = setdiff(1:n_PDE*NN,dof_D);
dof_matrix = reshape(1:n_PDE*NN,NN,n_PDE);
dof_Din = reshape(dof_matrix(Nodes_D,:).*prop.Didx,1,n_PDE*length(Nodes_D)); dof_Din = dof_Din(dof_Din~=0);
dof_Nin = reshape(dof_matrix(Nodes_N,:).*prop.Nidx,1,n_PDE*length(Nodes_N)); dof_Nin = dof_Nin(dof_Nin~=0);
N_Omega = n_PDE*NN;
M_Omega = m_PDE*NE;
N_free  = n_PDE*(NN-length(Nodes_D));
N_Dir   = length(dof_D);

% Permutation matrix: r = MP*[rr;rD] (correlative sorted)
[~,idx_sort] = sort([dof_F dof_D']); MP = eye(N_Omega);
MP  = sparse(MP(idx_sort,:));

% Input map matrix of Neumann inputs
[is_free,loc_Nin] = ismember(dof_Nin,dof_F);
loc_Nin = loc_Nin(is_free); 
dof_Nin_active = dof_Nin(is_free);
BT = sparse(loc_Nin(:),(1:length(loc_Nin))',1,N_free,length(loc_Nin));

I_NF  = sparse(eye(N_free));
I_ND  = sparse(eye(N_Dir));
Ind_D = sparse(I_ND(:,logical(prop.Didx)));

if not(isequal(length(Nodes_D)+length(Nodes_N),2))  
    error('All boundary points must be specified as Dirichlet or Neumann')
end

%% Constant matrices 

xi_g = [-sqrt(3/5) 0 sqrt(3/5)];
W_g  = [5/9 8/9 5/9];
N_gauss = length(xi_g);

M_local = zeros(NE,9);
V_local = zeros(NE,3);

for p = 1:N_gauss
    xi = xi_g(p);
    % Shape functions
    N_xi   = [0.5*(xi^2-xi)   1-xi^2   0.5*(xi^2+xi)]; 
    % Matrices to be integrated
    M_mat  = N_xi'*N_xi; 
    % Integration
    M_local  = M_local + W_g(p)*((Le/2).*M_mat(:)');
    V_local  = V_local + W_g(p)*((Le/2).*N_xi);
end
% Assembly:
I_mat = [conec(:,1), conec(:,2), conec(:,3), conec(:,1), conec(:,2), conec(:,3), conec(:,1), conec(:,2), conec(:,3)];
J_mat = [conec(:,1), conec(:,1), conec(:,1), conec(:,2), conec(:,2), conec(:,2), conec(:,3), conec(:,3), conec(:,3)]; 
I_vec = [conec(:,1), conec(:,2), conec(:,3)];

M_glob = sparse(I_mat(:),J_mat(:),M_local(:),NN,NN); 
V_glob = sparse(I_vec(:),1,V_local(:),NN,1);         

%% Stiffness matrix 

k_bend_base  = [7 -8 1 -8 16 -8 1 -8 7];
k_bend       = (E*I0./(3*Le)).*k_bend_base;                  
k_sh_dw      = (kappa*G*A0./(3*Le)).*k_bend_base;        

k_mass_base  = [4 2 -1 2 16 2 -1 2 4];
k_sh_psi     = (kappa*G*A0.*Le/30).*k_mass_base; 

k_cross_base = [-3 -4 1 4 0 -4 -1 4 3];
k_sh_psiw    = (kappa*G*A0/6*ones(NE,1)).*k_cross_base; 
k_sh_wpsi    = (kappa*G*A0/6*ones(NE,1)).*[-3 4 -1 -4 0 4 1 -4 3];

% Assembly
I_psi = I_mat; I_w = I_mat + NN; 
J_psi = J_mat; J_w = J_mat + NN;      
I_Kr  = [I_psi(:); I_w(:); I_psi(:); I_w(:)];
J_Kr  = [J_psi(:); J_w(:); J_w(:); J_psi(:)];
V_Kr  = [k_bend(:) + k_sh_psi(:); k_sh_dw(:); k_sh_psiw(:); k_sh_wpsi(:)];
       
Kr = sparse(I_Kr,J_Kr,V_Kr,N_Omega,N_Omega);

%% Matrices of the model
I_nPDE = eye(n_PDE);

Mn = kron(I_nPDE,M_glob);
b  = kron(prop.b,V_glob);

Mpr = Mn(dof_F,dof_F);    iMpr = sparse(eye(N_free)/Mpr);
Krr = Kr(dof_F,dof_F);
br  = b(dof_F);
bD  = b(dof_D);

%% Static solution (Theorem 1, P2)

Force  = 1; 
factor = 0:1:5;
rr = zeros(N_free,length(factor)); 
rD = zeros(N_Dir,1);

for i = 1:length(factor)
    P_star = Force*factor(i);
    F_ext  = iMpr'*BT*P_star;
    rr_sol = (iMpr'*Krr)\F_ext;
    rr(:,i) = rr_sol;
end

r = MP*[rr;repmat(rD,1,length(factor))];
psi_static = r(1:NN,:);       
w_static   = r(NN+1:2*NN,:);

%% Static configuration 

figure('color','w','Position',[440 278 560 308])
hold on, box on
plot(coor,w_static,'LineWidth',1.5)
text(0.29,0.37,'Full integration','Interpreter','latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
grid on
xlabel('$X$','Interpreter','latex','FontSize',16)
ylabel('Beam configurations','Interpreter','latex','FontSize',16)
hold on

% Exact analytical solution
XX = linspace(0,L,100);
wx = @(P) (P/(E*I0))*(((L*XX.^2)/2)-(XX.^3/6)) + ((P*XX)/(kappa*G*A0));
title(strcat({'N$^\circ$ of elements: '}, num2str(NE), {', Shape functions: P2'}),'Interpreter','latex','FontSize',14)
for i=1:length(factor)  
    P = (Force*factor(i));             
    plot(XX,wx(P),'k--','LineWidth',1.5)
end

legend('$P^\star = 0$','$P^\star = 1$','$P^\star = 2$','$P^\star = 3$','$P^\star = 4$','$P^\star = 5$','Exact','Interpreter','latex','FontSize',14,'location','northwest')
xlim([0 L]*1.1)
ylim([-0.05 0.45])