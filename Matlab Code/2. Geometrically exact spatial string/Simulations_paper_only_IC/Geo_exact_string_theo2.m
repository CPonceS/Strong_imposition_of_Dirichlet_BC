clc , close all, clear all, format shortG 
tic

%% Physical properties
rho   = 2800;         % Density of the material
E     = 2e6;          % Young's Modulus
g     = 9.806;        % Acceleration of gravity
L     = 150/100;      % Length of the string
A0    = 1e-6;         % Cross section area

prop.dom = [0 L];

prop.M =  rho*A0*eye(3);
prop.K =  E*A0;
prop.b =  [0;-rho*g*A0;0];

n_PDE  = 3; % Dimension of the generalized displacement vector
m_PDE  = 1; % Dimension of the generalized straing vector

% Boundary conditions
prop.Dirichlet = prop.dom(1);  % Coordinate where Dirichlet BC are applied
prop.Neumann   = prop.dom(2);  % Coordinate where Neumann BC are applied
prop.Didx = [0 0 0];           % Logical index relative to which generalized velocities are applied as control inputs.
prop.Nidx = [0 0 0];           % Logical index relative to which generalized tractions are applied as control inputs.

% Initial configuration
prop.r0  = @(s) s.*[1;-1;0]*sqrt(2)/2;
prop.dr0 = [1;-1;0]*sqrt(2)/2;    % ||dr0|| = 1

Nelements = 3;

%% Discretization (Theorem 2, P1-P0)

NE = Nelements;
NN = NE+1;
coor  = linspace(prop.dom(1),prop.dom(2),NN); 
conec = (0:NN-2)'+(1:2);
Le    = diff(coor)';
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
MPu = sparse([eye(NN) zeros(NN,2*NN)])*MP;      % to extract rx
MPv = sparse([zeros(NN) eye(NN) zeros(NN)])*MP; % to extract ry
MPw = sparse([zeros(NN,2*NN) eye(NN)])*MP;      % to extract rz

% Input map matrix of Neumann inputs
[is_free, loc_Nin] = ismember(dof_Nin,dof_F);
loc_Nin = loc_Nin(is_free); 
BT = sparse(loc_Nin(:),(1:length(loc_Nin))',1,N_free,length(loc_Nin));

I_NF  = sparse(eye(N_free));
I_ND  = sparse(eye(N_Dir));
Ind_D = sparse(I_ND(:,logical(prop.Didx)));

% Shape functions (P1)
x1  = coor(conec(:,1))'; x2 = coor(conec(:,2))';
NL  = @(el,x) [(x2(el)-x) (x-x1(el))]/Le(el);
dNL = [-1 1]./Le;

if not(isequal(length(Nodes_D)+length(Nodes_N),2))  
    error('All boundary points must be specified as Dirichlet or Neumann')
end

%% Constant matrices

xi_g = [-1/sqrt(3) 1/sqrt(3)];
W_g  = [1 1];
N_gauss = length(xi_g);

M_local  = zeros(NE,4);
V_local  = zeros(NE,2);

for p = 1:N_gauss
    % Shape functions 
    N_xi   = [0.5*(1-xi_g(p))  0.5*(1+xi_g(p))]; 
    % Matrices to be integrated
    M_mat  = N_xi'*N_xi; 
    % Integration
    M_local  = M_local  + W_g(p)*((Le/2).*M_mat(:)');
    V_local  = V_local  + W_g(p)*((Le/2).*N_xi);
end
% Assembly:
I_mat = [conec(:,1), conec(:,2), conec(:,1), conec(:,2)];
J_mat = [conec(:,1), conec(:,1), conec(:,2), conec(:,2)]; 
I_vec = [conec(:,1), conec(:,2)];
J_mv  = repmat((1:NE)',1,2);

M_glob  = sparse(I_mat(:),J_mat(:),M_local(:),NN,NN); % discretization of 1 (matrix) 
V_glob  = sparse(I_vec(:),1,V_local(:),NN,1);         % discretization of 1 (vector)

% Matrices of the model
I_nPDE = eye(n_PDE);

Mi  = kron(I_nPDE,(1/(rho*A0))*M_glob);
Mn  = kron(I_nPDE,M_glob);
Me  = sparse(diag(Le));   
iMe = sparse(diag(1./Le));
Ke  = E*A0*Me;
b   = kron(prop.b,V_glob);

Mpr = Mn(dof_F,dof_F);     iMpr = sparse(eye(N_free)/Mpr);
Mpp = Mi(dof_F,dof_F);
MpD = Mn(dof_F,dof_D);

br = b(dof_F);
bD = b(dof_D);

% Nonlinear terms
Fr_handle = make_handle(conec,NN,NE,Le,prop,dof_F,MP);

%% Time simulation 

tspan = [0 2]; 
dt    = 5e-3;
time  = tspan(1):dt:tspan(2);

% Initial configuration/conditions
vel0 = linspace(0,2,NN);
p0   = zeros(N_free,1);    p0((2*N_free/3)+1:end,1) = rho*A0*vel0(2:end); 
eps0 = zeros(M_Omega,1);
rr0  = zeros(N_free,1);
rD0  = zeros(N_Dir,1);
z0   = [p0;eps0;rr0];

p    = z0(1:N_free); 
eps  = z0(N_free+1:N_free+M_Omega);
rr   = z0(N_free+M_Omega+1:end);
rD   = rD0; 
z    = z0;

conf_r0 = prop.r0(coor)'; conf_r0 = conf_r0(:);
vr0     = conf_r0(dof_F);    vrD0 = conf_r0(dof_D);

Epot_ref = vr0'*br;

E_kin = zeros(1,length(time)); E_ela = E_kin; E_pot = E_kin;

opciones_fsolve = optimoptions('fsolve', 'Display', 'off'); % off, iter
for k=1:length(time)

    progress = 100*(time(k)-time(1))/(time(end)-time(1));  clc
    disp('Start Simulation')
    fprintf('Progress: %.2f%%',progress);

    Fun = @(zm) [  p(:,k) + (dt/2)*(-iMpr'*Fr_handle(zm(N_free+M_Omega+1:end),rD(:,k))'*iMe'*Ke*zm(N_free+1:N_free+M_Omega) + iMpr'*br) - zm(1:N_free);
                 eps(:,k) + (dt/2)*(iMe*Fr_handle(zm(N_free+M_Omega+1:end),rD(:,k))*iMpr*Mpp*zm(1:N_free)) - zm(N_free+1:N_free+M_Omega);
                  rr(:,k) + (dt/2)*(iMpr*Mpp*zm(1:N_free)) - zm(N_free+M_Omega+1:end)];

    zm = fsolve(Fun,z0,opciones_fsolve); z0 = zm; 
    z  = 2*zm - z;    

    p(:,k+1)   = z(1:N_free);
    eps(:,k+1) = z(N_free+1:N_free+M_Omega);
    rr(:,k+1)  = z(N_free+M_Omega+1:end);
    rD(:,k+1)  = rD(:,k); 

    E_kin(k) = 0.5*p(:,k)'*Mpp*p(:,k);
    E_ela(k) = 0.5*eps(:,k)'*Ke*eps(:,k);
    E_pot(k) = Epot_ref - rr(:,k)'*br - rD(:,k)'*bD;

end

toc

%% Energies plot

H = E_kin + E_ela + E_pot;

figure('color','w','Position',[440 278 560 308])
title(strcat({'N$^\circ$ of elements: '}, num2str(NE), {', Shape functions: P1-P0'}),'Interpreter','latex','FontSize',14)
hold on, box on
plot(time,H,'k',time,E_kin,'b',time,E_ela,'r',time,E_pot,'g','LineWidth',1.5)
hold off
xlabel('Time [s]','Interpreter','latex','FontSize',16)
ylabel('Energy','Interpreter','latex','FontSize',16)
legend('$\hat{H}$','$\hat{E}_{kin}$', '$\hat{E}_{elas}$', '$\hat{E}_{g}$','Interpreter','latex','FontSize',14,'location','northwest')

all_E = [H,E_kin,E_ela,E_pot];
ylim([min(all_E) max(all_E)]*1.1)

pause(1)

%%  Animation

anim = 1;

if isequal(anim,1)

    figure('color','w','Position',[440 380 390 318])
    for i=1:1:k
    
        x_actual(:,i) = MPu*[vr0+rr(:,i);rD(:,i)];
        y_actual(:,i) = MPv*[vr0+rr(:,i);rD(:,i)];
        z_actual(:,i) = MPw*[vr0+rr(:,i);rD(:,i)];
    
        hold on, box on, grid minor
        title(strcat({'N$^\circ$ of elements: '}, num2str(NE), {', Shape functions: P1-P0'}),'Interpreter','latex','FontSize',14)
        plot3(x_actual(:,i),z_actual(:,i),y_actual(:,i),'b-','LineWidth',2);
        plot3(x_actual(NN,1:i),z_actual(NN,1:i),y_actual(NN,1:i),'c.')
        plot3(x_actual(1,i),z_actual(1,i),y_actual(1,i),'r.','MarkerSize',20)
        text(1.45,1.4,-0.2,strcat({'Time = '},num2str(round(time(i),1)),{' [s]'}),'Interpreter','latex','FontSize',12)
    
        xlim([-1 1]*L)    , xl = xlim;
        ylim([-1 1]*L/1)  , yl = ylim;
        zlim([-1 0.0]*L)  , zl = zlim;
        
        GrayColor = [0.5 0.5 0.5]; 
        patch([xl(1) xl(2) xl(2) xl(1)], [0 0 0 0], [zl(1) zl(1) zl(2) zl(2)], GrayColor, 'FaceAlpha', 0.075, 'EdgeColor', 'none');
        hold off
    
        view([-110 35])
        pause(0)
    
        if isequal(i,k)
            break
        end
        clf
    end

end



%% Auxiliar functions

function Fr = NL_terms_fun(rr,rD,conec,NN,NE,Le,prop,dof_F,MP)
    
    r = MP*[rr;rD];
    
    r1 = r(conec); 
    r2 = r(conec+NN); 
    r3 = r(conec+2*NN);
    
    dr1 = (r1(:,2)-r1(:,1))./Le; 
    dr2 = (r2(:,2)-r2(:,1))./Le;  
    dr3 = (r3(:,2)-r3(:,1))./Le;
    
    F1NL_1 = prop.dr0(1) + dr1;  
    F1NL_2 = prop.dr0(2) + dr2; 
    F1NL_3 = prop.dr0(3) + dr3;
    
    FNp_local = [-F1NL_1 F1NL_1 -F1NL_2 F1NL_2 -F1NL_3 F1NL_3]; 
    idx_total = [conec conec+NN conec+2*NN];
    I_idx_NL  = repmat((1:NE)',1,6);
    
    F_NL = sparse(I_idx_NL(:),idx_total(:),FNp_local(:),NE,3*NN);
    Fr   = F_NL(:,dof_F);

end

function Fr_handle = make_handle(conec,NN,NE,Le,prop,dof_F,MP)
    Fr_handle = @(rr,rD) NL_terms_fun(rr,rD,conec,NN,NE,Le,prop,dof_F,MP);
end