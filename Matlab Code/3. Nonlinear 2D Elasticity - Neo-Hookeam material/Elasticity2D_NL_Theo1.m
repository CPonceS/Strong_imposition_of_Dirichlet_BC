clc , clear all, format shortG
close all 

%% Generate mesh 
DomainObj = load('Domain5.mat');
Hmax = 0.015;  % Intended maximum mesh edge length: Hmax > 0. Automatic -> Hmax  = [];
tag_Nodes = 'off'; tag_Edges = 'off'; tag_Elements = 'on';  % To be displayed in the Figure

meshFEM = gMesh2D(DomainObj.g,Hmax,tag_Nodes,tag_Edges,tag_Elements);
close(1);        % Comment to view the Figure with the mesh numbering
pause(0.5)

%% Physical properties
g  = 9.806;    % Acceleration of gravity
ht = 10/1000;  % thickness of the body
rho = 1000;    % density of the  material 
v = 0.4;       % Poisson's ratio material 
E = 50e3;      % Young's modulus material

% Coefficients of Generalized strain energy density function: Neo-Hookean
prop.lam_L = (E*v)/((1+v)*(1-(2*v)));
prop.mu_L  = E/(2*(1+v));
prop.h     = ht;

% Infinite-dimensional model properties
prop.M   = rho*diag([ht,ht]);
prop.b   = 0*[0;-rho*g*ht];
prop.F0L = zeros(3,2);
prop.F1L = [1 0;0 0;0 1];
prop.F2L = [0 0;0 1;1 0];

% Boundary conditions
prop.Dirichlet_Edges = [60:64 90:94]';    % Numbers of Global edges where Dirichlet BC are applied (See the Figure)
prop.Din  = prop.Dirichlet_Edges;         % Edges where Dirichlet BC are applied as control inputs.
prop.Didx = [1 1];                        % Logical index relative to which generalized velocities are applied as control inputs.
prop.Nin  = [];                           % Edges where Neumann BC are applied as control inputs.
prop.Nidx = [0 0];                        % Logical index relative to which generalized tractions are applied as control inputs.


%% Discretization (Theorem 1, P1)

n_PDE = size(prop.F0L,2);     % number of generalized displacements 
m_PDE = size(prop.F0L,1);     % number of generalized strains

NN  = size(meshFEM.coor,1);   % number of nodes
NE  = size(meshFEM.conec,1);  % number of elements

N_Omega = n_PDE*NN;
M_Omega = m_PDE*NE;

% Boundary conditions
Nodes_D = unique(meshFEM.Bconec(prop.Dirichlet_Edges,:),'stable');     
prop.Neumann_Edges = setdiff(1:size(meshFEM.Bconec,1),prop.Dirichlet_Edges')';  % Numbers of global edges where Neumann BC are applied
nDin = meshFEM.Bconec(prop.Din,:)'; Nodes_Din = unique(nDin(:)','stable');      % Nodes where are applied the Boundary Dirichlet inputs
nNin = meshFEM.Bconec(prop.Nin,:)'; Nodes_Nin = unique(nNin(:)','stable');      % Nodes where are applied the Boundary Neumann inputs

Nodes_F = setdiff(1:NN,Nodes_D);
N_free  = n_PDE*(NN-length(Nodes_D));

dof_D = reshape(Nodes_D(:)+(0:n_PDE-1)*NN,[],1);
dof_F = setdiff(1:N_Omega, dof_D)';
N_D   = length(dof_D);

dof_matrix = reshape(1:N_Omega,NN,n_PDE);
dof_Din = reshape(dof_matrix(Nodes_Din,:).*prop.Didx,1,[]);
dof_Din = dof_Din(dof_Din ~= 0); 
dof_Nin = reshape(dof_matrix(Nodes_Nin,:).*prop.Nidx,1,[]);
dof_Nin = dof_Nin(dof_Nin ~= 0);

% Sorting
Matrix_indD = kron(diag(prop.Didx),eye(length(Nodes_D)));
Ind_D = sparse(Matrix_indD(:,any(Matrix_indD)));

[~,idx_sort] = sort([dof_F' dof_D']);
MP = eye(N_Omega); % Permutation matrix: r = MP*[rr;rD] (correlative sorted)
MP = sparse(MP(idx_sort,:));


%% Constant matrices

Ae = meshFEM.area(:);
conec = meshFEM.conec;

I_M = conec(:,[1 1 1 2 2 2 3 3 3]); 
J_M = conec(:,[1 2 3 1 2 3 1 2 3]);
M_local = (Ae/12).*[2 1 1 1 2 1 1 1 2];
M_glob = sparse(I_M(:),J_M(:),M_local(:),NN,NN);

V_local = (Ae/3).*[1 1 1];
V_glob = sparse(conec(:),1,V_local(:),NN,1);

I_D = conec;
J_D = repmat((1:NE)',1,3);
Dx_local = meshFEM.shapeFcn.dNdx.*Ae;
Dy_local = meshFEM.shapeFcn.dNdy.*Ae;
Dx_glob  = sparse(I_D(:),J_D(:),Dx_local(:),NN,NE);
Dy_glob  = sparse(I_D(:),J_D(:),Dy_local(:),NN,NE);

edge_N = meshFEM.Bconec(prop.Neumann_Edges, :); % Neumann edges (Ne_N x 2)
Lb = vecnorm(meshFEM.coor(edge_N(:,2),:) - meshFEM.coor(edge_N(:,1),:),2,2);
I_Mb = edge_N(:,[1 1 2 2]);
J_Mb = edge_N(:,[1 2 1 2]);
Mb_local = (Lb/6).*[2 1 1 2];
Mb_glob  = sparse(I_Mb(:),J_Mb(:),Mb_local(:),NN,NN);

Me_scalar = sparse(1:NE,1:NE,Ae,NE,NE);

%%  Matrices of the model
I_nPDE = speye(n_PDE);  I_mPDE = speye(m_PDE);

Mn = kron(I_nPDE,M_glob);
Mi = (1/(rho*ht))*kron(I_nPDE,M_glob);
Me = kron(I_mPDE,Me_scalar);  
Ft = kron(prop.F1L',Dx_glob) + kron(prop.F2L',Dy_glob);
b  = kron(prop.b,V_glob);
BT = kron(I_nPDE, Mb_glob);

if isempty(dof_Nin);  BN   = []; else; BN   = BT(dof_F,dof_Nin); end

Mpr = Mn(dof_F,dof_F);     iMpr = sparse(eye(N_free)/Mpr);
Mpp = Mi(dof_F,dof_F);
MpD = Mn(dof_F,dof_D);
iMe = sparse(diag(1./diag(Me)));

br = b(dof_F);
bD = b(dof_D);

Bv  = iMpr*MpD;
BvD = Bv*Ind_D;

% Nonlinear terms
er_handle = @(x) NLterm(x,MP,dof_F,NN,N_free,N_D,meshFEM,prop,Ae);
U_handle  = @(x) E_elast_fun(x,MP,N_free,N_D,NN,meshFEM,prop,Ae);

%% Linear stiffness matrix 

% Constitutive matrix (Saint-Venant Kircchoff)
D11 = prop.h*(prop.lam_L + 2*prop.mu_L);
D12 = prop.h*prop.lam_L;
D33 = prop.h*prop.mu_L;

D_mat  = [D11 D12 0; D12 D11 0; 0 0 D33];
D_glob = kron(D_mat,speye(NE)); 

Dx_B = sparse(J_D(:),I_D(:),meshFEM.shapeFcn.dNdx(:),NE,NN);
Dy_B = sparse(J_D(:),I_D(:),meshFEM.shapeFcn.dNdy(:),NE,NN);

FNp_glob = kron(prop.F1L,Dx_B) + kron(prop.F2L,Dy_B);
K_r      = Ft*D_glob*FNp_glob;
Krr      = K_r(dof_F,dof_F);

JAC_mat  = sparse([zeros(N_free)     -iMpr'*Krr         zeros(N_free,N_D);
                   iMpr*Mpp           zeros(N_free)     zeros(N_free,N_D);
                   zeros(N_D,N_free)  zeros(N_D,N_free) zeros(N_D)]); 

Jx = sparse([zeros(N_free) -iMpr'  zeros(N_free,N_D); 
      iMpr           zeros(N_free,N_free+N_D); 
      zeros(N_D,2*N_free+N_D)]);

Gx = sparse([zeros(N_free,N_D)*Ind_D; -BvD; Ind_D]);

%% Time simulation

tspan = [0 1];
dt    = 2e-4;

% SV: 10 [s] app, IMR: 30 [s] app
Integrator = 'SV'; % 'SV' or 'IMR'

% Boundary velocity
AMP = 0.2; freq = 2*pi*10; tt = 0.5;
vD  = @(t) ones(size(BvD,2),1).*(AMP*sin(freq*t).*(t<tt));

% Stormer-Verlet method + IMR (quasi-Newton method)
x0 = zeros(2*N_free+N_D,1);
time = tspan(1):dt:tspan(2);   

JAC_approx = speye(2*N_free+N_D)-0.5*dt*JAC_mat;
LU_JAC = decomposition(JAC_approx,'lu');

x = x0; x_k = x0; 
E_kin = zeros(length(time),1); E_elas = E_kin; Eg = E_kin; H = E_kin;

tic
for k=1:length(time)
    progress = 100*(time(k)-time(1))/(time(end)-time(1)); clc;
    disp('Start Time Simulation')
    fprintf('Progress: %.2f%%\n',progress);
   
    % -------- Stormer-Verlet ---------------
    p_n   = x(1:N_free,1);
    rr_n  = x(N_free+1:2*N_free,1);
    rD_n  = x(2*N_free+1:2*N_free+N_D,1);
    
    er_n = er_handle(x);

    % FIRST STEP (PREDICTION)
    p_m    = (p_n + 0.5*dt*(-iMpr'*er_n + iMpr'*br));
    rr_m   = rr_n + 0.5*dt*iMpr*Mpp*p_m - 0.5*dt*BvD*vD(time(k));
    rD_m   = rD_n + 0.5*dt*Ind_D*vD(time(k));
    x_m    = [p_m;rr_m;rD_m];

    er_m = er_handle(x_m);

    % SECOND STEP (UPDATE/CORRECTION)
    if k<length(time)
        p_np   = (p_m + 0.5*dt*(-iMpr'*er_m + iMpr'*br));
        rr_np  = rr_m + 0.5*dt*iMpr*Mpp*p_np - 0.5*dt*BvD*vD(time(k+1));
        rD_np  = rD_m + 0.5*dt*Ind_D*vD(time(k+1));
        
        x = [p_np;rr_np;rD_np]; 
        x_k(:,k+1) = x;
        
        if isequal(Integrator,'IMR')
        % ------- Implicit midpoint rule ---------
            max_it = 5;   tol = 1e-6;
            x_mp   = (x_k(:,k)+x_k(:,k+1))/2; 
            p_mp   = x_mp(1:N_free,1);
            er_mp  = er_handle(x_mp);
            t_mp   = time(k)+dt/2;
    
            gxm    = x_mp - x_k(:,k) - 0.5*dt*(Jx*[Mpp*p_mp; er_mp-br; -bD] + Gx*vD(t_mp));
    
            j=0; z = x_mp;
            error_it = zeros(max_it,1);
            w = zeros(length(x),max_it);
            while norm(gxm)>tol && j<max_it
                j = j+1;
                z = z - (LU_JAC\gxm);  
                w(:,j) = z;
    
                p_mp = z(1:N_free,1);
                er_z = er_handle(z);
                gxm  = z - x_k(:,k) - 0.5*dt*(Jx*[Mpp*p_mp; er_z-br; -bD] + Gx*vD(t_mp));
                
                error_it(j) = norm(gxm);
                if isequal(j,max_it)
                    [~,idx] = min(error_it(1:j));
                    z = w(:,idx(1));
                end
            end
            x = 2*z - x_k(:,k);  x_k(:,k+1) = x;
            % % To check step error
            % norm(gxm) 
            % pause(1)
            % ----------------------------------------
        end
    end
    p_k   = x_k(1:N_free,k);
    rr_k  = x_k(N_free+1:2*N_free,k);
    rD_k  = x_k(2*N_free+1:end,k);
    
    % Evaluation of energy
    E_kin(k)  = 0.5*p_k'*Mpp*p_k;
    E_elas(k) = U_handle(x_k(:,k));
    Eg(k)     = -rr_k'*br - rD_k'*bD;
    H(k)      = E_kin(k) + E_elas(k) + Eg(k);
end

p  = x_k(1:N_free,:);
rr = x_k(N_free+1:2*N_free,:); 
rD = x_k(2*N_free+1:end,:);
r  = MP*[rr;rD];

toc


%% Energies plot

if isequal(Integrator,'SV')
    TITLE = 'St\"ormer-Verlet method';
elseif isequal(Integrator,'IMR')
    TITLE = 'Implicit midpoint rule';
end

figure('color','w')
title(TITLE,'Interpreter','latex','FontSize',14)
hold on, box on
plot(time,H,'linewidth',1)
plot(time,E_kin,'linewidth',1)
plot(time,E_elas,'linewidth',1)
% plot(time,Eg,'linewidth',1)
hold off
xlabel('Time [s]','Interpreter','latex','FontSize',16)
ylabel('Energy','Interpreter','latex','FontSize',16)

legend('$\hat{H}$','$\hat{E}_{kin}$','$\hat{E}_{elas}$','Interpreter','latex','FontSize',14,'location','northwest')
% ylim([0 0.15])
pause(2)


%% Animation

anim = 1;

if isequal(anim,1)

    u1  = r(1:end/2,:); 
    u2  = r(end/2+1:end,:);
    
    u_mag = sqrt((u1*100).^2 + (u2*100).^2);    
    cmap = custom_colormap(256,0.165*0.4);
    
    set(groot, 'defaultTextInterpreter', 'latex')
    figure('Color','w');  % set(gcf,'renderer','painters');
    Boundary_ref = 'on';
    View = [0 90];     
    
    coordx = ones(1,size(u1,2)).*meshFEM.coor(:,1)+u1;  coordy = ones(1,size(u2,2)).*meshFEM.coor(:,2)+u2;  
    for k=1:50:length(time)
    
         coorDX = meshFEM.coor(:,1) + u1(:,k);  coorDY = meshFEM.coor(:,2) + u2(:,k);  
         coorDZ = u_mag(:,k);
    
         TITLE = strcat({'Time: '},num2str(round(time(k),2)),{' [s]'});
         hold on , box on
         trimesh(meshFEM.conec,coorDX(1:size(meshFEM.coor,1)), coorDY(1:size(meshFEM.coor,1)), coorDZ(1:size(meshFEM.coor,1)),'FaceAlpha',0);
         if isequal(Boundary_ref,'on')
         trimesh(meshFEM.conec,coorDX(1:size(meshFEM.coor,1)), coorDY(1:size(meshFEM.coor,1)), coorDZ(1:size(meshFEM.coor,1)),'FaceAlpha',0);
         pdegplot(meshFEM.g)
         end
         hold off
       
         title(TITLE,'interpreter','latex','Fontsize',10)
    
         view(View);      axis off;     %colorbar     
         shading interp;  axis equal;   colormap(cmap);  caxis([min(min(u_mag)) max(max(u_mag))]*1);
    
         xlim([min(min(coordx))-0.1*max(max(coordx))  max(max(coordx))*1.1])
         ylim([min(min(coordy))-0.1*max(max(coordy))  max(max(coordy))*1.1])
        
    
        % % Colorbar settings
        % cb = colorbar('southoutside');             % horizontal and below the plot
        % cb.Label.String = '$\| r(\texttt{X},t)\|_2$ [cm]'; % title of the colorbar
        % cb.Label.Interpreter = 'latex';            % latex formatting
        % cb.Label.FontSize = 16; 
    
         pause(0)
         if k==length(time)
             break;
         end
         clf
    end
end



%% Auxiliar functions (for this specific model)

function er = NLterm(x,MP,dof_F,NN,N_free,N_D,meshFEM,prop,Ae)
    
    rr  = x(N_free+1:2*N_free,1);
    rD  = x(2*N_free+1:2*N_free+N_D,1);
    r   = MP*[rr;rD];
    
    edof_u = meshFEM.conec;       
    edof_v = meshFEM.conec + NN;  

    u_loc = r(edof_u); 
    v_loc = r(edof_v); 
    
    dNdx = meshFEM.shapeFcn.dNdx; 
    dNdy = meshFEM.shapeFcn.dNdy; 
    
    ux = sum(dNdx.*u_loc,2); 
    uy = sum(dNdy.*u_loc,2);
    vx = sum(dNdx.*v_loc,2); 
    vy = sum(dNdy.*v_loc,2);
    
    % Green-Lagrange strain tensor
    eps1 = ux + 0.5*(ux.^2 + vx.^2);
    eps2 = vy + 0.5*(uy.^2 + vy.^2);
    eps3 = uy + vx + ux.*uy + vx.*vy;

    IIIC = 4*eps1.*eps2 + 2*eps1 + 2*eps2 + 1 - eps3.^2;
    factor = prop.mu_L./(2*IIIC) - prop.lam_L/2 + prop.lam_L./(2*sqrt(IIIC));
    
    % 2PK-tensor
    e_eps1 = prop.h*(prop.mu_L - (4*eps2+2).*factor);
    e_eps2 = prop.h*(prop.mu_L - (4*eps1+2).*factor);
    e_eps3 = prop.h*(2*eps3.*factor);
    
    % Deformation gradient tensor
    F11 = 1 + ux;  
    F21 = vx;
    F12 = uy;      
    F22 = 1 + vy;
    
    % 1PK tensor
    P11 = F11.*e_eps1 + F12.*e_eps3;
    P12 = F11.*e_eps3 + F12.*e_eps2;
    P21 = F21.*e_eps1 + F22.*e_eps3;
    P22 = F21.*e_eps3 + F22.*e_eps2;
    
    Fu_local = (dNdx.*P11 + dNdy.*P12).*Ae; 
    Fv_local = (dNdx.*P21 + dNdy.*P22).*Ae; 
    
    er_global = sparse([edof_u(:);edof_v(:)],1,[Fu_local(:);Fv_local(:)],2*NN,1);
    er = er_global(dof_F);

end

function E_elas = E_elast_fun(x,MP,N_free,N_D,NN,meshFEM,prop,Ae)

    rr  = x(N_free+1:2*N_free,1);
    rD  = x(2*N_free+1:2*N_free+N_D,1);
    r   = MP*[rr;rD];
    
    edof_u = meshFEM.conec; 
    edof_v = meshFEM.conec + NN;

    u_loc = r(edof_u);
    v_loc = r(edof_v);
    
    dNdx = meshFEM.shapeFcn.dNdx; 
    dNdy = meshFEM.shapeFcn.dNdy;
    
    ux = sum(dNdx.*u_loc,2); 
    uy = sum(dNdy.*u_loc,2);
    vx = sum(dNdx.*v_loc,2); 
    vy = sum(dNdy.*v_loc,2);
    
    % Green-Lagrange strain tensor
    eps1 = ux + 0.5*(ux.^2 + vx.^2);
    eps2 = vy + 0.5*(uy.^2 + vy.^2);
    eps3 = uy + vx + ux.*uy + vx.*vy;
    
    IC   = 2*eps1 + 2*eps2 + 3;
    IIIC = 4*eps1.*eps2 + 2*eps1 + 2*eps2 + 1 - eps3.^2;
    
    Psi = 0.5*prop.h*prop.mu_L*(IC - 3 - log(IIIC)) + 0.5*prop.h*prop.lam_L*(sqrt(IIIC)-1).^2;
    E_elas = sum(Ae.*Psi);
end