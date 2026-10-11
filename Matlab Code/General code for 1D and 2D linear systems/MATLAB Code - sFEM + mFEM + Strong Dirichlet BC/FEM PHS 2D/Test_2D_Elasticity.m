clc , clear all, format shortG
close all 

%% Generate mesh 
DomainObj = load('Domain5.mat');
Hmax = 0.01;  % Intended maximum mesh edge length: Hmax > 0. Automatic -> Hmax  = [];
tag_Nodes = 'off'; tag_Edges = 'on'; tag_Elements = 'off';  % To be displayed in the Figure

meshFEM = gMesh2D_aux(DomainObj.g,Hmax,tag_Nodes,tag_Edges,tag_Elements);
% close(1);        % Comment to view the Figure with the mesh numbering
pause(0.5)

%% Physical properties
ht  = 10/1000;   % thickness of the body
rho = 1000;      % density of the  material 
v   = 0.4;       % Poisson's ratio material 
E   = 50e3;      % Young's modulus material

Cb = (E/(1-(v)^2))*[1 v 0; v 1 0; 0 0 (1-v)/2]; % Constitutive matrix for plane stress

% Infinite-dimensional model (2D Elasticity)
prop.M  = @(x1,x2) rho*diag([ht,ht]);
prop.Ke = @(x1,x2) ht*Cb;
prop.F0 = zeros(3,2);
prop.F1 = [1 0;0 0;0 1];
prop.F2 = [0 0;0 1;1 0]; 

% Boundary conditions
prop.Dirichlet_Edges = [93:100 141:148]';     % Numbers of Global edges where Dirichlet BC are applied (See the Figure)
prop.Din  = prop.Dirichlet_Edges;             % Edges where Dirichlet BC are applied as control inputs.
prop.Didx = [1 1];                            % Logical index relative to which generalized velocities are applied as control inputs.
prop.Nin  = [13 14];                          % Edges where Neumann BC are applied as control inputs.
prop.Nidx = [0 0];                            % Logical index relative to which generalized tractions  are applied as control inputs.


%% Finite-dimensional PHS model

method   = 'sFEM';   % sFEM: Hamilton's principle  |  mFEM: Hu-Washizu principle          
TAG_mFEM = 'C-D';    % 'C-C': strains continuous   |  'C-D': strains discontinuous

model = femPHS2D(prop,meshFEM,TAG_mFEM);

J  = model.(method).J;
Q  = model.(method).Q;
G  = model.(method).G;

%% Eigenvalues and "eigenmodes"

lam = eig(full(J*Q))/2/pi; % Eigenvalues  [Hz]

imag_lam = round(imag(lam),2);  
wn_fem   = sort(imag_lam(imag_lam>0))';

N_modes = 9;
modes2DElasticity(model,meshFEM,method,N_modes);



%% Time simulation

tspan = [0 1]; 
dt    = 2e-4;  
time  = tspan(1):dt:tspan(2);

% Input
AMP  = 0.2; 
freq = 2*pi*10; 
tref = 0.5;
u    = @(t) ones(size(G,2),1).*(AMP*sin(freq*t).*(t<tref));


% Initial conditions
x0 = zeros(size(J,1),1);
x  = zeros(length(x0),length(time));
H  = zeros(length(time),1);

% Implicit midpoint rule 
Id = speye(size(J,1)); 
A = (dt/2)*J*Q;
LHS_dec = decomposition(Id-A); 
RHS_mat = Id + A;

for k = 1:length(time)

    H(k) = 0.5*x(:,k)'*Q*x(:,k); 
    if k < length(time)
        u_mid = u(time(k) + dt/2);
        rhs = RHS_mat*x(:,k) + dt*(G*u_mid);
        x(:,k+1) = LHS_dec\rhs;
    end
end

figure('color','w')
plot(time,H,'r','linewidth',2)
xlabel('Time [s]')
ylabel('H(x)')

pause(1)



%% Animation

coor     = model.coor;
NN       = model.NN;
N_free   = model.N_free;
M_Omega  = model.M_Omega;
idx_sort = model.idx_sort;

if isequal(method,'mFEM'); factor = 1; else; factor = 0; end

r_unsorted = x(N_free + factor*M_Omega + 1:end,:);
r = r_unsorted(idx_sort,:);


u1  = r(1:end/2,:); 
u2  = r(end/2+1:end,:);

u_mag = sqrt((u1*100).^2 + (u2*100).^2);    
cmap  = custom_colormap(256,0.165*0.4);

set(groot, 'defaultTextInterpreter', 'latex')
% figure('Color','w');  % set(gcf,'renderer','painters');


figure('color','w','Position',[440 125 645 575])

Boundary_ref = 'on';
View = [0 90];     

coordx = ones(1,size(u1,2)).*meshFEM.coor(:,1)+u1;  coordy = ones(1,size(u2,2)).*meshFEM.coor(:,2)+u2;  
for k=1:25:length(time)

     coorDX = meshFEM.coor(:,1) + u1(:,k);  coorDY = meshFEM.coor(:,2) + u2(:,k);  
     coorDZ = u_mag(:,k);

     TITLE = strcat({'Time: '},num2str(round(time(k),2)),{' [s]'});
     hold on , box on
     trimesh(meshFEM.conec,coorDX(1:size(meshFEM.coor,1)), coorDY(1:size(meshFEM.coor,1)), coorDZ(1:size(meshFEM.coor,1)),'FaceAlpha',0);
     if isequal(Boundary_ref,'on')
     trimesh(meshFEM.conec,coorDX(1:size(meshFEM.coor,1)), coorDY(1:size(meshFEM.coor,1)), coorDZ(1:size(meshFEM.coor,1)),'FaceAlpha',0);
     meshFEM.Bplot(); 
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








