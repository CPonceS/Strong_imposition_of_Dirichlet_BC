clc , clear all, format shortG
close all 

%% Generate mesh 
DomainObj = load('Domain3.mat');
Hmax = 0.08;  % "Maximum" edge length: Hmax > 0. Automatic -> Hmax  = [];
tag_Nodes = 'off'; tag_Edges = 'on'; tag_Elements = 'off'; 

meshFEM = gMesh2D(DomainObj.g,Hmax,tag_Nodes,tag_Edges,tag_Elements);
% close(1);        % Comment to view the Figure with the mesh numbering
pause(0.5)

%% Physical properties
h   = 1/1000;       % thickness
rho = 7730;         % density of the material
I2  = (h^3)/12;     % second moment of inertia of the thickness
E   = 2.05e11;      % Young's modulus
v   = 0.3;          % Poisson's ratio
G   = E/(2*(1+v));  % Shear modulus

Cb = (E/(1-(v^2)))*[1 v 0; v 1 0; 0 0 (1-v)/2];  % Constitutive matrix for plane stress
Cs = 5/6*[G 0;0 G];                              % Constitutive matrix for shearing stress

% Infinite-dimensional model (Mindlin Plate)
prop.M  = @(x1,x2) rho*diag([I2,I2,h]);
prop.Ke = @(x1,x2) blkdiag(I2*Cb,h*Cs);
prop.F0 = [zeros(3,2) zeros(3,1);-eye(2) zeros(2,1)];
prop.F1 = [1 0 0;0 0 0;0 1 0;0 0 1;0 0 0];
prop.F2 = [0 0 0;0 1 0;1 0 0;0 0 0;0 0 1]; 

% Boundary conditions  
prop.Dirichlet_Edges = (19:23)';   % Numbers of Global edges where Dirichlet BC are applied (See the Figure)
prop.Din  = [];                    % Edges where Dirichlet BC are applied as control inputs.
prop.Didx = [0 0 1];               % Logical index relative to which generalized velocities are applied as control inputs.
prop.Nin  = [1 2]';                % Edges where Neumann BC are applied as control inputs.
prop.Nidx = [0 0 1];               % Logical index relative to which generalized tractions are applied as control inputs.

%% Finite-dimensional PHS model

method   = 'mFEM';   %  sFEM: Hamilton's principle |   mFEM: Hu-Washizu principle        
TAG_mFEM = 'C-C';    % 'C-C': strains continuous   |  'C-D': strains discontinuous

model = femPHS2D(prop,meshFEM,TAG_mFEM);

J  = model.(method).J;
Q  = model.(method).Q;
G  = model.(method).G;

%% Eigenvalues and "eigenmodes"

lam = eig(full(J*Q))/2/pi; % Eigenvalues  [Hz]

imag_lam = round(imag(lam),2);  
wn_fem   = sort(imag_lam(imag_lam>0))';

N_modes = 9;
modesMindlinPlate(model,meshFEM,method,N_modes);

%% Time simulation

tspan = [0 1];   dt = 2e-4;  
time  = tspan(1):dt:tspan(2);

% Input
AMP  = 10*0.01;   freq = 2*pi*15;    tref = 0.25;
u    = @(t) ones(size(G,2),1).*(AMP*freq*sin(freq*t).*(t<tref));

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

w = r(2*NN+1:end,:);

u_mag = w;    
cmap  = custom_colormap(256,0.165*0.4);

set(groot, 'defaultTextInterpreter', 'latex')


figure('color','w','Position',[440 125 645 575])

Boundary_ref = 'on';
View = [35 30];     
  
coorDX = meshFEM.coor(:,1);  coorDY = meshFEM.coor(:,2); 
for k=1:25:length(time)
    
     coorDZ = w(:,k);

     TITLE = strcat({'Time: '},num2str(round(time(k),2)),{' [s]'});
     hold on , box on
     trisurf(meshFEM.conec,coorDX(1:size(meshFEM.coor,1)), coorDY(1:size(meshFEM.coor,1)), coorDZ(1:size(meshFEM.coor,1)),'FaceAlpha',1);
     if isequal(Boundary_ref,'on')
     trimesh(meshFEM.conec,coorDX(1:size(meshFEM.coor,1)), coorDY(1:size(meshFEM.coor,1)), coorDZ(1:size(meshFEM.coor,1)),'FaceAlpha',0);
     meshFEM.Bplot();
     end
     hold off

     title(TITLE,'interpreter','latex','Fontsize',10)

     view(View);      axis off;     %colorbar     
     shading interp;  axis equal;   colormap(cmap);  caxis([min(min(u_mag)) max(max(u_mag))]*1);

     xlim([min(min(coorDX))-0.1*max(max(coorDX))  max(max(coorDX))*1.1])
     ylim([min(min(coorDY))-0.1*max(max(coorDY))  max(max(coorDY))*1.1])
     zlim([min(min(w))-0.1*max(max(w))            max(max(w))*1.1])

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




