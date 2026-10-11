clc , clear all, format shortG
close all 

%% Generate mesh 
DomainObj = load('Domain2.mat');
Hmax = 0.03;  % Intended maximum mesh edge length: Hmax > 0. Automatic -> Hmax  = [];
tag_Nodes = 'off'; tag_Edges = 'on'; tag_Elements = 'off';  % To be displayed in the Figure

meshFEM = gMesh2D_aux(DomainObj.g,Hmax,tag_Nodes,tag_Edges,tag_Elements);
% close(1);        % Comment to view the Figure with the mesh numbering
pause(0.5)

%% Physical properties
h   = 5/1000;       % thickness
rho = 1250;         % density of the material
T   = 15;           % Tension of the membrane
                          

% Infinite-dimensional model (2D Wave)
prop.M  = @(x1,x2) rho*h;
prop.Ke = @(x1,x2) T*eye(2);
prop.F0 = zeros(2,1);
prop.F1 = [1;0];
prop.F2 = [0;1]; 

% Boundary conditions  
prop.Dirichlet_Edges = [45:54 1:17]';         % Numbers of Global edges where Dirichlet BC are applied (See the Figure)
prop.Din  = prop.Dirichlet_Edges;        % Edges where Dirichlet BC are applied as control inputs.
prop.Didx = [1];                         % Logical index relative to which generalized velocities are applied as control inputs.
prop.Nin  = []';                         % Edges where Neumann BC are applied as control inputs.
prop.Nidx = [];                          % Logical index relative to which generalized tractions are applied as control inputs.

%% Finite-dimensional PHS model

method   = 'mFEM';   % sFEM: Hamilton's principle  |  mFEM: Hu-Washizu principle          
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
modes2DMembrane(model,meshFEM,method,N_modes);


% Analytic eigenvalues
Lx1 = max(meshFEM.coor(:,1))-min(meshFEM.coor(:,1));
Lx2 = max(meshFEM.coor(:,2))-min(meshFEM.coor(:,2));

% n: modes in x1 | m: modes in x2
% Indices (n,m) start at 1 for Fixed edges, and 0 for Free edges

% x1: Fixed-Fixed or Free-Free |      x2: Fixed-Fixed or Free-Free
% wn = @(n,m) sort(reshape((sqrt(T/(rho*h))/2)*sqrt(((n(:)/Lx1).^2) + ((m(:)'/Lx2).^2)),[],1))'; % [Hz]

% x1: Fixed-Free               |      x2: Fixed-Fixed or Free-Free
% wn = @(n,m) sort(reshape((sqrt(T/(rho*h))/2)*sqrt((((2*n(:)-1)/(2*Lx1)).^2) + ((m(:)'/Lx2).^2)),[],1))'; % [Hz]

% x1: Fixed-Free               |      x2: Fixed-Free
wn = @(n,m) sort(reshape((sqrt(T/(rho*h))/2)*sqrt((((2*n(:)-1)/(2*Lx1)).^2) + (((2*m(:)'-1)/(2*Lx2)).^2)),[],1))'; % [Hz]

wn_analytic = wn(1:5,1:5);

%% Time simulation

tspan = [0 1]; 
dt    = 5e-4;  
time  = tspan(1):dt:tspan(2);

% Input
AMP  = 0.02; 
freq = 2*pi*4; 
tref = 0.5;
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

w = r(1:NN,:);

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



