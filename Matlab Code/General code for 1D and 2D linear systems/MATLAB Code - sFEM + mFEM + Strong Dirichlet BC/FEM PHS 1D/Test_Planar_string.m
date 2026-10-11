clc, close all , clear all , format shortg

%% Physical properties
rho = 2800;      % Density of the material
T   = 2;         % Tension of the string
L   = 1.5;       % Length of the string
A0  = 1e-6;      % Cross section area

% Infinite-dimensional model properties (Planar String)
prop.dom = [0 L];
prop.M  = @(x) rho*A0;
prop.Ke = @(x) T;
prop.F0 = 0; 
prop.F1 = 1;

% Boundary conditions
prop.Dirichlet = prop.dom;  % Coordinate where Dirichlet BC are applied
prop.Neumann   = [];        % Coordinate where Neumann BC are applied
prop.Didx = 1;              % Logic index relative to which generalized velocities are applied as control inputs.
prop.Nidx = 1;              % Logic index relative to which generalized tractions  are applied as control inputs.


%% Finite-dimensional PHS model

Nelements = 20;       % Number of elements
method    = 'sFEM';   % sFEM: Hamilton's principle  |   mFEM: Hu-Washizu principle          
TAG_mFEM  = 'C-D';    % 'C-C': strains continuous   |  'C-D': strains discontinuous

model = femPHS1D(prop,Nelements,TAG_mFEM);

J  = model.(method).J;
Q  = model.(method).Q;
G  = model.(method).G;

%% Eigenvalues and "eigenmodes"

lam = eig(full(J*Q))/2/pi; % Eigenvalues  [Hz]

imag_lam = round(imag(lam),2);  
wn_fem   = sort(imag_lam(imag_lam>0))';

N_modes  = 6; 
modesPlanarString(model,method,N_modes)

% Analytic eigenvalues
% wn  = @(n) (n/(2*L))*sqrt(T/(rho*A0));         % Fixed-Fixed, n = 1,2,3,... [Hz]
% wn  = @(n) ((2*n-1)/(4*L))*sqrt(T/(rho*A0));   % Mixed,       n = 1,2,3,... [Hz]
% wn  = @(n) (n/(2*L))*sqrt(T/(rho*A0));         % Free-Free,   n = 0,1,2,... [Hz]

%% Time simulation

tspan = [0 1];   dt = 2e-4;  
time  = tspan(1):dt:tspan(2);

% Input
Amp  = 0.01;  freq = 2*pi*15;
u    = @(t) [-freq*Amp*sin(freq*t);0].*(t<0.5);

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
N_free   = model.N_free;
M_Omega  = model.M_Omega;
idx_sort = model.idx_sort;

if isequal(method,'mFEM'); factor = 1; else; factor = 0; end

r_unsorted = x(N_free + factor*M_Omega + 1:end,:);
r = r_unsorted(idx_sort,:);

figure('color','w')
for i=1:10:length(time) 
    
    TITLE = strcat({'Time: '},num2str(round(time(i),2)),{' [s]'});

    plot(coor,r(:,i),'k','linewidth',2)
    xlabel('$x$','interpreter','latex','Fontsize',15)
    ylabel('$w_0(x,t)$','interpreter','latex','Fontsize',15)
    title(TITLE,'interpreter','latex','Fontsize',10)

    xlim([prop.dom(1) prop.dom(2)])
    ylim([-1 1]*(L/2))

    pause(0)
    clf
end




