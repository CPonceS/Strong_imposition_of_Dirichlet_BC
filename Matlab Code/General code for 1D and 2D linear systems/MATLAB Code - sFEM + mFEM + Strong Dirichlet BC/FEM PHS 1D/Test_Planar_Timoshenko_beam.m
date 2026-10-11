clc, close all , clear all , format shortg

%% Physical properties
rho   = 7800;         % Density of the material
E     = 2.1e11;       % Young's modulus
v     = 0.3;          % Poisson ratio
G     = E/(2*(1+v));  % Shear modulus
kappa = 5/6;          % Correction factor
L     = 80/100;       % Length of the beam
b     = 3/100;        % Width of the beam
h     = 1/1000;       % Thickness of the beam
A0    = b*h;          % Cross section area
I0    = b*(h^3)/12;   % Moment of inertia of the cross section

% Infinite-dimensional model (Planar Timoshenko beam)
prop.dom = [0 L]; 
prop.M  = @(x) [rho*I0 0;0 rho*A0];
prop.Ke = @(x) [E*I0 0;0 kappa*G*A0];
prop.F0 = [0 0;-1 0]; 
prop.F1 = [1 0;0 1];

% Boundary conditions
prop.Dirichlet = prop.dom(1);  % Coordinate where Dirichlet BC are applied
prop.Neumann   = prop.dom(2);  % Coordinate where Neumann BC are applied
prop.Didx = [0 0];             % Logical index relative to which generalized
                             ... velocities are applied as control inputs.
prop.Nidx = [0 1];             % Logical index relative to which generalized 
                             ... tractions are applied as control inputs.


%% Finite-dimensional PHS model

Nelements = 30;       % Number of elements
method    = 'mFEM';   % sFEM: Hamilton's principle  |   mFEM: Hu-Washizu principle          
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
modesPlanarTimBeam(model,method,N_modes);


%% Time simulation

tspan = [0 2];  dt = 2e-4;  
time  = tspan(1):dt:tspan(2);

% Input
Amp  = 1;  freq = 2*pi*8;
u    = @(t) Amp*sin(freq*t).*(t<0.5);

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
w = r(NN+1:end,:);

figure('color','w')
for i=1:10:length(time) 

    plot(coor,w(:,i),'k','linewidth',2)
    xlabel('$x$','interpreter','latex','Fontsize',15)
    ylabel('$w_0(x,t)$','interpreter','latex','Fontsize',15)

    xlim([prop.dom(1) prop.dom(2)])
    ylim([-1 1]*(L/2))

    pause(0)
    clf
end










