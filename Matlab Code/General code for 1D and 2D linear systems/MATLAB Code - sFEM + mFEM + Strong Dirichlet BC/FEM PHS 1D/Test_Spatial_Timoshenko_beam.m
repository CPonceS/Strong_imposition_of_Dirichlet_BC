clc, close all , clear all , format shortg

%% Physical properties
rho   = 7800;         % Density of the material
E     = 2.1e11;       % Young's Modulus
v     = 0.3;          % Poisson Ratio
G     = E/(2*(1+v));  % Shear modulus
kappa = 5/6;          % Correction factor
L     = 30/100;       % Length of the beam
b     = 1/100;        % Width of the beam
h     = 1/1000;       % Thickness of the beam
A0    = b*h;          % Cross section area
I2    = b*(h^3)/12;   % Second moment of inertia of the cross section

% Infinite-dimensional model (Spatial Timoshenko beam)
prop.dom = [0 L];
prop.M  = @(x) [rho*I2 0 0 0 0; 0 rho*I2 0 0 0; 0 0 rho*A0 0 0; 0 0 0 rho*A0 0; 0 0 0 0 rho*A0];
prop.Ke = @(x) [E*A0 0 0 0 0; 0 E*I2 0 0 0; 0 0 E*I2 0 0; 0 0 0 kappa*G*A0 0; 0 0 0 0 kappa*G*A0];
prop.F0 = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 -1 0 0 0; -1 0 0 0 0]; 
prop.F1 = [0 0 1 0 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0; 0 0 0 0 1];

% Boundary conditions
prop.Dirichlet = prop.dom(1);  % Coordinate where Dirichlet BC are applied
prop.Neumann   = prop.dom(2);  % Coordinate where Neumann BC are applied
prop.Didx = [0 0 0 0 0];       % Logical index relative to which generalized velocities are applied as control inputs.
prop.Nidx = [0 0 0 1 1];       % Logical index relative to which generalized tractions  are applied as control inputs.

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
modesSpatialTimBeam(model,method,N_modes);


%% Time simulation

tspan = [0 2]; 
dt    = 2e-4;  
time  = tspan(1):dt:tspan(2);

% Input
time_ref = 1;
Amp1 = 2; freq1 = 2; 
Amp2 = 1; freq2 = 2;
u    = @(t) [Amp1*sin(2*pi*freq1*t).*(t<time_ref);
             Amp2*cos(2*pi*freq2*t).*(t<time_ref)];

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

u0 = r(2*NN+1:3*NN,:);
v0 = r(3*NN+1:4*NN,:);
w0 = r(4*NN+1:5*NN,:);


figure('color','w')
for i=1:20:length(time)
    
    
    j = round(time_ref/(time(2)-time(1)));

    hold on, box on, grid minor, view([60 30])
    plot3(coor+u0(:,i),v0(:,i),w0(:,i),'k','linewidth',4)
    plot3(L+u0(end,1:i),v0(end,1:i),w0(end,1:i),'r','linewidth',1)
    if time(i)>time_ref
        plot3(L+u0(end,j:i),v0(end,j:i),w0(end,j:i),'b','linewidth',1)
    end
    hold off
    xlim([0 L+max(max(u0))])
    ylim([min(min(v0)) max(max(v0))])
    zlim([min(min(w0)) max(max(w0))])
    xlabel('L [m]'), ylabel('v'), zlabel('w')
    title(strcat({'Time = '},num2str(round(time(i),2)), {' [s]'}))
    pause(0)
    if time(i)<time(end)
        clf
    end
end

