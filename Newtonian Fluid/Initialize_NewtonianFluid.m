%%
% Lucas Caparini 53547155 May 28 2020
%
% Initialize the OTM nodal (nd) and material point (mp) configurations and
% properties at time t=0.
%
% This test will be for the Newtonian fluid material model
clear, close
tic;
%% Initialize Solver Parameters
% Grid variables
Solver.domain.dim = 2; dim = Solver.domain.dim; % Solver.dimension of problem
Solver.domain.Centre = [0,0]; Centre = Solver.domain.Centre; % Boundary definitions
Solver.domain.Radius = 0.005; Radius = Solver.domain.Radius;
Solver.domain.Nr = 5; Nr = Solver.domain.Nr; % Number of nodes in the radial direction
Solver.domain.Nc = ceil(2*pi*Nr); Nc = Solver.domain.Nc; % # nodes around the circumference
Solver.domain.Omega = 2*pi; Omega = Solver.domain.Omega; % Angular velocity of boundary

% Mesh the circle
[x,y] = MeshCircle(Centre,Radius,Nc,Nr);

% LME parameters
Solver.LMEParam.gamma = 4.0; % Nondimensional spacing parameter
Solver.LMEParam.h = sqrt(pi*Radius^2 / (Nc*(Nr-1))); h = Solver.LMEParam.h; % Measure of nodal spacing used for LME computation
Solver.LMEParam.beta = Solver.LMEParam.gamma/Solver.LMEParam.h.^2;
Solver.LMEParam.Tol_support = 10^-12; % Cutoff tolerance
Solver.LMEParam.R_cutoff = sqrt(-log(Solver.LMEParam.Tol_support)./Solver.LMEParam.beta); % Support radius
Solver.LMEParam.stab = 0; % Stabilization Parameters (epsilon in Weibenfels paper)


% Gravity
Solver.gravity = [0 0]; 

% Material Property Details
% Details of constitutive model used
Solver.Material.ConstitutiveEq = ["NewtonianFluid", "Isothermal", "Inviscid"]; % a Newtonian Fluid using isothermal pressure relation and w/Inviscid flag
% Parameters
Solver.Material.dens0 = 7.9*10^3; rho = Solver.Material.dens0; % nominal material density [kg/m^3]
Solver.Material.n = 7; % Coefficient for Tait's equation of state
Solver.Material.BulkMod = 1.0*10^8; K = Solver.Material.BulkMod; % Bulk Modulus of water [Pa] 2.2*10^9
Solver.Material.visc = 0.00642; % dynamic viscosity [Pa*s]
Solver.Material.damping = 0.0; % Damping coefficient. [min max] = [0 1]

% Timestepping details
sf = 1/20;
Umax = Omega*Radius;
Solver.time.dt = sf*h/Umax;%0.01; % size of timestep [s] computed from stability condition used by Weibenfels
% Let's try a conservative timestep estimate based on wavespeed...
c = sqrt(K/rho); % Speed of sound waves in fluid
Solver.time.dt = sf*h/c;
Solver.time.T = 4; % duration of time integration
Solver.time.t = 0:Solver.time.dt:Solver.time.T;
Solver.time.Nt = length(Solver.time.t);

%% Nodal Quantities

DT = delaunayTriangulation([x(:) y(:)]); % Create delaunay 
ConHull = convexHull(DT);

nd.x0 = DT.Points; % Nodal positions at previous timestep
nd.x_start = nd.x0; % Reference config
nd.x1 = nd.x0; % Nodal positions at current timestep
nd.f = zeros(size(nd.x0)); % Nodal forces
nd.l = nd.f; % Nodal "momentum"
nd.mass = ones(size(nd.x0,1),1); % For a lumped mass matrix

% ID's of constrained nodes
nd.Dirichlet.Nodes = ConHull; % Left edge is stationary
nd.Dirichlet.Values = nd.x0(nd.Dirichlet.Nodes(:),:);

% Free Nodes
Nnodes = size(nd.x0,1);
nd.Free = setdiff([1:Nnodes]',nd.Dirichlet.Nodes(:));

% Solver Mass Matrix Style
nd.MassMatrixStyle = "Lumped"; % String either "Lumped" or "Distributed". "Distributed" hasn't been properly implemented yet.
%% Material Point Quantities
mp.x1 = incenter(DT); % Material Point positions
mp.x0 = mp.x1; % mp positions from previous timestep 
mp.x_start = mp.x1; % Reference config
% Volume of the mp
    a_x = nd.x1(DT.ConnectivityList(:,1),1)-nd.x1(DT.ConnectivityList(:,2),1);
    a_y = nd.x1(DT.ConnectivityList(:,1),2)-nd.x1(DT.ConnectivityList(:,2),2);
    b_x = nd.x1(DT.ConnectivityList(:,1),1)-nd.x1(DT.ConnectivityList(:,3),1);
    b_y = nd.x1(DT.ConnectivityList(:,1),2)-nd.x1(DT.ConnectivityList(:,3),2);
mp.vol = 1/2*abs(a_x.*b_y-a_y.*b_x);
mp.dens = Solver.Material.dens0*ones(size(mp.vol)); % Density of each mp
mp.mass = mp.dens.*mp.vol; % mass of each mp. Should be a constant (no splitting)
mp.J = ones(size(mp.dens)); % Jacobian of deformation matrix

% Deformation Measures: Note, stress and strain ones may be reduced by
% symmetry, but are not yet.
F = eye(dim);
mp.F = repmat(F(:)',length(mp.x1),1); % Deformation gradient F = [Fxx Fxy ; Fyx Fyy], F(:) = [Fxx;Fyx;Fxy;Fyy]. Initialize as identity
mp.Fincr = mp.F; % Increment of deformation gradient.
mp.Fdot = zeros(size(mp.F));
mp.strain = zeros(size(mp.F)); % Strain measure 
mp.strainrate = zeros(size(mp.F)); % strain rate measure
mp.stress = zeros(size(mp.F)); % Cauchy stress at the mp (or whatever stress measure used)

%% Compute Initial Shape Functions

Shape = LME_Reg_Neigh(Solver,nd.x1,mp.x1); % Shape functions
RoI = NodalNeighbours(Shape,size(nd.x1,1),nd.x1,mp.x1,Solver); % Nodal Sphere of influence

toc
%% Plotting
%%%%COMMENT OUT BEFORE RUNNING OTM_main.m%%%%
% Plot to check everything looks okay

% Plot discretization
triplot(DT); % Plot nodal discretization
hold on
plot(mp.x1(:,1),mp.x1(:,2),'*r'); % mp's
plot(nd.x1(nd.Dirichlet.Nodes(:),1),nd.x1(nd.Dirichlet.Nodes(:),2),'kx'); % Nodes on the convex hull
axis equal;
legend("Nodal Triangulation","Material Points","Stationary Boundaries");

%% Cleanup
% Most of the variables are not desirable afterwards
clearvars nu sf Umax Omega h rho K
clearvars DT x y a_x a_y b_x b_y ConHull point test F
clearvars dim Nr Nc Radius Centre Nnodes 
