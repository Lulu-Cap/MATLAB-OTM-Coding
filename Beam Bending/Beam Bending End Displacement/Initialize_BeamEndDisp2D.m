%%
% Lucas Caparini 53547155 May 28 2020
%
% Initialize the OTM nodal (nd) and material point (mp) configurations and
% properties at time t=0.
%
% This configuration will be a long rectangular beam which will have a
% fixed displacement applied to one end and the other end held stationary
% - No gravity will be applied
% - Plane strain and Lagrangian strain measures used
% - Viscous dissipation will be implemented to damp out vibrations and
% compare to steady state solution
clear, close
tic;
%% Initialize Solver Parameters
% Grid variables
Solver.domain.dim = 2; dim = Solver.domain.dim; % Solver.dimension of problem
Solver.domain.Width = 0.35; Width = Solver.domain.Width; % Boundary definitions
Solver.domain.Height = 0.02; Height = Solver.domain.Height;
Solver.domain.Ny = 5; Ny = Solver.domain.Ny;
Solver.domain.Nx = ceil((Ny-1)*Width/Height)+1; Nx = Solver.domain.Nx; % # nodes in each direction


x = linspace(-Width/2,Width/2,Nx)';
y = linspace(0,Height,Ny)';
[x,y] = meshgrid(x,y);

% LME parameters
Solver.LMEParam.gamma = 4.0; % Nondimensional spacing parameter
Solver.LMEParam.h = sqrt(Width*Height / ((Nx-1)*(Ny-1))); % Measure of nodal spacing used for LME computation
Solver.LMEParam.beta = Solver.LMEParam.gamma/Solver.LMEParam.h.^2;
Solver.LMEParam.Tol_support = 10^-12; % Cutoff tolerance
Solver.LMEParam.R_cutoff = sqrt(-log(Solver.LMEParam.Tol_support)./Solver.LMEParam.beta); % Support radius
Solver.LMEParam.stab = 10^4; % Stabilization Parameters (epsilon in Weibenfels paper)


% Gravity
Solver.gravity = [0 0]; 

% Material Property Details
% Details of constitutive model used
Solver.Material.ConstitutiveEq = ["SolidLinearElastic", "LagStrain", "PlaneStrain"]; % a linear elastic solid using the Lagrangian strain matrix, and plane strain assumption.
% Parameters
Solver.Material.dens0 = 1*10^8; rho = Solver.Material.dens0; % nominal material density [kg/m^3]
Solver.Material.E = 1.0*10^6; E = Solver.Material.E; % Young's Modulus [Pa]
Solver.Material.poisson = 0.4; nu = Solver.Material.poisson; % Poisson's ratio
Solver.Material.visc = 0.0; % dynamic viscosity
Solver.Material.damping = 0.0; % Damping coefficient. [min max] = [0 1]

% Timestepping details
K = E/(3*(1-2*nu));
Cl = sqrt(3*K*(1-nu)/(rho*(1+nu)));
dl = Height/Ny;
sf = 1/20;
Solver.time.dt = sf*dl/Cl;%0.01; % size of timestep [s] computed from stability condition used by Weibenfels
Solver.time.T = 300; % duration of time integration
Solver.time.t = 0:Solver.time.dt:Solver.time.T;
Solver.time.Nt = length(Solver.time.t);
%% Nodal Quantities
nd.x0 = [x(:) y(:)]; % Nodal positions at previous timestep
nd.x_start = nd.x0; % Reference config
nd.x1 = nd.x0; % Nodal positions at current timestep
nd.f = zeros(size(nd.x0)); % Nodal forces
nd.l = nd.f; % Nodal "momentum"
nd.mass = ones(size(nd.x0,1),1); % For a lumped mass matrix

DT = delaunayTriangulation(nd.x1); % Create delaunay 
ConHull = convexHull(DT);


% ID's of constrained nodes
left = 1:Ny; % Stationary left boundary nodes
right = ((Nx-1)*Ny+1):Nx*Ny; % Displaced right boundary nodes
nd.Dirichlet.Nodes = [left';right']; % Left edge is stationary
nd.Dirichlet.Values = nd.x0(nd.Dirichlet.Nodes(:),:);

% Free Nodes
nd.Free = setdiff([1:Nx*Ny]',nd.Dirichlet.Nodes(:));

% Solver Mass Matrix Style
nd.MassMatrixStyle = "Lumped"; % String either "Lumped" or "Distributed"
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
mp.mass = mp.dens.*mp.vol; % mass of each mp. Should be a constant in initial formulation
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
% % Plot to check everything looks okay
% 
% % Plot discretization
% triplot(DT); % Plot nodal discretization
% hold on
% plot(mp.x1(:,1),mp.x1(:,2),'*r'); % mp's
% plot(nd.x1(nd.Dirichlet.Nodes(:),1),nd.x1(nd.Dirichlet.Nodes(:),2),'kx'); % Nodes on the convex hull
% axis equal;
% legend("Nodal Triangulation","Material Points","Stationary Boundaries");
% 
% % Plot shape functions (make sure they behave correctly at a sample point)
% point = 11; % test point to visualize
% 
% % Plot shape functions w/neigh search
% figure
% test = zeros(size(x));
% test(Shape(point).neigh) = Shape(point).p;
% surfl(x,y,test); hold on;
% plot(mp.x1(point,1),mp.x1(point,2),'ro','MarkerSize',10);
% ylim([0 Height]); axis equal;
% 
% % Plot shape function derivatives w/neigh search:
% figure
% test = zeros(size(x));
% subplot(1,2,1);
% test(Shape(point).neigh) = Shape(point).gradp(:,1);
% surfl(x,y,test); hold on;
% plot(mp.x1(point,1),mp.x1(point,2),'ro','MarkerSize',10);
% ylim([-0.5,1.5]); title("X-derivatives");
% xlabel("x");ylabel("y");zlabel("p(x_p)");
% subplot(1,2,2);
% test(Shape(point).neigh) = Shape(point).gradp(:,2);
% surfl(x,y,test); hold on;
% plot(mp.x1(point,1),mp.x1(point,2),'ro','MarkerSize',10);
% ylim([-0.5,1.5]); title("X-derivatives");
% xlabel("x");ylabel("y");zlabel("p(x_p)");
%% Cleanup
% Most of the variables are not desirable afterwards
clearvars E rho nu sf
clearvars DT x y a_x a_y b_x b_y Indent bnd_bott bnd_left bnd_right ConHull point test F
clearvars dim Nx Ny Width Height

%close all;