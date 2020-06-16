%%
% Lucas Caparini 53547155 March 11 2020
%
% Initialize the OTM nodal (nd) and material point (mp) configurations and
% properties at time t=0.
%
% This configuration will be a long rectangular block. Boussinesq solution
% with applied displacements or tractions due to an indentor will be
% investigated with the simple goal of validating OTM. 
clear, close
tic;
%% Initialize Solver Parameters
% Grid variables
Solver.domain.dim = 2; dim = Solver.domain.dim; % Solver.dimension of problem
Solver.domain.Width = 0.35; Width = Solver.domain.Width; % Boundary definitions
Solver.domain.Height = 0.02; Height = Solver.domain.Height;
Solver.domain.Ny = 10; Ny = Solver.domain.Ny; % 150 corresponds to 89400 mps and takes ~90s/shape function refresh, compared to Ny=15, this appears to show p_a(x_p) computation scales almost linearly with mps
Solver.domain.Nx = Ny*ceil(Width/Height); Nx = Solver.domain.Nx; % # nodes in each direction


x = linspace(-Width/2,Width/2,Nx)';
y = linspace(0,Height,Ny)';
[x,y] = meshgrid(x,y);

% LME parameters
Solver.LMEParam.gamma = 4.0; % Nondimensional spacing parameter
Solver.LMEParam.h = sqrt(2*Width*Height / (Nx*Ny)); % Measure of nodal spacing used for LME computation
Solver.LMEParam.beta = Solver.LMEParam.gamma/Solver.LMEParam.h.^2;
Solver.LMEParam.Tol_support = 10^-12; % Cutoff tolerance
Solver.LMEParam.R_cutoff = sqrt(-log(Solver.LMEParam.Tol_support)./Solver.LMEParam.beta); % Support radius

% Gravity
Solver.gravity = [0 -2]; 

% Timestepping details
Solver.time.dt = 0.01; % size of timestep [s]
Solver.time.T = 10; % duration of time integration
Solver.time.t = 0:Solver.time.dt:Solver.time.T;
Solver.time.Nt = length(Solver.time.t);

% Material Property Details
% Details of constitutive model used
Solver.Material.ConstitutiveEq = ["SolidLinearElastic", "InfStrain", "PlaneStrain"]; % a linear elastic solid using the Lagrangian strain matrix, and plane strain assumption.
% Parameters
Solver.Material.dens0 = 1000; % nominal material density [kg/m^3]
Solver.Material.E = 1.4*10^6; % Young's Modulus [Pa]
Solver.Material.poisson = 0.4; % Poisson's ratio
Solver.Material.visc = 0.0; % dynamic viscosity

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
nd.Dirichlet.Nodes = [1:Ny]'; % Left edge is stationary
nd.Dirichlet.Values = nd.x0(nd.Dirichlet.Nodes,:);

% Free Nodes
nd.Free = setdiff([1:Nx*Ny]',nd.Dirichlet.Nodes);

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
mp.Fdot = zeros(size(mp.F));
mp.strain = zeros(size(mp.F)); % Strain measure 
mp.strainrate = zeros(size(mp.F)); % strain rate measure
mp.stress = zeros(size(mp.F)); % Cauchy stress at the mp (or whatever stress measure used)

%% Compute Initial Shape Functions

Shape = LME_Reg_Neigh(Solver,nd.x1,mp.x1); % Shape functions
RoI = NodalNeighbours(Shape,size(nd.x1,1)); % Nodal Sphere of influence

toc
%% Plotting
% % Plot to check everything looks okay
% 
% % Plot discretization
% triplot(DT); % Plot nodal discretization
% hold on
% plot(mp.x1(:,1),mp.x1(:,2),'*r'); % mp's
% plot(nd.x1(nd.Dirichlet.Nodes,1),nd.x1(nd.Dirichlet.Nodes,2),'kx'); % Nodes on the convex hull
% axis equal;
% legend("Nodal Triangulation","Material Points","Stationary Boundaries");
% 
% % Plot shape functions (make sure they behave correctly at a sample point)
% point = 11; % test point to visualize
% % figure
% % surfl(x,y, reshape(nd.p(:,point),size(x)));
% % ylim([-0.5,1.5]);
% 
% % Plot shape functions w/neigh search
% figure
% test = zeros(size(x));
% test(Shape(point).neigh) = Shape(point).p;
% surfl(x,y,test); hold on;
% plot(mp.x1(point,1),mp.x1(point,2),'ro','MarkerSize',10);
% ylim([-0.5,1.5]);
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
clearvars DT x y a_x a_y b_x b_y Indent bnd_bott bnd_left bnd_right ConHull point test F
clearvars dim Nx Ny Width Height

%close all;