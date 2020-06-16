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
Solver.domain.Width = 1; Width = Solver.domain.Width; % Boundary definitions
Solver.domain.Height = 1; Height = Solver.domain.Height;
Solver.domain.Ny = 15; Ny = Solver.domain.Ny; % 150 corresponds to 89400 mps and takes ~90s/shape function refresh, compared to Ny=15, this appears to show p_a(x_p) computation scales almost linearly with mps
Solver.domain.Nx = 2*Ny+1; Nx = Solver.domain.Nx; % # nodes in each direction


x = linspace(-Width,Width,Nx)';
y = linspace(0,Height,Ny)';
[x,y] = meshgrid(x,y);

% LME parameters
Solver.LMEParam.gamma = 0.8; % Nondimensional spacing parameter
Solver.LMEParam.h = sqrt(2*Width*Height / (Nx*Ny)); % Measure of nodal spacing used for LME computation
Solver.LMEParam.beta = Solver.LMEParam.gamma/Solver.LMEParam.h.^2;
Solver.LMEParam.Tol_support = 10^-12; % Cutoff tolerance
Solver.LMEParam.R_cutoff = sqrt(-log(Solver.LMEParam.Tol_support)./Solver.LMEParam.beta); % Support radius

% Gravity
Solver.gravity = [0 0]; 

% Timestepping details
Solver.time.dt = 0.01; % size of timestep [s]
Solver.time.T = 10; % duration of time integration
Solver.time.t = 0:Solver.time.dt:Solver.time.T;
Solver.time.Nt = length(Solver.time.t);

% Material Property Details
% Details of constitutive model used
Solver.Material.ConstitutiveEq = ["SolidLinearElastic", "LagStrain", "PlaneStrain"]; % a linear elastic solid using the Lagrangian strain matrix, and plane strain assumption.
% Parameters
Solver.Material.dens0 = 115; % nominal material density [kg/m^3]
Solver.Material.E = 25*10^6; % Young's Modulus [Pa]
Solver.Material.poisson = 0.4; % Poisson's ratio
Solver.Material.visc = 0.0; % dynamic viscosity

%% Nodal Quantities
nd.x0 = [x(:) y(:)]; % Nodal positions at previous timestep
nd.x1 = nd.x0; % Nodal positions at current timestep
nd.f = zeros(size(nd.x0)); % Nodal forces
nd.l = nd.f; % Nodal "momentum"
nd.mass = ones(size(nd.x0,1),1); % For a lumped mass matrix

DT = delaunayTriangulation(nd.x1); % Create delaunay 
ConHull = convexHull(DT);

Indent = Ny*[ceil(Nx/2)-1 ceil(Nx/2) ceil(Nx/2)+1]'; % Nodes to suffer the indentor (point force application for now)
% ID's of constrained nodes (Far-Field nodes are stationary)
    bnd_left = [1:Ny]'; % left edge
    bnd_right = [(Nx-1)*Ny+1:Nx*Ny]'; % right edge
    bnd_bott = [1:Ny:(Nx-1)*Ny+1]'; % bottom edge
nd.Dirichlet.Nodes = [unique([bnd_left;bnd_right;bnd_bott]);Indent];
nd.Dirichlet.Values = nd.x0(nd.Dirichlet.Nodes,:);

% Free Nodes
nd.Free = setdiff([1:Nx*Ny]',nd.Dirichlet.Nodes);

% Solver Mass Matrix Style
nd.MassMatrixStyle = "Lumped"; % String either "Lumped" or "Distributed"
%% Material Point Quantities
mp.x1 = incenter(DT); % Material Point positions
mp.x0 = mp.x1; % mp positions from previous timestep (note: only used when splitting!)
% Volume of the mp
    a_x = nd.x1(DT.ConnectivityList(:,1),1)-nd.x1(DT.ConnectivityList(:,2),1);
    a_y = nd.x1(DT.ConnectivityList(:,1),2)-nd.x1(DT.ConnectivityList(:,2),2);
    b_x = nd.x1(DT.ConnectivityList(:,1),1)-nd.x1(DT.ConnectivityList(:,3),1);
    b_y = nd.x1(DT.ConnectivityList(:,1),2)-nd.x1(DT.ConnectivityList(:,3),2);
mp.vol = 1/2*abs(a_x.*b_y-a_y.*b_x);
mp.dens = Solver.Material.dens0*ones(size(mp.vol)); % Density of each mp
mp.mass = mp.dens.*mp.vol; % mass of each mp. Should be a constant in initial formulation

% Deformation Measures: Note, stress and strain ones may be reduced by
% symmetry, but are not yet.
F = eye(dim);
mp.F = repmat(F(:)',length(mp.x1),1); % Deformation gradient F = [Fxx Fxy ; Fyx Fyy], F(:) = [Fxx;Fyx;Fxy;Fyy]. Initialize as identity
mp.Fdot = zeros(size(mp.F));
mp.strain = zeros(size(mp.F)); % Strain measure 
mp.strainrate = zeros(size(mp.F)); % strain rate measure
mp.stress = zeros(size(mp.F)); % Cauchy stress at the mp (or whatever stress measure used)

%% Compute Initial Shape Functions

% Version w/out neighbour search
%[nd.p,nd.gradp] = LME_Reg(nd.x1,mp.x1,Solver.gamma,Solver.h,Solver.dim); %

% Version w/ neighbour search
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
% plot(nd.x1(nd.Dirichlet.Nodes(1:end-3),1),nd.x1(nd.Dirichlet.Nodes(1:end-3),2),'kx'); % Nodes on the convex hull
% plot(nd.x1(nd.Dirichlet.Nodes(end-2:end),1),nd.x1(nd.Dirichlet.Nodes(end-2:end),2),'kv','MarkerSize',10);
% axis equal;
% legend("Nodal Triangulation","Material Points","Stationary Boundaries","Indentor Site");
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