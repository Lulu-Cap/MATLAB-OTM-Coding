%%
% Lucas Caparini 53547155 March 11 2020
%
% Initialize the OTM nodal (nd) and material point (mp) configurations and
% properties at time t=0.
%
% This configuration will be a long rectangular beam which will have a
% fixed displacement applied to one end and the other end held stationary
% - No gravity will be applied
% - Plane strain and Lagrangian strain measures used
% - This will be the first 3D test
% - Viscous dissipation will be implemented to damp out vibrations and
% compare to steady state solution
clear, close
tic;
%% Initialize Solver Parameters
% Grid variables
Solver.domain.dim = 3; dim = Solver.domain.dim; % Solver.dimension of problem
Solver.domain.Length = 0.35; Length = Solver.domain.Length; % Boundary definitions, x dim
Solver.domain.Width = 0.02; Width = Solver.domain.Width; % y dim
Solver.domain.Height = 0.02; Height = Solver.domain.Height; % z dim
Solver.domain.Nz = 3; Nz = Solver.domain.Nz; 
Solver.domain.Ny = Solver.domain.Nz; Ny = Solver.domain.Ny;
Solver.domain.Nx = ceil((Ny-1)*Length/Height)+1; Nx = Solver.domain.Nx; % # nodes in each direction


x = linspace(0,Length,Nx)';
y = linspace(0,Width,Ny)';
z = linspace(0,Height,Nz)';
[x,y,z] = meshgrid(x,y,z);

% LME parameters
Solver.LMEParam.gamma = 4.0; % Nondimensional spacing parameter
Solver.LMEParam.h = (Length*Width*Height / ((Nx-1)*(Ny-1)*(Nz-1)))^(1/3); % Measure of nodal spacing used for LME computation
Solver.LMEParam.beta = Solver.LMEParam.gamma/Solver.LMEParam.h.^2;
Solver.LMEParam.Tol_support = 10^-12; % Cutoff tolerance
Solver.LMEParam.R_cutoff = sqrt(-log(Solver.LMEParam.Tol_support)./Solver.LMEParam.beta); % Support radius
Solver.LMEParam.stab = 0; % Stabilization Parameters (epsilon in Weibenfels paper)


% Gravity
Solver.gravity = [0 0 0]; 

% Material Property Details
% Details of constitutive model used
Solver.Material.ConstitutiveEq = ["SolidLinearElastic", "LagStrain", " "]; % a linear elastic solid using the Lagrangian strain matrix, and plane strain assumption.
% Parameters
Solver.Material.dens0 = 1*10^8; rho = Solver.Material.dens0; % nominal material density [kg/m^3]
Solver.Material.E = 1.0*10^6; E = Solver.Material.E; % Young's Modulus [Pa]
Solver.Material.poisson = 0.4; nu = Solver.Material.poisson; % Poisson's ratio
Solver.Material.visc = 0.0; % dynamic viscosity
Solver.Material.damping = 0.0; % Damping factor. [min max] = [0 1]

% Timestepping details
K = E/(3*(1-2*nu));
Cl = sqrt(3*K*(1-nu)/(rho*(1+nu)));
dl = Height/(Ny-1);
sf = 1/10;
Solver.time.dt = sf*dl/Cl;%0.01; % size of timestep [s] computed from stability condition used by Weibenfels
Solver.time.T = 600; % duration of time integration
Solver.time.t = 0:Solver.time.dt:Solver.time.T;
Solver.time.Nt = length(Solver.time.t);
%% Nodal Quantities
nd.x0 = [x(:) y(:) z(:)]; % Nodal positions at previous timestep
nd.x_start = nd.x0; % Reference config
nd.x1 = nd.x0; % Nodal positions at current timestep
nd.f = zeros(size(nd.x0)); % Nodal forces
nd.l = nd.f; % Nodal "momentum"
nd.mass = ones(size(nd.x0,1),1); % For a lumped mass matrix

DT = delaunayTriangulation(nd.x1); % Create delaunay 
%ConHull = convexHull(DT);


% ID's of constrained nodes
Left = find(nd.x0(:,1)==0); % Left edge is stationary
Right = find(nd.x0(:,1)==Length); % Right edge prescribed displacment
nd.Dirichlet.Nodes = [Left;Right];
nd.Dirichlet.Values = nd.x0(nd.Dirichlet.Nodes,:);


% Free Nodes
nd.Free = setdiff([1:Nx*Ny*Nz]',nd.Dirichlet.Nodes);

% Solver Mass Matrix Style
nd.MassMatrixStyle = "Lumped"; % String either "Lumped" or "Distributed"
%% Material Point Quantities
mp.x1 = incenter(DT); % Material Point positions
mp.x0 = mp.x1; % mp positions from previous timestep 
mp.x_start = mp.x1; % Reference config
% Volume of the mp
mp.vol = zeros(size(mp.x0,1),1);
for ii = 1:length(mp.vol)
    s1 = nd.x1(DT.ConnectivityList(ii,1),:) - nd.x1(DT.ConnectivityList(ii,2),:);
    s2 = nd.x1(DT.ConnectivityList(ii,1),:) - nd.x1(DT.ConnectivityList(ii,3),:);
    s3 = nd.x1(DT.ConnectivityList(ii,1),:) - nd.x1(DT.ConnectivityList(ii,4),:);
    mp.vol(ii) = 1/6*abs(dot(s1,cross(s2,s3)));
end
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
% %trisurf(DT); % Plot nodal discretization
% plot3(nd.x1(:,1),nd.x1(:,2),nd.x1(:,3),'bo');hold on
% plot3(mp.x1(:,1),mp.x1(:,2),mp.x1(:,3),'*r'); % mp's
% plot3(nd.x1(nd.Dirichlet.Nodes,1),nd.x1(nd.Dirichlet.Nodes,2),nd.x1(nd.Dirichlet.Nodes,3),'kx'); % Nodes on the convex hull
% axis equal; hold off;
% legend("Nodal Triangulation","Material Points","Stationary Boundaries");
% 
% % Plot shape functions (make sure they behave correctly at a sample point)
% point = 11; % test point to visualize
% 
% figure
% neigh = Shape(point).neigh;
% plot3(mp.x1(point,1),mp.x1(point,2),mp.x1(point,3),'rx','MarkerSize',10); hold on;
% plot3(nd.x1(:,1),nd.x1(:,2),nd.x1(:,3),'bo','MarkerSize',5);
% scatter3(nd.x1(neigh,1),nd.x1(neigh,2),nd.x1(neigh,3),100,log(Shape(point).p),'filled');
% axis equal; colorbar
% xlabel('x');ylabel('y');zlabel('z');
% title("N_a(x_p) around test mp");
% 
% % Shape function gradient plots
% figure
% neigh = Shape(point).neigh;
% for ii = 1:length(neigh)
%     s = sign(Shape(point).gradp(ii,1));
%     gpx(ii) = s*log(abs(Shape(point).gradp(ii,1)));
%     s = sign(Shape(point).gradp(ii,2));
%     gpy(ii) = s*log(abs(Shape(point).gradp(ii,2)));
%     s = sign(Shape(point).gradp(ii,3));
%     gpz(ii) = s*log(abs(Shape(point).gradp(ii,3)));
% end
% subplot(3,1,1)
% plot3(mp.x1(point,1),mp.x1(point,2),mp.x1(point,3),'rx','MarkerSize',10); hold on;
% plot3(nd.x1(:,1),nd.x1(:,2),nd.x1(:,3),'bo','MarkerSize',5);
% scatter3(nd.x1(neigh,1),nd.x1(neigh,2),nd.x1(neigh,3),100,gpx,'filled');
% axis equal; colorbar
% xlabel('x');ylabel('y');zlabel('z');
% title("dN_a(x_p)/dx around test mp");
% 
% subplot(3,1,2)
% plot3(mp.x1(point,1),mp.x1(point,2),mp.x1(point,3),'rx','MarkerSize',10); hold on;
% plot3(nd.x1(:,1),nd.x1(:,2),nd.x1(:,3),'bo','MarkerSize',5);
% scatter3(nd.x1(neigh,1),nd.x1(neigh,2),nd.x1(neigh,3),100,gpy,'filled');
% axis equal; colorbar
% xlabel('x');ylabel('y');zlabel('z');
% title("dN_a(x_p)/dy around test mp");
% 
% subplot(3,1,3)
% plot3(mp.x1(point,1),mp.x1(point,2),mp.x1(point,3),'rx','MarkerSize',10); hold on;
% plot3(nd.x1(:,1),nd.x1(:,2),nd.x1(:,3),'bo','MarkerSize',5);
% scatter3(nd.x1(neigh,1),nd.x1(neigh,2),nd.x1(neigh,3),100,gpz,'filled');
% axis equal; colorbar
% xlabel('x');ylabel('y');zlabel('z');
% title("dN_a(x_p)/dz around test mp");
%% Cleanup
% Most of the variables are not desirable afterwards
clearvars E rho nu sf Cl dl ii K
clearvars DT x y s1 s2 s3 Left Right ConHull point test F
clearvars dim Nx Ny Nz Length Width Height
clearvars neigh gpx gpy gpz z s 

%close all;