%% OTM Algorithm
% This document will be contain the main timestepping algorithm of the OTM
% in its most basic solid mechanics formualation. 
%
% Supplementary components, such as LME copmutations, force calculations,
% etc. will be included in separate files

%% Setup
% Initialize geometry and node/mp properties
Initialize_BeamBending;

% Material Properties
% Consider an ideal isotropic elastic material for now
% (stress and all that should already be initialized in Boussinesq.m)

% Timestepping details
Solver.time.dt = 0.000018; %0.00005; % size of timestep [s]. Note that it is constant (not adaptive)
Solver.time.T = 1; % duration of time integration
Solver.time.t = 0:Solver.time.dt:Solver.time.T;
Solver.time.Nt = length(Solver.time.t);

% Mass Matrix Method
nd.MassMatrixStyle = "Lumped"; % Simpler than "Distributed", which isn't implemented yet.

if nd.MassMatrixStyle == "Distributed"
    % DO NOT USE: Current implementation is far too slow!!!
    input("WARNING: The current implementation of the distributed/consistent mass matrix is extremely slow. \nIt is recommended you restart the simulation with a lumped mass matrix style\nPress enter to continue or Ctrl+C to exit\n");
end

%% Timestepping
% Main OTM algorithm included here.

% Set up movie
writerObj = VideoWriter('out.mp4','MPEG-4'); % Name it.
writerObj.FrameRate = 400; % How many frames per second.
open(writerObj); 

for kk = 1:Solver.time.Nt
    
    disp(sprintf('Iteration %i of %i',kk,Solver.time.Nt)); %#ok<NOPTS>
    
    % Compute nodal mass matrix, nodal momentum, nodal forces
    nd.mass = Mass_Matrix(nd,mp,Shape,RoI); % Calculates the mass matrix
    nd.l = Momentum(Solver,nd,mp,Shape); % Calculates nodal momentum
    [mp,nd.f] = Nodal_Forces(Solver,nd,mp,Shape,RoI); % Calculates nodal forces ***Not yet implemented for solids, liquids, inelasticity, damage, etc.***
    
    % Update nd.x1 (next nodal position --> explicit 1st order integration)
    nd.x0 = nd.x1;
    nd.x1 = update_nd(Solver,nd); 
    
    % Compute mp volumes, densities, Deformation gradient (+rate), shape functions
    mp.x0 = mp.x1;
    mp.x1 = move_mp(Solver,nd.x1,Shape); % Update material point position
    mp = update_properties(Solver,nd,mp,Shape); % Updates vol,F,Fdot
    Shape = LME_Reg_Neigh(Solver,nd.x1,mp.x1); % Shape functions
    RoI = NodalNeighbours(Shape,size(nd.x1,1)); % Nodal Sphere of influence
    
    % Check if splitting is required. Splitting algorithm implementation to
    % come
    
    % Record Movie
    plot(nd.x_start(:,1),nd.x_start(:,2),'ro'); 
    hold on;
    plot(nd.x1(:,1),nd.x1(:,2),'bo');
    plot(mp.x_start(:,1),mp.x_start(:,2),'m*');
    plot(mp.x1(:,1),mp.x1(:,2),'g^');
    axis equal;
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end
    hold off



end

close(writerObj); % Saves the movie.
