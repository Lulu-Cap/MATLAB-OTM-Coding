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
% Solver.time.dt = 0.001; %0.00005; % size of timestep [s]. Note that it is constant (not adaptive)
% Solver.time.T = 1; % duration of time integration
% Solver.time.t = 0:Solver.time.dt:Solver.time.T;
% Solver.time.Nt = length(Solver.time.t);

% Mass Matrix Method
nd.MassMatrixStyle = "Lumped"; % Simpler than "Distributed", which isn't implemented yet.

if nd.MassMatrixStyle == "Distributed"
    % DO NOT USE: Current implementation is far too slow!!!
    input("WARNING: The current implementation of the distributed/consistent mass matrix is extremely slow. \nIt is recommended you restart the simulation with a lumped mass matrix style\nPress enter to continue or Ctrl+C to exit\n");
end

%% Timestepping
% Main OTM algorithm included here.

%% Set up movie
writerObj1 = VideoWriter('out.mp4','MPEG-4'); % Name it.
writerObj1.FrameRate = 1000; % How many frames per second.
open(writerObj1); 
figure('units','normalized','outerposition',[0.5 0.5 1 1])

Ny = Solver.domain.Ny;
dx = round(10^6*(nd.x_start(Ny+1,1)-nd.x_start(1,1)))/10^6;
xv = [min(nd.x_start(:,1)):dx:max(nd.x_start(:,1))]';
yv = [min(nd.x_start(:,2)):dx:max(nd.x_start(:,2))]';
[X,Y] = meshgrid(xv,yv);
clear dx xv yv Ny

%%
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
    if mod(kk,100)==0     
        %% Displacement plot
        subplot(3,3,[1,2,3])
        plot(nd.x_start(:,1),nd.x_start(:,2),'ro');
        hold on;
        u_nd = nd.x1 - nd.x_start; amp = 2; % Amplify the displacment
        plot(nd.x_start(:,1)+amp*u_nd(:,1),nd.x_start(:,2)+amp*u_nd(:,2),'bo');
        plot(mp.x_start(:,1),mp.x_start(:,2),'m*');
        u_mp = mp.x1 - mp.x_start;
        plot(mp.x_start(:,1)+amp*u_mp(:,1),mp.x_start(:,2)+amp*u_mp(:,2),'g^');
        axis equal;
        legend('node start','mp start','node','mp');
        hold off;
        %% Stress Plots
        % Calculate stress values at nodes
        Fxx = scatteredInterpolant(mp.x_start(:,:),mp.stress(:,1),'linear','linear');
        Fxy = scatteredInterpolant(mp.x_start(:,:),mp.stress(:,2),'linear','linear');
        Fyy = scatteredInterpolant(mp.x_start(:,:),mp.stress(:,4),'linear','linear');
        
        Nd_stress_xx = Fxx(nd.x_start(:,:));
        Nd_stress_xy = Fxy(nd.x_start(:,:));
        Nd_stress_yy = Fyy(nd.x_start(:,:));
        
        Zxx = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_stress_xx,X,Y);
        Zxy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_stress_xy,X,Y);
        Zyy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_stress_yy,X,Y);
        
        % Create the plots
        subplot(3,3,4)
        contourf(X,Y,Zxx); axis equal;
        title("\sigma_{xx}");
        colorbar('southoutside');
        
        subplot(3,3,5)
        contourf(X,Y,Zxy); axis equal;
        title("\sigma_{xy}");
        colorbar('southoutside');
        
        subplot(3,3,6)
        contourf(X,Y,Zyy); axis equal;
        title("\sigma_{yy}");
        colorbar('southoutside');
        %% Strain Plots   
        % Calculate strain at nodes
        Fxx = scatteredInterpolant(mp.x_start(:,:),mp.strain(:,1),'linear','linear');
        Fxy = scatteredInterpolant(mp.x_start(:,:),mp.strain(:,2),'linear','linear');
        Fyy = scatteredInterpolant(mp.x_start(:,:),mp.strain(:,4),'linear','linear');
        
        Nd_strain_xx = Fxx(nd.x_start(:,:));
        Nd_strain_xy = Fxy(nd.x_start(:,:));
        Nd_strain_yy = Fyy(nd.x_start(:,:));
        
        Zxx = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_strain_xx,X,Y);
        Zxy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_strain_xy,X,Y);
        Zyy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_strain_yy,X,Y);

        subplot(3,3,7)
        contourf(X,Y,Zxx); axis equal;
        title("\epsilon_{xx}");
        colorbar('southoutside');
        
        subplot(3,3,8)
        contourf(X,Y,Zxy); axis equal;
        title("\epsilon_{xy}");
        colorbar('southoutside');
        
        subplot(3,3,9)
        contourf(X,Y,Zyy); axis equal;
        title("\epsilon_{yy}");
        colorbar('southoutside');
        
        frame = getframe(gcf);
        writeVideo(writerObj1, frame);
    end


end

close(writerObj1); % Saves the movie.










