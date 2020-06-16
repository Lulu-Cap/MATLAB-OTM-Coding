function [mp,force] = LinearElasticForces(Solver,nd,mp,Shape,RoI)
% Computes the force on a set of nodes based on a linear elastic model. 
%   Solver: Struct containing solver details
%   nd: struct containing nodal details
%   mp: struct containing material point details
%   Shape: struct array containing shape function evaluation and details
%   for the neighbourhood
    dim = Solver.domain.dim;
    dt = Solver.time.dt;
    StrainType = Solver.Material.ConstitutiveEq(2);
    StressType = Solver.Material.ConstitutiveEq(3);
    
    % Compute the strain
    if StrainType == "LagStrain"
        % Lagrangian strain measure
        [mp.strain,mp.strainrate] = LagStrain(mp,dim,dt);
    elseif StrainType == "InfStrain"
        % Infinitesimal strain measure
        [mp.strain, mp.strainrate] = InfStrain(mp,dim,dt);
    else
        % Not recognized strain measure
        input("Not a recognized strain measure. Please exit and select a known measure (LagStrain, InfStrain, etc)");
        return;
    end
    
    % Compute linear elastic stress 
    if dim == 1
        mp.stress = Solver.Material.E * mp.strain;
    elseif dim ==2
        if StressType == "PlaneStrain"
            mp.stress = StressPlaneStrain(mp,Solver); % Plane strain formulation
        elseif StressType == "PlaneStress"
            mp.stress = StressPlaneStress(mp,Solver); % Plane stress formulation
        elseif StressType == "W"
            mp.stress = StressWeibenfels(mp,Solver); % the formulation from Weibenfels
        else
            mp.stress = Stress2D(mp,Solver); % Simply ignores the 3rd dimension
        end
    elseif dim == 3
        mp.stress = Stress3D(mp,Solver);
    else
        % Not a recognized mass matrix type
        disp("ERROR: Invalid dimension of problem");
        return;
    end
    
    % Compute forces on the nodes
    force = zeros(size(nd.x1,1),dim);
    grav = Solver.gravity;
    for node = 1:size(nd.x1,1)
        for jj = 1:RoI(node).numneigh
            matpnt = RoI(node).neigh(jj);
            index = RoI(node).index(jj); % The index of the node within the mp's neighbour list
            force(node,:) = force(node,:) + ...
                mp.mass(matpnt)*grav*Shape(matpnt).p(index) ... body force contribution
                - mp.vol(matpnt)*TensMult(mp.stress(matpnt,:),Shape(matpnt).gradp(index,:),dim); % Stress forces
        end
    end
       
end


%%%%%%%%%%%% Minor Functions %%%%%%%%%%%%%%%%%%

% Strain Functions
function [strain, strainrate] = LagStrain(mp,dim,dt)
% Computes the Lagrangian Strain tensor and time derivative for each mp
% E = 1/2*(F^T * F - I)
% Hard-coded to avoid looping. May be unnecessary though
    strain = zeros(size(mp.x1,1),dim^2);
    
    if dim == 1
        strain = (mp.F.^2 - 1) / 2;
    elseif dim == 2
        strain = 1/2 * [mp.F(:,1).^2 + mp.F(:,2).^2 - 1, ...
            mp.F(:,1).*mp.F(:,3) + mp.F(:,2).*mp.F(:,4), ...
            mp.F(:,1).*mp.F(:,3) + mp.F(:,2).*mp.F(:,4), ...
            mp.F(:,3).^2 + mp.F(:,4).^2 - 1];
    elseif dim == 3
        strain = 1/2 * [mp.F(:,1).^2 + mp.F(:,2).^2 + mp.F(:,3).^2- 1, ...
            mp.F(:,1).*mp.F(:,4) + mp.F(:,2).*mp.F(:,5) + mp.F(:,3).*mp.F(:,6), ...
            mp.F(:,1).*mp.F(:,7) + mp.F(:,2).*mp.F(:,8) + mp.F(:,3).*mp.F(:,9), ...
            mp.F(:,1).*mp.F(:,4) + mp.F(:,2).*mp.F(:,5) + mp.F(:,3).*mp.F(:,6), ...
            mp.F(:,4).^2 + mp.F(:,5).^2 + mp.F(:,6).^2- 1, ...
            mp.F(:,4).*mp.F(:,7) + mp.F(:,5).*mp.F(:,8) + mp.F(:,6).*mp.F(:,9), ...
            mp.F(:,1).*mp.F(:,7) + mp.F(:,2).*mp.F(:,8) + mp.F(:,3).*mp.F(:,9), ...
            mp.F(:,4).*mp.F(:,7) + mp.F(:,5).*mp.F(:,8) + mp.F(:,6).*mp.F(:,9), ...
            mp.F(:,7).^2 + mp.F(:,8).^2 + mp.F(:,9).^2- 1];
    end
    
    strainrate = (strain - mp.strain)/dt; % Strain rate update

end

function [strain, strainrate] = InfStrain(mp,dim, dt)
% Computes the infinitesimal strain tensor and time derivative for each mp
% E = 1/2(F + F^T) - I
    strain = zeros(size(mp.x1,1),dim^2);
    
    if dim == 1
        strain = mp.F - 1;
    elseif dim == 2
        strain = [mp.F(:,1) - 1, ...
            1/2*(mp.F(:,2) + mp.F(:,3)), ...
            1/2*(mp.F(:,2) + mp.F(:,3)), ...
            mp.F(:,4) - 1];
    elseif dim == 3
        strain = [mp.F(:,1)-1, ...
            1/2*(mp.F(:,2)+mp.F(:,4)), ...
            1/2*(mp.F(:,3)+mp.F(:,7)), ...
            1/2*(mp.F(:,2)+mp.F(:,4)), ...
            mp.F(:,5)-1, ...
            1/2*(mp.F(:,6)+mp.F(:,8)), ...
            1/2*(mp.F(:,3)+mp.F(:,7)), ...
            1/2*(mp.F(:,6)+mp.F(:,8)), ...
            mp.F(:,9)-1];
    end
    
    strainrate = (strain - mp.strain)/dt; % Strain rate update

end

% Stress Functions
function PK1 = StressPlaneStrain(mp,Solver)
% Computes the plane strain formulation of the stress for a 2D scenario
    E = Solver.Material.E;
    nu = Solver.Material.poisson;
    C = E/((1+nu)*(1-2*nu));
    
    S = zeros(size(mp.strain,1),4); % Cauchy Stress
    PK1 = zeros(size(mp.strain,1),4);  
    
    S(:,1) = C*( (1-nu)*mp.strain(:,1) + nu*mp.strain(:,4) ); % sigma_11
    S(:,2) = C*( (1-2*nu) * mp.strain(:,2) ); % sigma_12
    S(:,3) = S(:,2); % sigma_21
    S(:,4) = C*( nu*mp.strain(:,1) + (1-nu)*mp.strain(:,4) ); % sigma_22
    
    PK1 = S;
%     % Convert to 1st Poila-Kirchoff Stress
%     PK1(:,1) = S(:,1).*mp.F(:,4) - S(:,2).*mp.F(:,3);
%     PK1(:,2) = -S(:,1).*mp.F(:,2) + S(:,2).*mp.F(:,1);
%     PK1(:,3) = S(:,3).*mp.F(:,4) - S(:,4).*mp.F(:,3);
%     PK1(:,4) = -S(:,3).*mp.F(:,2) + S(:,4).*mp.F(:,1);

end

function PK1 = StressPlaneStress(mp,Solver)
% Computes the plane stress formulation of the stress for a 2D scenario
    E = Solver.Material.E;
    nu = Solver.Material.poisson;
    C = E/(1-nu^2);
    
    S = zeros(size(mp.strain,1),4); % Cauchy Stress
    PK1 = zeros(size(mp.strain,1),4);
    
    S(:,1) = C*( (1)*mp.strain(:,1) + nu*mp.strain(:,4) ); % sigma_11
    S(:,2) = C*( (1-nu) * mp.strain(:,2) ); % sigma_12
    S(:,3) = S(:,2); % sigma_21
    S(:,4) = C*( nu*mp.strain(:,1) + (1)*mp.strain(:,4) ); % sigma_22
    
    PK1 = S;
%     % Convert to 1st Poila-Kirchoff Stress
%     PK1(:,1) = S(:,1).*mp.F(:,4) - S(:,2).*mp.F(:,3);
%     PK1(:,2) = -S(:,1).*mp.F(:,2) + S(:,2).*mp.F(:,1);
%     PK1(:,3) = S(:,3).*mp.F(:,4) - S(:,4).*mp.F(:,3);
%     PK1(:,4) = -S(:,3).*mp.F(:,2) + S(:,4).*mp.F(:,1);
end

function PK1 = Stress2D(mp,Solver)
% Computes the stress for a 2D scenario as if 3rd dimension didn't exist.
% I'm not even sure this has a physical meaning, but it only took a second
    E = Solver.Material.E;
    nu = Solver.Material.poisson;
    
    S = zeros(size(mp.strain,1),4); % Cauchy Stress
    PK1 = zeros(size(mp.strain,1),4);
    
    S(:,1) = E*( (1-nu)*mp.strain(:,1) + nu*mp.strain(:,4) ); % sigma_11
    S(:,2) = E*( (1-2*nu) * mp.strain(:,2) ); % sigma_12
    S(:,3) = S(:,2); % sigma_21
    S(:,4) = E*( nu*mp.strain(:,1) + (1-nu)*mp.strain(:,4) ); % sigma_22
    
    PK1 = S;
%     % Convert to 1st Poila-Kirchoff Stress
%     PK1(:,1) = S(:,1).*mp.F(:,4) - S(:,2).*mp.F(:,3);
%     PK1(:,2) = -S(:,1).*mp.F(:,2) + S(:,2).*mp.F(:,1);
%     PK1(:,3) = S(:,3).*mp.F(:,4) - S(:,4).*mp.F(:,3);
%     PK1(:,4) = -S(:,3).*mp.F(:,2) + S(:,4).*mp.F(:,1);
    
end

function PK1 = Stress3D(mp,Solver)
% Computes the stress for a full 3D simulation
    E = Solver.Material.E;
    nu = Solver.Material.poisson;
    
    S = zeros(size(mp.strain,1),9); % Cauchy Stress
    PK1 = zeros(size(mp.strain,1),9);
    
    S(:,1) = E*( (1-nu)*mp.strain(:,1) + nu*mp.strain(:,5) + nu*mp.strain(:,9)); % sigma_11
    S(:,2) = E*( (1-2*nu) * mp.strain(:,2) ); % sigma_12
    S(:,3) = E*( (1-2*nu) * mp.strain(:,3) ); % sigma_13
    S(:,4) = S(:,2); % sigma_21
    S(:,5) = E*( (nu)*mp.strain(:,1) + (1-nu)*mp.strain(:,5) + nu*mp.strain(:,9)); % sigma_22
    S(:,6) = E*( (1-2*nu) * mp.strain(:,6) ); % sigma_23
    S(:,7) = S(:,3); % sigma_31
    S(:,8) = S(:,6); % sigma_32
    S(:,9) = E*( (nu)*mp.strain(:,1) + nu*mp.strain(:,5) + (1-nu)*mp.strain(:,9)); % sigma_33
    
    PK1 = S;
%     % Convert to 1st Poila-Kirchoff Stress
%     invF = inverse_F(mp.J,mp.F);
%     PK1(:,1) = mp.J(:).*( S(:,1).*invF(:,1) + S(:,2).*invF(:,4) + S(:,3).*invF(:,7) );
%     PK1(:,2) = mp.J(:).*( S(:,1).*invF(:,2) + S(:,2).*invF(:,5) + S(:,3).*invF(:,8) );
%     PK1(:,3) = mp.J(:).*( S(:,1).*invF(:,3) + S(:,2).*invF(:,6) + S(:,3).*invF(:,9) );
%     PK1(:,4) = mp.J(:).*( S(:,4).*invF(:,1) + S(:,5).*invF(:,4) + S(:,6).*invF(:,7) );
%     PK1(:,5) = mp.J(:).*( S(:,4).*invF(:,2) + S(:,5).*invF(:,5) + S(:,6).*invF(:,8) );
%     PK1(:,6) = mp.J(:).*( S(:,4).*invF(:,3) + S(:,5).*invF(:,6) + S(:,6).*invF(:,9) );
%     PK1(:,7) = mp.J(:).*( S(:,7).*invF(:,1) + S(:,8).*invF(:,4) + S(:,9).*invF(:,7) );
%     PK1(:,8) = mp.J(:).*( S(:,7).*invF(:,2) + S(:,8).*invF(:,5) + S(:,9).*invF(:,8) );
%     PK1(:,9) = mp.J(:).*( S(:,7).*invF(:,3) + S(:,8).*invF(:,6) + S(:,9).*invF(:,9) );
    
end

function PK1 = StressWeibenfels(mp,Solver)
% Solves for linear elastic stress based strictly on the deformation
% gradient via Weibenfels equation 25
E = Solver.Material.E;
nu = Solver.Material.poisson;
dim = Solver.domain.dim;
lmbda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

S = zeros(size(mp.x0,1),dim*dim);
for ii = 1:size(mp.x0,1)
    F = reshape(mp.F(ii,:)',dim,dim);
    s = 1./mp.J(ii) * ( lmbda/2*(mp.J(ii).^2-1) + mu*(F*F'-eye(dim)) );
    S(ii,:) = s(:)';
end

    PK1 = S;

end

% Inverses 3x3 deformation gradient explicitly
function invF = inverse_F(detF,F)
invF = zeros(size(F));

invF(:,1) = (F(:,5).*F(:,9) - F(:,6).*F(:,8))./detF;
invF(:,2) = (F(:,3).*F(:,8) - F(:,2).*F(:,9))./detF;
invF(:,3) = (F(:,2).*F(:,6) - F(:,3).*F(:,5))./detF;
invF(:,4) = (F(:,6).*F(:,7) - F(:,4).*F(:,9))./detF;
invF(:,5) = (F(:,1).*F(:,9) - F(:,3).*F(:,7))./detF;
invF(:,6) = (F(:,3).*F(:,4) - F(:,1).*F(:,6))./detF;
invF(:,7) = (F(:,4).*F(:,8) - F(:,5).*F(:,7))./detF;
invF(:,8) = (F(:,2).*F(:,7) - F(:,1).*F(:,8))./detF;
invF(:,9) = (F(:,1).*F(:,5) - F(:,2).*F(:,4))./detF;

end

% Tensor multiplication function
function b = TensMult(A,x,dim)
% Explicitly multiplies matrix equation Ax=b when A has been flattened in
% the typical matlab manner.
A = reshape(A,dim,dim);
x = x';

b = (A*x)';

end







