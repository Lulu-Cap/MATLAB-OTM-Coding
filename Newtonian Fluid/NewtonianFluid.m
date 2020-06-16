function [mp,force] = NewtonianFluid(Solver,nd,mp,Shape,RoI)
% Computes the force on a set of nodes based on a Newtonian Fluid model. 
%   Solver: Struct containing solver details
%   nd: struct containing nodal details
%   mp: struct containing material point details
%   Shape: struct array containing shape function evaluation and details
%   for the neighbourhood
    dim = Solver.domain.dim;
    dt = Solver.time.dt; %#ok<NASGU>
    PressureType = Solver.Material.ConstitutiveEq(2);
    ViscosityType = Solver.Material.ConstitutiveEq(3);
    
    % Compute the stress tensor's pressure component
    if PressureType == "Tait"
        stress_elastic = Pressure_Tait(Solver,mp);
    elseif PressureType == "Adiabatic"
        stress_elastic = Pressure_Adiabatic(Solver,mp);
    elseif PressureType == "Isothermal"
        stress_elastic = Pressure_Isothermal(Solver,mp);
    else
        % Not recognized strain measure
        input("Not a recognized pressure relation. Please exit and select a known measure (Taite, Adiabatic, Isothermal, etc)");
        return;
    end
    
    % Viscous term of fluid stress
    if ViscosityType == "Inviscid"
        stress_viscous = 0; % No viscous terms
    elseif ViscosityType == "Viscous"
        stress_viscous = ViscousStress(Solver,mp); 
    end
    mp.stress = stress_elastic + stress_viscous;
    
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

function stress = Pressure_Tait(Solver,mp)
% Computes the pressure term of the Cauchy stress due to Taite's relation. For
% liquid water n=7

K = Solver.Material.BulkMod; % Bulk Modulus
n = Solver.Material.n; % Coefficient of Taite's equation
dim = Solver.domain.dim; 
p0 = 0; % starting pressure of the problem

I = eye(dim);
stress = repmat(I(:)',size(mp.x1,1),1); % Identity matrices, flattened
pressure = -(p0 + K/n*(mp.J.^(-n) - 1)); % Pressure at each mp (using solid pressure convention).
stress = repmat(pressure,1,dim*dim).*stress; % multiplied by identity to get Cauchy stress

detF = mp.F(:,1).*mp.F(:,4) - mp.F(:,2).*mp.F(:,3);
invF = inverse_F(detF,mp.F); % inverse of deformation gradient
tmp = invF(:,2); invF(:,2) = invF(:,3); % Transpose 
invF(:,3) = tmp;

stress = repmat(mp.J,1,dim^2).*pressure.*invF;

end

function stress = Pressure_Adiabatic(Solver,mp)
% Computes the pressure term of the Cauchy stress due to an Adiabatic relation

C = Solver.Material.PressureConstant; % Initialization constant for pressure
gamma = Solver.Material.IGRatio; % gamma = c_p/c_v
dim = Solver.domain.dim; 

I = eye(dim);
stress = repmat(I(:)',size(mp.x1,1),1); % Identity matrices, flattened
pressure = C*mp.dens.^gamma; % Pressure at each mp
stress = repmat(pressure,1,dim*dim).*stress; % multiplied by identity

end

function stress = Pressure_Isothermal(Solver,mp)
% Computes the pressure term of the Cauchy stress due to an Isothermal relation

K = Solver.Material.BulkMod; % Bulk Modulus
dim = Solver.domain.dim; 
p0 = 0; % starting pressure of the problem

I = eye(dim);
stress = repmat(I(:)',size(mp.x1,1),1); % Identity matrices, flattened
pressure = p0 + K/2*(mp.J - mp.J.^(-1)); % Pressure at each mp
stress = repmat(pressure,1,dim*dim).*stress; % multiplied by identity

end

function stress = ViscousStress(Solver,mp)
% Computes the viscous component of the Cauchy stress 

mu = Solver.Material.visc; % dynamic viscosity
dim = Solver.domain.dim; 

stress = zeros(size(mp.F));
for ii = 1:size(mp.x1,1)
    F = reshape(mp.F(ii,:),dim,dim);
    Fdot = reshape(mp.Fdot(ii,:),dim,dim);
    if cond(F)==Inf % Prevents stationary fluid from returning inf stress
        input("Deformation gradient tensor is singular... Exiting");
        return;
    else
        L = Fdot/F; % velocity gradient tensor
        D = 1/2*(L + L'); % stretch rate tensor
        D_dev = D - 1/3*trace(D)*eye(dim); % deviatoric portion of stretch rate tensor
        stress(ii,:) = 2*mu*D_dev(:);
    end
end

end

% Inverses 3x3 deformation gradient explicitly
function invF = inverse_F(detF,F)
invF = zeros(size(F));

if size(F,2) == 9
    invF(:,1) = (F(:,5).*F(:,9) - F(:,6).*F(:,8))./detF;
    invF(:,2) = (F(:,3).*F(:,8) - F(:,2).*F(:,9))./detF;
    invF(:,3) = (F(:,2).*F(:,6) - F(:,3).*F(:,5))./detF;
    invF(:,4) = (F(:,6).*F(:,7) - F(:,4).*F(:,9))./detF;
    invF(:,5) = (F(:,1).*F(:,9) - F(:,3).*F(:,7))./detF;
    invF(:,6) = (F(:,3).*F(:,4) - F(:,1).*F(:,6))./detF;
    invF(:,7) = (F(:,4).*F(:,8) - F(:,5).*F(:,7))./detF;
    invF(:,8) = (F(:,2).*F(:,7) - F(:,1).*F(:,8))./detF;
    invF(:,9) = (F(:,1).*F(:,5) - F(:,2).*F(:,4))./detF;

elseif size(F,2) == 4
    invF(:,1) = F(:,4)./detF;
    invF(:,2) = -F(:,2)./detF;
    invF(:,3) = -F(:,3)./detF;
    invF(:,4) = F(:,1)./detF;
    
end

end

% Tensor multiplication function
function b = TensMult(A,x,dim)
% Explicitly multiplies matrix equation Ax=b when A has been flattened in
% the typical matlab manner.
A = reshape(A,dim,dim);
if size(x,1) ~= dim
    x = x';
end

b = (A*x)';

end









