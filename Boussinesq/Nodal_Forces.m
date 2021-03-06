function [mp,force] = Nodal_Forces(Solver,nd,mp,Shape,RoI)
% Compute the nodal forces based on the consitutive equation specified
%   Solver: Struct containing solver details
%   nd: struct containing nodal details
%   mp: struct containing material point details
%   Shape: struct array containing shape function evaluation and details
%   for the neighbourhood

% Solver.Mat contains the parameters for the specified constitutive
% relation.
switch Solver.Material.ConstitutiveEq(1)
    case "SolidLinearElastic"
        [mp,force] = LinearElasticForces(Solver,nd,mp,Shape,RoI);
    case "SolidNonlinearElastic"
        input("No nonlinear solid constitutive equations are defined yet. \nPlease exit this session and add the desired relation.");
        return;
    case "LiquidNewtonian"
        input("No Newtonian liquid constitutive equations are defined yet. \nPlease exit this session and add the desired relation.");
        return;
    case "LiquidNonNewtonian"
        input("No Non-Newtonian liquid constitutive equations are defined yet. \nPlease exit this session and add the desired relation.");
        return;
    otherwise 
        input("Not a known constitutive relation. \nPlease exit this session and add the desired relation.");
        return;
end

end

