function force_s = Hourglass_Stabilization(Solver,nd,mp,Shape,RoI)
%Implements the hourglass stabilization algorithm proposed by Weibenfels
%and Wriggers.
% Anisotropic stabilization is an interesting future option.
%   Solver: Struct containing solver details
%   nd: struct containing nodal details
%   mp: struct containing material point details
%   Shape: struct array containing shape function evaluation and details
%   for the neighbourhood
%   RoI: struct containing info on Range of Influence of nodes
%   force_s: Stabilized forces

dim = Solver.domain.dim;
eps = Solver.LMEParam.stab; % Stabilizing parameter

df = zeros(size(nd.x1,1),dim);
for node = 1:size(nd.x0,1) % for every node
    for jj = 1:RoI(node).numneigh % for every mp in the neighbourhood
        matpnt = RoI(node).neigh(jj); % material point number
        index = RoI(node).index(jj); % position of this node in Shape structure
        Fincr = reshape(mp.Fincr(matpnt,:),dim,dim);
        dx0 = nd.x0(node,:)' - mp.x0(matpnt,:)';
        dx1 = nd.x1(node,:)' - mp.x1(matpnt,:)';
        err = (dx1 - Fincr*dx0)./norm(dx0); % Error term: deviation from linear shape functions
        
        df(node,:) = df(node,:) - eps*Shape(matpnt).p(index)*err';
        
    end
end

force_s = nd.force + df;

end

