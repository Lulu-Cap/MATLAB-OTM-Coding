function mom = Momentum(Solver,nd,mp,Shape)
%Computes nodal momentum using first order approximation. Exact formulation
%changes depending on mass matrix form
%   Solver: Struct containing solver details
%   nd: struct containing nodal details
%   mp: struct containing material point details
%   Shape: struct array containing shape function evaluation and details
%   for the neighbourhood

dim = Solver.domain.dim;
dt = Solver.time.dt;

mom = zeros(length(nd.x1),dim);

if nd.MassMatrixStyle == "Lumped"
    % Compute momentum using lumped mass matrix style. Refer to Bo's paper
    % Eq 42. (Some of his later stuff doesn't make sense...)
    mom = repmat(nd.mass,1,dim).*(nd.x1-nd.x0)/dt;
    
elseif nd.MassMatrixStyle == "Distributed"
    % Compute momentum for use w/ distributed mass matrix. Slower, but most
    % papers use it
    for ii = 1:length(nd.x1) % For each node ii 
        m = 0;
        for jj = 1:RoI(ii).numneigh
            matpnt = RoI(ii).neigh(jj); % Material point w/in range
            index = RoI(ii).index(jj); % Index of node w/in Shape
            m = m + mp.mass(matpnt)*Shape(matpnt).p(index) * ...
                    (mp.x1(matpnt)-mp.x0(matpnt)) / dt; % sum contributions from each mp w/in range
        end
        
        mom(ii,:) = m;

    end
    
else
    % Not a recognized mass matrix type
    disp("ERROR: Not a recognized mass matrix type.");
    return;
end

end

