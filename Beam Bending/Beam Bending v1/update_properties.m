function [mp] = update_properties(Solver,nd,mp,Shape)
% Updates assorted properties of the material points. It was convenient to
% do these at the same time to avoid storing additional variables
%   Solver: Struct containing solver details
%   nd: struct containing nodal details
%   mp: struct containing material point details
%   Shape: struct array containing shape function evaluation and details
%   for the neighbourhood
dim = Solver.domain.dim;
dt = Solver.time.dt;

% Calculate gradient of the deformation map \grad*\phi_{k \rightarrow k+1}
% (Elsewhere referred to as incremental deformation gradient, F_{p,k->k+1})
for ii = 1:size(Shape,2) % for each mp
    Fincr = zeros(dim); % stores the incremental deformation gradient mapping for this mp
    for jj = 1:Shape(ii).numneigh % Sum effects of nodes in neighbourhood
        node = Shape(ii).neigh(jj);
        Fincr = Fincr + nd.x1(node,:)'*Shape(ii).gradp(jj,:);
    end
    
    % Calculate determinant and volume update
    mp.vol(ii) = mp.vol(ii)*det(Fincr);
    Fold = reshape(mp.F(ii,:),dim,dim);
    Fnext = Fincr*Fold; % **Is this the problem? Try different formulation** 
    mp.Fdot(ii,:) = (Fnext(:)' - Fold(:)') ./ dt; % Fdot = (Fnext-F)/dt
    mp.F(ii,:) = Fnext(:)';
    mp.J(ii) = det(Fnext);
    
end

end



%%%%%%%%%%%%%%% Minor Functions %%%%%%%%%%%%%%%%



