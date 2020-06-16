function x_new = move_mp(Solver,x_a,Shape)
% Computes the next position of the material points based on the updated
% nodal coordinates
%   Solver: Structure containing solver information
%   x_a: Nodal positions
%   Shape: An array of shape function structures which contain relevant
%       information

dim = Solver.domain.dim;

x_new = zeros(size(Shape,2),dim);

for ii = 1:size(Shape,2) % every mp
    for jj = 1:Shape(ii).numneigh % each nd in neighbourhood
        node = Shape(ii).neigh(jj); % node number
        x_new(ii,:) = x_new(ii,:) + x_a(node,:)*Shape(ii).p(jj);
    end
end

end

