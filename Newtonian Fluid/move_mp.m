function x_new = move_mp(Solver,x_a,Shape)
% Computes the next position of the material points based on the updated
% nodal coordinates
%   Solver: Structure containing solver information
%   x_a: Nodal positions (updated)
%   Shape: An array of shape function structures which contain relevant
%       information

dim = Solver.domain.dim;
N_mp = size(Shape,2);

x_new = zeros(N_mp,dim);
min = 1000;
for ii = 1:N_mp % every mp
    
    neigh = Shape(ii).neigh;
    p = Shape(ii).p;
    x = x_a(neigh,:).*repmat(p,1,dim);
    x_new(ii,:) = sum(x);
%     for jj = 1:Shape(ii).numneigh % each nd in neighbourhood
%         node = Shape(ii).neigh(jj); % node number
%         x_new(ii,:) = x_new(ii,:) + x_a(node,:)*Shape(ii).p(jj);
%         
        if Shape(ii).numneigh < min
            min = Shape(ii).numneigh;
        end
%     end
end

end

