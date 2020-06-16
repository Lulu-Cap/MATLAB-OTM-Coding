function x_new = update_nd(Solver,nd)
% Returns the updated nodal position after explicit time integration.
%   Solver: Struct containing solver details
%   nd: struct containing nodal details
%   Shape: struct array containing shape function evaluation and details
%       for the neighbourhood

FreeNd = nd.Free;
DirichletNd = nd.Dirichlet.Nodes;
dt = Solver.time.dt;
dim = Solver.domain.dim;
c = 1-Solver.Material.damping;

x_new = zeros(size(nd.x1,1),dim);

% Impose boundary conditions
x_new(DirichletNd,:) = nd.Dirichlet.Values;

if nd.MassMatrixStyle == "Lumped"
    x_new(FreeNd,:) = nd.x0(FreeNd,:) + dt./repmat(nd.mass(FreeNd),1,dim) .* (c*nd.l(FreeNd,:) + dt*nd.f(FreeNd,:));
elseif nd.MassMatrixStyle == "Distributed"
    % LU Factorization b/c of repeated matrix inversion
    [L,U,P] = lu(nd.mass);
    for d = 1:dim
        tmp = L\(P*(nd.l(:,d) + dt/2 * nd.f(:,d)));
        x_new(:,d) = nd.x0(:,d) + dt * U\tmp;
    end
else
    % Not a recognized mass matrix type
    disp("ERROR: Not a recognized mass matrix type.");
    return;
end

end

