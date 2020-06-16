function nd = Dirichlet_Update(Solver,nd,iter)
%Updates the positions of the Dirichlet end nodes (end of bar). Movement is
%quite slow to avoid annoying stuff

% time
t = Solver.time.t(iter);

% Select the displaced nodes
dim = Solver.domain.dim;
Ny = Solver.domain.Ny;
if dim == 2
    start = Ny-1;
elseif dim == 3
    Nz = Solver.domain.Nz;
    start = Ny*Nz-1;
end
DispNodes = nd.Dirichlet.Nodes(end-start:end);

% Reaches displacment of 1cm in 10 seconds
if t>=0 && t<=10
    nd.Dirichlet.Values(end-start:end,end) = nd.x_start(DispNodes,end) + 0.005*(1-cos(pi/10*t));
elseif t>10
    nd.Dirichlet.Values(end-start:end,end) = nd.x_start(DispNodes,end) + 0.01;
else
    disp("Time is outside of allowable bounds (likely negative). Terminating computation");
    return;
end
end

