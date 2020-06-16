function nd = Dirichlet_Update(Solver,nd,iter)
%Updates the positions of the Dirichlet end nodes (end of bar). Movement is
%quite slow to avoid an exploding solution

% time
t = Solver.time.t(iter); 
%t = 0; % Solution explodes even when nodes are stationary!
Omega = Solver.domain.Omega; % angular velocity
Theta = Omega*t; % angle change from original position
DirichletNd = nd.Dirichlet.Nodes; % displaced nodes

% Outer nodes move with a constant angular velocity (rotation matrix about origin)
x = nd.x_start(DirichletNd,:);
nd.Dirichlet.Values(:,1) = x(:,1)*cos(Theta) - x(:,2)*sin(Theta);
nd.Dirichlet.Values(:,2) = x(:,1)*sin(Theta) + x(:,2)*cos(Theta);
end

