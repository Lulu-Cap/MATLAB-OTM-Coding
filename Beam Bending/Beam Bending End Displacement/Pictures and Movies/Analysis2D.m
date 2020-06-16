%% Compare the 2D beam displacement data to the analytical solution
% June 1st 2020
% Lucas Caparini

%% Displacement
% The displacement is analytically given as y(x)=delta*x^2/L^3*(3L - 2x)
L = Solver.domain.Width;
Ny = Solver.domain.Ny;
Nx = Solver.domain.Nx;
y_start = nd.x_start(1:Ny,2);


%plotting
x = linspace(0,L);
disp = 0.01;
y = disp*x.^2/L^3 .* (3*L - 2*x); % Analytical displacement
y = repmat(y_start,1,length(y))' + repmat(y,length(y_start),1)';
x = repmat(x,Ny,1)' - L/2;

plot(x,y,'b'); hold on;
plot(nd.x1(:,1),nd.x1(:,2),'ro'); hold off
axis equal

% Error analysis
x = unique(nd.x_start(:,1)) + L/2;
y = disp*x.^2/L^3 .* (3*L - 2*x); % Analytical displacement

b = repmat(y',Ny,1);
b = b(:);
err = b - (nd.x1(:,2) - nd.x_start(:,2));
e = (sum(err.^2)/(Nx*Ny))^(1/2); % Lines up pretty well even with this coarse a discretization - e is O(10^-5)
dx = (sum(sum((nd.x1(:,:) - nd.x0(:,:)).^2))/(Nx*Ny))^(1/2); % O(10^-12) which suggests I didn't need to run it for this long. Should I do a convergence analysis?

% We don't observe linear convergence in the coarse to medium
% discretizations (e=2.4e-5 and e=1.2e-5 respectively). I believe this is because the medium discretization had not equilibrated as closely by the time it was finished. 
% It still had dx=O(10^-10) vs O(10^-12). This is hardly a fair comparison.
% I should think of better ways to compare quasi-static test examples.
%% Force & Moment Analysis
I = 1/12*1*Solver.domain.Height.^3; 
L = Solver.domain.Width;
E = Solver.Material.E;
delta = 0.01;

V_theory = 12*delta/(E*I*L^3);

V_act1 = sum(nd.f(end-2:end,2));
V_act2 = sum(nd.f(1:3,2));

% This didn't work so well... I think some of the assumptions going towards
% the I calculations are fucked up.







