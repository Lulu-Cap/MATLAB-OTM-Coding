%% Test Problem for OTM Newtonian Fluid Formulation
% A circular domain of radius R is filled with a Newtonian fluid.
% Fluid is initially at rest for t<0
% at t >= 0 the outer boundary moves at tangential velocity u = omega*R
% An analytical solution has been found (SLP)

% Fluid properties
rho = 998; % Density
mu = .001; % Dynamic viscosity
nu = mu/rho; % Kinematic viscosity
omega = 2*pi; % rotational speed [rad/s]

% Find zeros of 1st bessel function of first kind
n = 1000; % Number of zeros to compute 
z = besselzero(1,n,1); % zeros of J_1(z)

% Domain
R = 0.005;
r = linspace(0,R,1000); % radius vector
T = 4;
t = [linspace(0,floor(T/4)),linspace(floor(T/4),T),1000]'; % time vector
u = zeros(length(t),length(r)); % Velocity [t,r]

% Transient solution
for ii = 1:n
    du = 2*omega*R/(z(ii)*besselj(0,z(ii))) * exp(-nu*(z(ii)/R)^2*t) * ...
        besselj(1,z(ii)*r/R);
    u = u + du;
end

% Add steady state solution - Actually, for initialization reasons no SS
% will be used
%u = u + omega*repmat(r,length(t),1);
u = -u;

plot(r,u(1,:)); hold on;
plot(r,u(5,:));
for ii = 10:20:201
    plot(r,u(ii,:)); hold on;
end
plot(r,u(201,:));
hold off;
title("Tangential Velocity Development over Time");
xlabel("Radial Position, r");
ylabel("Tangential Velocity [m/s]");
%legend("t=1","t=2","t=3","t=4","t=5","t=6","t=7","t=8","t=9","t=10");










