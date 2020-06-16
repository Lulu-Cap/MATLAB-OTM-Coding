% Lucas Caparini 53547155 March 13 2020

function Shape = LME_Reg_Neigh(Solver,x_a,x_p)
% LME_Reg: A regularized implementation of the Newton-Raphson iterations used
% in computing LME shape functions and their derivatives. Includes a
% neighbour search algorithm for faster computation.
%
% Introduced by Polyak (2009) and used for LME by Foca (2015) and Kumar et 
% al. (2018). This approach avoids the singularity of matrices in case of
% high locality, which is observed in standard iterations. 
%
%   Solver: Structure containing solver information
%   x_a: Nodal positions
%   x_p: Material point positions
%   h: a measure of the particle spacing (local or global? How does this influence spatial derivatives of shape functions?)
%   gamma: Nondimensional kernel width
%   dim: dimension of the problem
%   Shape: An array of shape function structures which contain relevant
%       information

ITER_MAX = 100; % Maximum Newton-Raphson's to do before giving up
N_mp = length(x_p);
dim = Solver.domain.dim;
R_cutoff = Solver.LMEParam.R_cutoff;
beta = Solver.LMEParam.beta;

% Compute Beta parameter (shape function width around the material point)
%beta = gamma/h^2; % Already contained in Solver, but may need to become
%more complex as algorithm evolves (e.g. beta->beta(x) or anisotropic
%tensor forms)

% Cutoff and neighbour search
NeighNodes = rangesearch(x_a,x_p,R_cutoff); 
for mp = N_mp:-1:1 % Backwards for preallocation
    Shape(mp).numneigh = length(NeighNodes{mp}); % Number of neighbours of mp
    Shape(mp).neigh = NeighNodes{mp}'; % Neighbour ids
    % Preallocate more
    Shape(mp).p = zeros(Shape(mp).numneigh,1); 
    Shape(mp).gradp = zeros(Shape(mp).numneigh,dim);
end

Tol_lambda = 10^-12;
for mp = 1:N_mp
    %mp
    lambda0 = zeros(dim,1);
    lambda = ones(dim,1); % Initial guess for lambda
    iter = 0;
    while norm(lambda-lambda0) > Tol_lambda
        iter = iter + 1;
        lambda0 = lambda; 
        p = P(x_p(mp,:),x_a,Shape(mp),lambda,beta); % I can eliminate the 'p' variable and incorporate into Shape later
        [H,r] = LME_Reg_Hessian(x_p(mp,:),x_a,Shape(mp),p,dim); % Hessian matrix and gradient wrt lambda
        
        H = inv_Hessian(H,dim); % Invert the Hessian Matrix
        
        lambda = lambda0 - H*r; % Newton-Raphson iterations with explicit matrix inversion
        
        if iter > ITER_MAX
            disp("Error, LME not converging");
            iter,mp
            return;
        end
    end
    
    % Send info to output
    Shape(mp).p = p; 
    Shape(mp).gradp = LME_Reg_Grad(x_p(mp,:),x_a,Shape(mp),p,H,dim); % May need adjustment near discontinuities
end
    
end

%%%%%%% Minor Functions %%%%%%%

function grad_p = LME_Reg_Grad(x,xa,Shape,p,H,dim)
% Computes the spatial gradient of shape functions around a material point
%   x: mp location to evaluate at [1 x dim]
%   xa: locations of all nodes w/in neighbourhood [N_h x dim]
%   Shape: Structure containing shape function and neighbour info for this
%           mp
%   p: Shape function evaluations for each neighbour node [N_h x 1]
%   H: Hessian matrix [dim x dim]
%   dim: dimension of problem

    %grad_beta = 0; K_a = 0; % Ignoring the gradient of beta for now. Keep it constant for simplicity
    
    grad_p = zeros(Shape.numneigh,dim);
    for nd = 1:Shape.numneigh
        grad_p(nd,:) = -p(nd)*H*(x-xa(Shape.neigh(nd),:))'; %+ p(nd)*K_a*grad_beta; % This term has been temporarily ignored
    end
    
end

function [H,r] = LME_Reg_Hessian(x,xa,Shape,p,dim)
% Computes the gradient of the shape functions w.r.t. lambda, and also the
% Hessian matrix. The regulization term is included here to improve
% convergence.
%   x: mp location to evaluate at [1 x dim]
%   xa: locations of all nodes w/in neighbourhood [N_h x dim]
%   Shape: Structure containing shape function and neighbour info for this
%           mp
%   p: Shape function evaluations for each neighbour node [N_h x 1]
%   dim: dimension of problem
%   H: Hessian matrix [dim x dim]
%   r: gradient [dim x 1]

    r = zeros(Shape.numneigh,dim); % gradient vector
    H = zeros(dim); % Hessian Matrix
    
    for nd = 1:Shape.numneigh
        r(nd,:) = (x-xa(Shape.neigh(nd),:));
        H = H + r(nd,:)'*r(nd,:)*p(nd);
        r(nd,:) = p(nd)*r(nd,:); 
    end
    
    r = sum(r,1)';
    H = H - r*r';
    % Regulization term: not giving good convergence. Not using for now.
    % Possibly implementing incorrectly
    %H = H + norm(r)*eye(dim); 
end

function invH = inv_Hessian(H,dim)
% Invert the Hessian matrix explicitly
    if dim == 1
        invH = 1/H; % inverse of 1D matrix
    elseif dim ==2
        invH = [H(2,2), -H(1,2); ...
            -H(2,1), H(1,1)] ...
            /(H(1,1)*H(2,2)-H(1,2)*H(2,1)); % Explicit inversion for 2D
    elseif dim == 3
        invH = inv3x3(H); % Explicit inversion for 3D
    else
        disp("ERROR: The LME computation is currently implemented only for 1, 2 and 3 dimensional systems!");
        return;
    end
end

function p = P(x,xa,Shape,lambda,beta)
% Evaluates the shape function given input parameters which may or may not 
% be optimal. 
%   x: point at which to evaluate (generally location of a mp) [dim x 1]
%   xa: nodal locations [#nodes x dim]
%   Shape: Structure containing shape function and neighbour info for this
%           mp
%   lambda: The lambda (Lagrange multiplier) being used. Generally not the
%   optimal one
%   Beta: Locality parameter
%   p: Shape function evaluations at each node 'a'

    % Evaluate partition function, Z(lambda,x)
    p = zeros(Shape.numneigh,1);
    for ii = 1:Shape.numneigh
        p(ii) = z(x,xa(Shape.neigh(ii),:),lambda,beta);
    end
    
    % Evaluate normalize shape function for node 'a'
    p = p/sum(p);
    
end

function f = z(x,xa,lambda,beta)
% Evaluates the exponential component of the shape function given input
% parameters which may or may not be optimal. 
%   x: point at which to evaluate (generally location of a mp) [dim x 1]
%   xa: nodal location [dim x 1]
%   lambda: The lambda (Lagrange multiplier) being used. Generally not the
%   optimal one [dim x 1]
%   Beta: Locality paramter [scalar]
    f = exp(-beta*dot(x-xa,x-xa) + dot(lambda,x-xa));
end

function B = inv3x3(a)
% Inverts the 3 x 3 matrix A explicitly
    determinant = a(1,1)*(a(3,3)*a(2,2)-a(3,2)*a(2,3)) ...
                 -a(2,1)*(a(3,3)*a(1,2)-a(3,2)*a(1,3)) ...
                 +a(3,1)*(a(2,3)*a(1,2)-a(2,2)*a(1,3)); % determinant of the matrix
             
    B = [+(a(3,3)*a(2,2)-a(3,2)*a(2,3)), -(a(3,3)*a(1,2)-a(3,2)*a(1,3)), +(a(2,3)*a(1,2)-a(2,2)*a(1,3)); ...
         -(a(3,3)*a(2,1)-a(3,1)*a(2,3)), +(a(3,3)*a(1,1)-a(3,1)*a(1,3)), -(a(2,3)*a(1,1)-a(2,1)*a(1,3)); ...
         +(a(3,2)*a(2,1)-a(3,1)*a(2,2)), -(a(3,2)*a(1,1)-a(3,1)*a(1,2)), +(a(2,2)*a(1,1)-a(2,1)*a(1,2))] ...
         / determinant; % Inverse
end















