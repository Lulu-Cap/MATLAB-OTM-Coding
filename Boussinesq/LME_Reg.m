function [phi, grad_phi] = LME_Reg(x_a,x_p,gamma,h,dim)
% LME_Reg: A regularized implementation of the Newton-Raphson iterations used
% in computing LME shape functions and their derivatives. 
%
% Introduced by Polyak (2009) and used for LME by Foca (2015) and Kumar et 
% al. (2018). This approach avoids the singularity of matrices in case of
% high locality, which is observed in standard iterations. 
%
%   x_a: Nodal positions
%   x_p: Material point positions
%   h: a measure of the particle spacing (local or global? How does this influence spatial derivatives of shape functions?)
%   gamma: Nondimensional kernel width
%   dim: dimension of the problem
%   phi: shape function matrix. Rows correspond to nodes, columns to mps
%   grad_phi: gradients of each shape function 

% Compute Beta parameter (shape function width around the material point)
beta = gamma/h^2;

% Cutoff and neighbour search
% I will improve it to use this soon enough
Solver.Tol_support = 10^-12; % Cutoff tolerance
Solver.R_cutoff = sqrt(-log(Solver.Tol_support)./beta); % Neighbour cutoff radius


% Compute shape functions
ITER_MAX = 100;
N_nd = length(x_a);
N_mp = length(x_p);
phi = zeros(N_nd,N_mp);
grad_phi = zeros(N_nd,N_mp,dim);

Tol_lambda = 10^-13;
for mp = 1:N_mp
    mp
    lambda0 = zeros(dim,1);
    lambda = ones(dim,1); % Initial guess for lambda
    iter = 0;
    while norm(lambda-lambda0) > Tol_lambda
        iter = iter + 1
        lambda0 = lambda; 
        p = P(x_p(mp,:),x_a,lambda,beta);
        [H,r] = LME_Reg_Hessian(x_p(mp,:),x_a,p,dim); % Hessian matrix and gradient wrt lambda
        
        H = inv_Hessian(H,dim); % Invert the Hessian Matrix
        
        lambda = lambda0 - H*r; % Newton-Raphson iterations with explicit matrix inversion
        
        if iter > ITER_MAX
            disp("Error, LME not converging");
            iter,mp
            return;
        end
    end
    phi(:,mp) = p;
    grad_phi(:,mp,:) = LME_Reg_Grad(x_p(mp,:),x_a,p,H,dim); % May need adjustment near discontinuities
end
    
    



end

function grad_p = LME_Reg_Grad(x,xa,p,H,dim)
% Compute the gradients of the shape functions around a material point

    grad_beta = 0; K_a = 0; % Ignoring the gradient of beta for now. Keep it constant for simplicity
    
    grad_p = zeros(length(xa),dim);
    for nd = 1:length(xa)
        grad_p(nd,:) = -p(nd)*H*(x-xa(nd,:))' + p(nd)*K_a*grad_beta;
    end
    
end

function [H,r] = LME_Reg_Hessian(x,xa,p,dim)
% Computes the gradient of the shape functions w.r.t. lambda, and also the
% Hessian matrix. The regulization term is included here to improve
% convergence.
%   x: mp location to evaluate at [1 x dim]
%   xa: locations of all nodes w/in neighbourhood [N_h x dim]
%   p: Shape function evaluations for each neighbour node [N_h x 1]
%   dim: dimension of problem
%   H: Hessian matrix [dim x dim]
%   r: gradient [dim x 1]

    r = zeros(length(xa),dim); % gradient vector
    H = zeros(dim); % Hessian Matrix
    
    for nd = 1:length(xa)
        r(nd,:) = (x-xa(nd,:));
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

function p = P(x,xa,lambda,beta)
% Evaluates the shape function given input parameters which may or may not 
% be optimal. 
%   x: point at which to evaluate (generally location of a mp) [dim x 1]
%   xa: nodal locations [#nodes x dim]
%   lambda: The lambda (Lagrange multiplier) being used. Generally not the
%   optimal one
%   Beta: Locality parameter
%   p: Shape function evaluations at each node 'a'

    % Evaluate partition function, Z(lambda,x)
    p = zeros(length(xa),1);
    for ii = 1:length(xa)
        p(ii) = z(x,xa(ii,:),lambda,beta);
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



















