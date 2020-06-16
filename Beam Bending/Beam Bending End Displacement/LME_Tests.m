% Lucas Caparini 53547155 May 24 2020
% Testing convergence of LME functions with Newton-Raphson iterations
%
% Results: Difficulty converging for larger values of gamma (>5.2
% doesn't work). Also difficult to converge on irregular grids
%
% The iterations can't seem to converge if the value of beta makes
% the function too local compared to the nodal spacing. This happens on
% irregular grids by accident commonly.
% The Lagrange multipliers are very large or small in these cases, which
% apparently makes convergence difficult. 
% 
% It may be worth using an additional, robust, non-gradient optimization 
% method initially to aid convergence robustness. Nelder-Mead simplex 
% algorithm has been used by some authors.
 
clear,clc

h = 0.2;
x = [-1:h:1];
[x,y] = meshgrid(x);
for ii = 1:length(x(:))
    x_a(ii,:) = [x(ii),y(ii)] ;%+ (rand(1,2)-0.5); %nodal positions
end
x_a(61,:) = [0 0];

x = [-0.5:0.01:0.5];
[x,y] = meshgrid(x);
for ii = 1:length(x(:))
    x_p(ii,:) = [x(ii),y(ii)]; %nodal positions
end

plot(x_a(:,1),x_a(:,2),'ro'); hold on;
plot(x_p(:,1),x_p(:,2),'bx');

ITER_MAX = 100; % Maximum Newton-Raphson's to do before giving up
N_mp = size(x_p,1);
dim = 2;
beta = 4.0/h^2;
R_cutoff = sqrt(-log(10^-12)./beta);

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
reg_flag = 0; % Flag for initializing Foca's (2015) regularization if necessary
maxIter = 0;
for mp = 1:N_mp
    X = x_p(mp,:);
    X_a = x_a(Shape(mp).neigh,:);
    
    lambda0 = zeros(dim,1);
    lambda = ones(dim,1); % Initial guess for lambda
    iter = 0;
    
    while norm(lambda-lambda0) > Tol_lambda
        iter = iter + 1;
        lambda0 = lambda;
        
        p = P(X,X_a,lambda,beta,h); % Guess at shape functions
        [H,r] = LME_Reg_Hessian(X,X_a,p,h,dim,reg_flag); % Hessian matrix and gradient wrt lambda
        
        Hinv = inv_Hessian(H,dim); % Invert the Hessian Matrix
        
        lambda = lambda0 - Hinv*r; % Newton-Raphson iterations with explicit matrix inversion
        
        if iter > ITER_MAX || isnan(lambda(1)) || isnan(lambda(2))
            if reg_flag == 1
                disp("Error, LME not converging");
                iter,mp
                pause;
                break;
            elseif reg_flag == 0
                reg_flag = 1;
                lambda0 = zeros(dim,1);
                lambda = ones(dim,1);
                iter = 0;
            end    
        end
    end
    
    if maxIter<iter
        maxIter = iter;
    end
    
    % Send info to output
    Shape(mp).p = p; 
    Shape(mp).gradp = LME_Reg_Grad(X,X_a,p,Hinv,h,dim); % May need adjustment near discontinuities
end
    

test = 61; %Node to test (right in the center)
numneigh = 0;
neigh = zeros(size(x_p,1),1);
index = zeros(size(x_p,1),1);
for ii = 1:size(x_p,1)
    for jj = 1:Shape(ii).numneigh
        if Shape(ii).neigh(jj) == test
            numneigh = numneigh + 1;
            neigh(numneigh) = ii;
            index(numneigh) = jj;
            p(numneigh) = Shape(ii).p(jj);
            gradp(numneigh,:) = Shape(ii).gradp(jj,:);
        end

    end
end
figure;
plot3(x(neigh(1:numneigh)),y(neigh(1:numneigh)),p,'bo','MarkerSize',1);
xlabel('x'),ylabel('y');zlabel('z');title("P");
figure;
plot3(x(neigh(1:numneigh)),y(neigh(1:numneigh)),gradp(:,1),'bo','MarkerSize',1);
xlabel('x'),ylabel('y');zlabel('z');title("dP/dx");
figure;
plot3(x(neigh(1:numneigh)),y(neigh(1:numneigh)),gradp(:,2),'bo','MarkerSize',1);
xlabel('x'),ylabel('y');zlabel('z');title("dP/dy");







%%%%%%% Minor Functions %%%%%%%

function grad_p = LME_Reg_Grad(x,xa,p,H,h,dim)
% Computes the spatial gradient of shape functions around a material point
%   x: mp location to evaluate at [1 x dim]
%   xa: locations of all nodes w/in neighbourhood [N_h x dim]
%   Shape: Structure containing shape function and neighbour info for this
%           mp
%   p: Shape function evaluations for each neighbour node [N_h x 1]
%   H: Inverse of Hessian matrix [dim x dim]
%   dim: dimension of problem

    %grad_beta = 0; K_a = 0; % Ignoring the gradient of beta for now. Keep it constant for simplicity
    NumNeigh = size(xa,1);
    
    grad_p = zeros(NumNeigh,dim);
    for nd = 1:NumNeigh
        grad_p(nd,:) = -p(nd)*H*(x-xa(nd,:))'/h^2; %+ p(nd)*K_a*grad_beta; % This term has been temporarily ignored
    end
    
end

function [H,r] = LME_Reg_Hessian(x,xa,p,h,dim,flag)
% Computes the gradient of the shape functions w.r.t. lambda, and also the
% Hessian matrix. The regulization term is included here to improve
% convergence.
%   x: mp location to evaluate at [1 x dim]
%   xa: locations of all nodes w/in neighbourhood [N_h x dim]
%   p: Shape function evaluations for each neighbour node [N_h x 1]
%   h: Nodal spacing parameter
%   dim: dimension of problem
%   H: Hessian matrix [dim x dim]
%   r: gradient [dim x 1]
    
    NumNeigh = size(xa,1);

    r = zeros(NumNeigh,dim); % gradient vector
    H = zeros(dim); % Hessian Matrix
    
    for nd = 1:NumNeigh
        r(nd,:) = (x-xa(nd,:))/h;
        H = H + r(nd,:)'*r(nd,:)*p(nd);
        r(nd,:) = p(nd)*r(nd,:); 
    end
    
    r = sum(r,1)';
    H = H - r*r';
    % Regulization term: Aids convergence. Proposed by Foca (2015), but a
    % bit slow
    if flag == 1
        H = H + norm(r)*eye(dim); 
    end
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

function p = P(x,xa,lambda,beta,h)
% Evaluates the shape function given input parameters which may or may not 
% be optimal. 
%   x: point at which to evaluate (generally location of a mp) [dim x 1]
%   xa: nodal locations w/in neighbourhood of x [#nodes x dim]
%   lambda: The lambda (Lagrange multiplier) being used. Generally not the
%   optimal one
%   Beta: Locality parameter
%   h: Nodal spacing parameter
%   p: Shape function evaluations at each node 'a'
    
    NumNeigh = size(xa,1);
    
    lambda = repmat(lambda',NumNeigh,1);
    x = repmat(x,NumNeigh,1);
    f = -beta*dot((x-xa)',(x-xa)')' + dot(lambda'/h,(x-xa)')';
    p = exp(f);
    p = p./sum(p);

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















