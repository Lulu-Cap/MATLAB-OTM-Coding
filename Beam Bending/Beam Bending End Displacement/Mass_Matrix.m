function M = Mass_Matrix(nd,mp,Shape,RoI)
% Computes the nodal mass matrix of the system (either lumped or distributed)
%   nd: struct containing nodal details
%   mp: struct containing material point details
%   Shape: struct array containing shape function evaluation and details
%   for the neighbourhood
%   RoI: struct containing range of influence info for each node

if nd.MassMatrixStyle == "Lumped"
    % Compute lumped mass matrix. Ref. Bo's paper
    % Note: This is the row sum of the consistent matrix below
    M = zeros(length(nd.x1),1);
    for ii = 1:length(nd.x1) % For each node ii 
        m = 0;
        for jj = 1:RoI(ii).numneigh
            matpnt = RoI(ii).neigh(jj); % Material point w/in range
            index = RoI(ii).index(jj); % Index of node w/in Shape
            m = m + mp.mass(matpnt)*Shape(matpnt).p(index); % sum contributions from each mp w/in range
        end
        
        M(ii,1) = m;

    end
    
elseif nd.MassMatrixStyle == "Distributed"
    % Compute distributed mass matrix. Slower, but most
    % papers use it
    % Row summing yields the lumped mass matrix
    
    % DO NOT USE: Current implementation is far too slow!!!
    
    M = zeros(length(nd.x1));
    for ii = 1:length(mp.x1) % Cycle through each mp
        for jj = 1:Shape(ii).numneigh
            for kk = jj:Shape(ii).numneigh % take advantage of symmetry
                row = Shape(ii).neigh(jj);
                col = Shape(ii).neigh(kk);
                M(row,col) = M(row,col) + mp.mass(ii)*Shape(ii).p(jj)*Shape(ii).p(kk);
                M(col,row) = M(row,col); % Symmetric
            end
            
        end
        %ii
    end
    M = sparse(M); % It will be sparse in most cases
else
    % Not a recognized mass matrix type
    disp("ERROR: Not a recognized mass matrix type.");
    return;
end

end

%%%%%%% Minor Functions %%%%%%%