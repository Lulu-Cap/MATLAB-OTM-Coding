% Lucas Caparini 53547155 March 13 2020

% This function needs some major improvement. It is quite slow, and bloody
% impossible to understand what is happening. Sort of like this entire shit
% code

function RoI = NodalNeighbours(Shape,numNodes,x_a,x_p,Solver)
% Finds the number and indices of material points within the range of
% influence of each node.
%   Shape: Shape function structure (defined in LME function)
%   RoI: Range of Influence structure. Contains # and indices of nearby mps
%       for each node

% Initialize structure
for ii = numNodes:-1:1
    RoI(ii).numneigh = 0;
    RoI(ii).neigh(3) = 0;
end

for mp = 1:size(Shape,2)
    for ii=1:Shape(mp).numneigh
        node = Shape(mp).neigh(ii); % The node being operated upon
        
        RoI(node).numneigh = RoI(node).numneigh + 1; % Increment the number of mp neighbours of this node
        RoI(node).neigh(RoI(node).numneigh) = mp; % Assign material point as nodal neighbour
        RoI(node).index(RoI(node).numneigh) = ii; % Shape function index for this node w/in Shape(mp).p(index)
    end
end

%Alternative formulation which should be faster - Actually, its only
%slightly faster. Sadness. 
% R_cutoff = Solver.LMEParam.R_cutoff;
% N_nd = size(x_a,1);
% 
% tstart2 = tic;
% NeighMPs = rangesearch(x_p,x_a,R_cutoff);
% for nd = N_nd:-1:1
%     RoI2(nd).numneigh = length(NeighMPs{nd});
%     RoI2(nd).neigh = NeighMPs{nd}';
%     
%     % Index this node corresponds to in the Shape data structure
%     for mp = 1:RoI2(nd).numneigh
%         j = RoI2(nd).neigh(mp);
%         RoI2(nd).index(mp) = find(Shape(j).neigh == nd);
%     end
%    
% end
% toc(tstart2)

end


