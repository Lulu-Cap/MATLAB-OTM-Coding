function Dirichlet_new = Dirichlet_Update(Dirichlet_old,t)
% Updates the Dirichlet boundaries with time. In this case the indentor
% gradually moves down into a final position
%   Dirichlet_old: Previous Dirichlet nodal values (not really used)
%   t: time
%   Dirichlet_new: New Dirichlet nodal values (returned)
Dirichlet_new = Dirichlet_old;
dispMax = 0.1;

if t<0
    % Negative time not allowed
    disp("ERROR: Negative values for time not allowed");
    return
elseif t<2
    Dirichlet_new(end-2:end,2) = 1 - sin(pi/2*t)*dispMax;
elseif t>=2
    Dirichlet_new(end-2:end,2) = 1  - dispMax;
end

end

