Nx = Solver.domain.Nx;
Ny = Solver.domain.Ny;

% stress plots
Fxx = scatteredInterpolant(mp.x_start(:,:),mp.stress(:,1),'linear','linear');
Fxy = scatteredInterpolant(mp.x_start(:,:),mp.stress(:,2),'linear','linear');
Fyy = scatteredInterpolant(mp.x_start(:,:),mp.stress(:,4),'linear','linear');

Nd_stress_xx = Fxx(nd.x_start(:,:));
Nd_stress_xy = Fxy(nd.x_start(:,:));
Nd_stress_yy = Fyy(nd.x_start(:,:));

dx = round(10^6*(nd.x_start(Ny+1,1)-nd.x_start(1,1)))/10^6;
xv = [min(nd.x_start(:,1)):dx:max(nd.x_start(:,1))]';
yv = [min(nd.x_start(:,2)):dx:max(nd.x_start(:,2))]';
[X,Y] = meshgrid(xv,yv);
Zxx = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_stress_xx,X,Y);
Zxy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_stress_xy,X,Y);
Zyy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_stress_yy,X,Y);

figure;
subplot(3,1,1)
contourf(X,Y,Zxx); axis equal;
title("\sigma_{xx}"); 
colorbar('southoutside');
subplot(3,1,2)
contourf(X,Y,Zxy); axis equal;
title("\sigma_{xy}");
colorbar('southoutside');
subplot(3,1,3)
contourf(X,Y,Zyy); axis equal;
title("\sigma_{yy}");
colorbar('southoutside');

% strain plots
Fxx = scatteredInterpolant(mp.x_start(:,:),mp.strain(:,1),'linear','linear');
Fxy = scatteredInterpolant(mp.x_start(:,:),mp.strain(:,2),'linear','linear');
Fyy = scatteredInterpolant(mp.x_start(:,:),mp.strain(:,4),'linear','linear');

Nd_strain_xx = Fxx(nd.x_start(:,:));
Nd_strain_xy = Fxy(nd.x_start(:,:));
Nd_strain_yy = Fyy(nd.x_start(:,:));

dx = round(10^6*(nd.x_start(Ny+1,1)-nd.x_start(1,1)))/10^6;
xv = [min(nd.x_start(:,1)):dx:max(nd.x_start(:,1))]';
yv = [min(nd.x_start(:,2)):dx:max(nd.x_start(:,2))]';
[X,Y] = meshgrid(xv,yv);
Zxx = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_strain_xx,X,Y);
Zxy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_strain_xy,X,Y);
Zyy = griddata(nd.x_start(:,1),nd.x_start(:,2),Nd_strain_yy,X,Y);
figure;
subplot(3,1,1)
contourf(X,Y,Zxx); axis equal;
title("\epsilon_{xx}"); 
colorbar('southoutside');
subplot(3,1,2)
contourf(X,Y,Zxy); axis equal;
title("\epsilon_{xy}");
colorbar('southoutside');
subplot(3,1,3)
contourf(X,Y,Zyy); axis equal;
title("\epsilon_{yy}");
colorbar('southoutside');

clear Nx Ny Fxx Fxy Fyy Nd_stress_xx Nd_stress_xy Nd_stress_yy dx xv yv X Y Zxx Zxy Zyy 
clear Nd_strain_xx Nd_strain_xy Nd_strain_yy


%%
force = zeros(size(nd.x1));
dim = 2;
m = mean(mp.vol);
Free = nd.Free; % free nodes
    for node = 1:size(nd.x1,1)
        for jj = 1:RoI(node).numneigh
            matpnt = RoI(node).neigh(jj);
            index = RoI(node).index(jj); % The index of the node within the mp's neighbour list
            force(node,1) = force(node,1) - Shape(matpnt).gradp(index(1),1)*mp.vol(matpnt)/m;
            force(node,2) = force(node,2) - Shape(matpnt).gradp(index(1),2)*mp.vol(matpnt)/m;
        end
    end

function b = TensMult(A,x,dim)
% Explicitly multiplies matrix equation Ax=b when A has been flattened in
% the typical matlab manner.
A = reshape(A,dim,dim);
if size(x,1) ~= dim
    x = x';
end

b = (A*x)';

end


