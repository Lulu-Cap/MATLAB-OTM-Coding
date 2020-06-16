function plotgrid(X,Y)

% plotgrid: To plot structured grid.
%
% plotgrid(X,Y)
%
% INPUT:
% X (matrix)    - matrix with x-coordinates of gridpoints
% Y (matrix)    - matrix with y-coordinates of gridpoints
% Coded by :    Siva Srinivas Kolukula, PhD      
%               Indian Tsunami Early Warning Centre (ITEWC)
%               Advisory Services and Satellite Oceanography Group (ASG)
%               Indian National Centre for Ocean Information Services (INCOIS)
%               Hyderabad, INDIA
% E-mail   :    allwayzitzme@gmail.com                                        
% web-link :    https://sites.google.com/site/kolukulasivasrinivas/   


if any(size(X)~=size(Y))
   error('Dimensions of X and Y must be equal');
end

[m,n]=size(X);

% Plot grid
% figure
set(gcf,'color','w') ;
hold on
axis equal
axis off
box on
% Plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'b','linewidth',1); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'b','linewidth',1); 
end
hold off
