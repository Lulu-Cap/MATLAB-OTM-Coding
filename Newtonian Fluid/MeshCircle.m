function [X,Y] = MeshCircle(C,R,M,N) 
% Function to mesh the given circle with quadrilaterals using two boundary
% Transfinite Interpolation (TFI)
% version 1 : 23/June/2016 
% MeshCircle(C,R,M,N)
%
% Input:
%      C  = center of circle
%      R  = Radius of the circle 
%      M  = Number of points along the circumference 
%      N  =  Number of points along the radius 
% Output:
%      plot of circle meshed with quadrilaterals 
%
% Coded by :    Siva Srinivas Kolukula, PhD      
%               Indian Tsunami Early Warning Centre (ITEWC)
%               Advisory Services and Satellite Oceanography Group (ASG)
%               Indian National Centre for Ocean Information Services (INCOIS)
%               Hyderabad, INDIA
% E-mail   :    allwayzitzme@gmail.com                                        
% web-link :    https://sites.google.com/site/kolukulasivasrinivas/   

%% Input check 
if nargin ~= 4
    error('Please check the inputs') ;
elseif length(C) ~= 2
    error('center of circle must have two coordinates') ;
end

%% Square at the centre 
S = R/4. ;
P1 = [C(1)-S C(2)+S] ;
P2 = [C(1)+S C(2)+S] ;
P3 = [C(1)-S C(2)-S] ;
P4 = [C(1)+S C(2)-S] ;
%% Mesh square at the center of circle
% Draw lines
NS = round(M/4) ; 
t = linspace(0,1,NS) ;
L1 = zeros(NS,2) ;
L2 = zeros(NS,2) ;
for i = 1:NS
    L1(i,:) = P1+t(i)*(P2-P1) ;
    L2(i,:) = P3+t(i)*(P4-P3) ;
end
% Generate mesh along square with Two boundary TFI
X1 = zeros(NS,NS) ;
Y1 = zeros(NS,NS) ;
for i = 1:NS
    for j = 1:NS
        X1(i,j) = L1(i,1)*(1-t(j))+L2(i,1)*t(j) ;
        Y1(i,j) = L1(i,2)*(1-t(j))+L2(i,2)*t(j) ;
    end
end
%% Mesh region between square and circle boundary 
% Get boundaries of square 
x1 = [X1(1,:)' ; X1(:,end) ; flipud(X1(end,:)') ;flipud(X1(:,1))] ;
y1 = [Y1(1,:)' ; Y1(:,end) ; flipud(Y1(end,:)') ;flipud(Y1(:,1))] ;
% Delete double points if any
[x1, y1] = Eliminate([x1 y1],10^-5) ;
% Circle coordinates 
M = length(x1) ; 
th0 = linspace(0,2*pi,M) ;
% Arrange th
idx = find(th0>=3*pi/4) ;
th = [th0(idx) th0(2:idx(1))] ;
x2 = C(1)+R*cos(th') ;
y2 = C(2)+R*sin(th') ;
% Mesh the region using two boundary TFI
t = linspace(0,1,N) ;
X2 = zeros(M,N) ;
Y2 = zeros(M,N) ;
for i = 1:M
    for j = 1:N
        X2(i,j) = x1(i)*(1-t(j))+x2(i)*t(j) ;
        Y2(i,j) = y1(i)*(1-t(j))+y2(i)*t(j) ;
    end
end

X = [X1(:);X2(:)];
Y = [Y1(:);Y2(:)];
%% plot
%plotgrid(X1,Y1) ;
%plotgrid(X2,Y2) ;