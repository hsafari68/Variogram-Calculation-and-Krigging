%% Start
% This program obtains a map based on Flow Unit 5 gross thickness data 
% using ordinary kriging
clc;
clear;
close all;

%% Load and plot the gross thickness data

FU5h=xlsread('Data.xlsx',1,'A2:C104');  % Load the gross thickness data
teta=0;

% Define variogram parameters
av=1400;
au=6000;
sill=16.25;
nugget=5;
C0=21.25;


h=FU5h(:,3);      % Thickness, ft
xcord=FU5h(:,1);  % X-coordinate, ft
ycord=FU5h(:,2);  % Y-coordinate, ft

% Plot the gross thickness data on scatter plot
figure;
scatter(FU5h(:,1),FU5h(:,2));
m=num2str(FU5h(:,3));
text(FU5h(:,1),FU5h(:,2),m);
grid on
xlim([-5000,25000]);
ylim([-10000,20000]);
title('\bfLocation Flow Unit 5 gross thickness data','FontSize',10)
xlabel('\bfX, ft','FontSize',10);
ylabel('\bfY, ft','FontSize',10);

%% Ordinary Kriging to Generate Gross Thickness Maps

x = xcord;
y = ycord;
z = h;
blocksize = 200;

% Discretizing the horizontal and vertical directions
[X1,X2] = meshgrid(x);
[Y1,Y2] = meshgrid(y);
LD = sqrt((((X1 - X2)*cos(teta)+(Y1 - Y2)*sin(teta))/au).^2 +...
((-(X1 - X2)*sin(teta)+(Y1 - Y2)*cos(teta))/av).^2);

% Defining the left side matrix of ordinary kriging
CovSamPoint = C0-(nugget + (sill*(1.5*LD - 0.5*(LD...
).^3).*(LD<=1) + sill*(LD>1)));
%Re-arranging the left side matrix
n = length(x);
CovSamPoint(:,n+1) = 1;
CovSamPoint(n+1,:) = 1;
CovSamPoint(n+1,n+1) = 0;

R1 = 5000 : blocksize : 15000-1;
R2 = 0 : blocksize : 10000-1;
[Xg1,Xg2] = meshgrid(R1,R2);
Xg = reshape(Xg1,[],1);
Yg = reshape(Xg2,[],1);
Zg = zeros(size(Xg));
Sigma2E = zeros(size(Xg));

for k = 1 : length(Xg)
    
    LD2 = sqrt((((x - Xg(k))*cos(teta)+(y - Yg(k))*sin(teta))/au).^2 +...
        ((-(x - Xg(k))*sin(teta)+(y - Yg(k))*cos(teta))/av).^2);
    
    %Defining the right side matrix of ordinary kriging
    CovSamUnsam = C0-(nugget + (sill*(1.5*LD2 - 0.5*(LD2...
        ).^3).*(LD2<=1) + sill*(LD2>1)));
    %Re-arranging the right side matrix
    CovSamUnsam(n+1) = 1;
    
    % Estimating the gross thickness by Ordinary Kriging Method
    lamda = CovSamPoint\ CovSamUnsam;
    Zg(k) = sum(lamda(1:n,1).*z);
    % Error Variance
    Sigma2E(k) = sill-(sum(lamda(1:n,1).*CovSamUnsam(1:n,1))-lamda(n+1,1));
end

%% Visualization

% Gross thickness map
figure;
r1 = length(R1);
r2 = length(R2);
Z = reshape(Zg,r2,r1);
h = pcolor(Xg1,Xg2,Z);
set(h,'LineStyle','none')
axis equal
ylim([0 10000])
xlim([5000 15000])
title('\bfOrdinary kriging map for Flow Unit 5 Gross thickness')
xlabel('\bfX, ft')
ylabel('\bfY, ft')
colormap(jet)
colorbar
% Variance map of gross thickness
figure;
r1 = length(R1);
r2 = length(R2);
SK = reshape(Sigma2E,r2,r1);
h = pcolor(Xg1,Xg2,SK);
set(h,'LineStyle','none')
axis equal
ylim([0 10000])
xlim([5000 15000])
title('\bfOrdinary kriging variance map for Flow Unit 5 Gross thickness')
xlabel('\bfX, ft')
ylabel('\bfY, ft')
colormap(jet)
colorbar

save Zg