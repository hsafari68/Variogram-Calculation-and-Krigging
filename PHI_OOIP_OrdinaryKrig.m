%% Start
% This program obtains a map based on Flow Unit 5 porosity data 
% using ordinary kriging and then calculates and obtains another map for
% OOIP
clc;
clear;
close all;

% Run the OrdinaryKrigingMap.m file first
load Zg;

%% Load and plot the porosity data

FU5phi=xlsread('Data.xlsx',2,'A2:C69');  % Load the porosity data
teta=0;

% Define variogram parameters
av=4000;
au=10000;
sill=16.25;
nugget=5;
C0=21.25;

phi=FU5phi(:,3);    % Porosity, %
xcord=FU5phi(:,1);  % X-coordinate, ft
ycord=FU5phi(:,2);  % Y-coordinate, ft

% Plot the porosity data on scatter plot
figure;
scatter(FU5phi(:,1),FU5phi(:,2));
m=num2str(FU5phi(:,3));
text(FU5phi(:,1),FU5phi(:,2),m);
grid on
xlim([-5000,25000]);
ylim([-10000,20000]);
title('\bfLocation Flow Unit 5 porosity data','FontSize',10)
xlabel('\bfX, ft','FontSize',10);
ylabel('\bfY, ft','FontSize',10);

%% Ordinary Kriging to Generate Porosity Map

x = xcord;
y = ycord;

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
phig = zeros(size(Xg));
Sigma2E = zeros(size(Xg));

for k = 1 : length(Xg)
    
    LD2 = sqrt((((x - Xg(k))*cos(teta)+(y - Yg(k))*sin(teta))/au).^2 +...
        ((-(x - Xg(k))*sin(teta)+(y - Yg(k))*cos(teta))/av).^2);
    
    %Defining the right side matrix of ordinary kriging
    CovSamUnsam = C0-(nugget + (sill*(1.5*LD2 - 0.5*(LD2...
        ).^3).*(LD2<=1) + sill*(LD2>1)));
    %Re-arranging the right side matrix
    CovSamUnsam(n+1) = 1;
    
    % Estimating the porosity by Ordinary Kriging Method
    lamda = CovSamPoint\ CovSamUnsam;
    phig(k) = sum(lamda(1:n,1).*phi);
    % Error Variance
    Sigma2E(k) = sill-(sum(lamda(1:n,1).*CovSamUnsam(1:n,1))-lamda(n+1,1));
end

%% OOIP Calculation
NTG=0.7;
Sw=0.2;
Bo=1.2;  % bbl/STB
A=blocksize*blocksize; % ft2
OOIP=(A*NTG*(1-Sw)*(phig/100).*Zg)/(5.615*Bo); % STB
Total_OOIP = sum(OOIP) % STB

%% Visualization

% Porosity map
figure;
r1 = length(R1);
r2 = length(R2);
PHI = reshape(phig,r2,r1);
h = pcolor(Xg1,Xg2,PHI);
set(h,'LineStyle','none')
axis equal
ylim([0 10000])
xlim([5000 15000])
title('\bfOrdinary kriging map for Flow Unit 5 porosity')
xlabel('\bfX, ft')
ylabel('\bfY, ft')
colormap(jet)
colorbar
% Variance map of porosity
figure;
r1 = length(R1);
r2 = length(R2);
SK = reshape(Sigma2E,r2,r1);
h = pcolor(Xg1,Xg2,SK);
set(h,'LineStyle','none')
axis equal
ylim([0 10000])
xlim([5000 15000])
title('\bfOrdinary kriging variance map for Flow Unit 5 porosity')
xlabel('\bfX, ft')
ylabel('\bfY, ft')
colormap(jet)
colorbar

% OOIP map
figure;
r1 = length(R1);
r2 = length(R2);
OOIP1 = reshape(OOIP,r2,r1);
h = pcolor(Xg1,Xg2,OOIP1);
set(h,'LineStyle','none')
axis equal
ylim([0 10000])
xlim([5000 15000])
title('\bfOrdinary kriging map for Flow Unit 5 OOIP, STB')
xlabel('\bfX, ft')
ylabel('\bfY, ft')
colormap(jet)
colorbar
