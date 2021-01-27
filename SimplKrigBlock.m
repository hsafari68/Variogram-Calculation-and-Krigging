%% Start
% This program estimates block map of the flow unit 5 thickness data using
% simple kriging method
clc;
clear;
close all;

%% Load and plot the gross thickness data

FU5h=xlsread('Data.xlsx',1,'A2:C104');  % Load the gross thickness data
teta=0;

% Define variogram parameters
au=1400;
av=6000;
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

%% Simple Kriging for Block Estimation
x = xcord;
y = ycord;
z = h;
blocksize = 200;

% Mesh Grid Definition
[X1,X2] = meshgrid(x);
[Y1,Y2] = meshgrid(y);
LD = sqrt((((X1 - X2)*cos(teta)+(Y1 - Y2)*sin(teta))/au).^2 +...
((-(X1 - X2)*sin(teta)+(Y1 - Y2)*cos(teta))/av).^2);

CovSamPoint = C0-(nugget + (sill*(1.5*LD - 0.5*(LD...
).^3).*(LD<=1) + sill*(LD>1)));

R1 = 5000 : blocksize : 15000-1;
R2 = 0 : blocksize : 10000-1;
[Xg1,Xg2] = meshgrid(R1,R2);
Xg = reshape(Xg1,[],1);
Yg = reshape(Xg2,[],1);

% Block Kriging
x = xcord;
y = ycord;
InnerBlockSize = 50; % Defining 16 points inside each block

for ii=1:length(Xg)
    
    R11 = Xg(ii) : InnerBlockSize : Xg(ii)+blocksize-1;
    R22 = Yg(ii) : InnerBlockSize : Yg(ii)+blocksize-1;
    [XgB1,XgB2] = meshgrid(R11,R22);
    XgB = reshape(XgB1,[],1);
    YgB = reshape(XgB2,[],1);
    
    
    for jj=1:length(XgB)
        
        LD2(:,jj) = sqrt((((x - XgB(jj))*cos(teta)+(y - YgB(jj))*sin(teta))/au).^2 +...
            ((-(x - XgB(jj))*sin(teta)+(y - YgB(jj))*cos(teta))/av).^2);
        CovSamUnsam = C0-(nugget + (sill*(1.5*LD2 - 0.5*(LD2...
            ).^3).*(LD2<=1) + sill*(LD2>1)));
        
        % Calculating the mean of covariance between 16 unsample points
        % and all other sample points
        CovSamUnsam = (1/numel(XgB))*sum(CovSamUnsam,2);
        
        
        lamda(:,ii) = CovSamPoint\ CovSamUnsam;
        mean=sum(h)/numel(h);
        lamda0(:,ii)=mean*(1-sum(lamda(:,ii)));
        H_Krig(ii,1) = lamda0(:,ii)+sum(lamda(:,ii)'*h);
        Sigma2E(ii,1) = sill-(sum(lamda(:,ii)'*CovSamUnsam));

        
    end
 
end



%% Visualization

% Gross Thickness Map
figure;
r1 = length(R1);
r2 = length(R2);
Z = reshape(H_Krig,r2,r1);
hh = pcolor(Xg1,Xg2,Z);
set(hh,'LineStyle','none')
axis equal
ylim([0 10000])
xlim([5000 15000])
title('\bfSimple block kriging map for Flow Unit 5 gross thickness')
xlabel('\bfX, ft')
ylabel('\bfY, ft')
colormap(jet)
colorbar

% Error variance map of gross thickness
figure;
r1 = length(R1);
r2 = length(R2);
SK = reshape(Sigma2E,r2,r1);
hh = pcolor(Xg1,Xg2,SK);
set(hh,'LineStyle','none')
axis equal
ylim([0 10000])
xlim([5000 15000])
title('\bfSimple block kriging variance map for Flow Unit 5 gross thickness')
xlabel('\bfX, ft')
ylabel('\bfY, ft')
colormap(jet)
colorbar





