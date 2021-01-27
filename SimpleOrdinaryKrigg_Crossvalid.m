%% Start
% This program performs a cross validation on Flow Unit 5 of the gross
% thickness data using both simple and ordinary kriging techniques
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


%% Simple Kriging Cross Validation

% Putting the X and Y coordinates into temporary variables
x = xcord;
y = ycord;
hh = h;
m = mean(h);

% Cross Validation
% Cross validation is performaed based on the leave-one-out method on all
% the sample points

for ii=1:numel(h)
    
    % Putting the sample point to be removed in defined variables to be
    % used later in distance calculatoins
    Xunsam = x(ii);  
    Yunsam = y(ii);
	Hunsam = hh(ii);
    
    % Removing the variable from the X, Y, and H matrices
    x(ii)=[];
    y(ii)=[];
	hh(ii)=[];
    
    % Putting the new sample points into newly defined variables
    X=x;
    Y=y;
	H=hh;
    for j=1:numel(X)
        for k=1:numel(X)
            
            %Distance calculation between all the sample points
            LD(j,k) = sqrt((((X(j) - X(k)))/au).^2 + (((Y(j) - Y(k)))/av).^2);
            
        end
        
        % Distance calculation betweeen the sample points and the removed
        % unsample point
        LD2(j,1) = sqrt((((X(j) - Xunsam))/au).^2 + (((Y(j) - Yunsam))/av).^2);
        
    end
    
    % Covariance between the sample points
    CovSamPoint = C0-(nugget + (sill*(1.5*LD - 0.5*(LD...
        ).^3).*(LD<=1) + sill*(LD>1)));
    
    % Covariance between the sample and unscample points
    CovSamUnsam = C0-(nugget + (sill*(1.5*LD2 - 0.5*(LD2...
        ).^3).*(LD2<=1) + sill*(LD2>1)));
    
    % Estimating the gross thickness by Simple Kriging Method
    lamda(:,ii) = CovSamPoint\ CovSamUnsam;
	lamda0(ii,1) = m*(1-sum(lamda(:,ii)));
	H_Krig(ii,1) = lamda0(ii,1)+sum(lamda(:,ii)'*H);
    
	% Reseting the X, Y coordinates and h to be perform the cross
	% validation on other points as well
    x = xcord;
    y = ycord;
    hh=h;
       
end

%% Ordinary Kriging Cross Validation

% Putting the X and Y coordinates into temporary variables
xo = xcord;
yo = ycord;
hho = h;

% Cross Validation
% Cross validation is performaed based on the leave-one-out method on all
% the sample points

for ii=1:numel(h)
    
    % Putting the sample point to be removed in defined variables to be
    % used later in distance calculatoins
    Xunsamo = xo(ii);  
    Yunsamo = yo(ii);
	Hunsamo = hho(ii);
    
    % Removing the variable from the X, Y, and H matrices
    xo(ii)=[];
    yo(ii)=[];
	hho(ii)=[];
    
    % Putting the new sample points into newly defined variables
    Xo=xo;
    Yo=yo;
	Ho=hho;
    for j=1:numel(Xo)
        for k=1:numel(Xo)
            
            %Distance calculation between all the sample points
            LDo(j,k) = sqrt((((Xo(j) - Xo(k)))/au).^2 + (((Yo(j) - Yo(k)))/av).^2);
            
        end
        
        % Distance calculation betweeen the sample points and the removed
        % unsample point
        LD2o(j,1) = sqrt((((Xo(j) - Xunsamo))/au).^2 + (((Yo(j) - Yunsamo))/av).^2);
        
    end
    
    % Covariance between the sample points
    CovSamPointo = C0-(nugget + (sill*(1.5*LDo - 0.5*(LDo...
        ).^3).*(LDo<=1) + sill*(LDo>1)));
    % Defining the left side covariance matrix for ordinary kriging
    n = length(Ho);
    CovSamPointo(:,n+1) = 1;
    CovSamPointo(n+1,:) = 1;
    CovSamPointo(n+1,n+1) = 0;
    
    % Covariance between the sample and unscample points
    CovSamUnsamo = C0-(nugget + (sill*(1.5*LD2o - 0.5*(LD2o...
        ).^3).*(LD2o<=1) + sill*(LD2o>1)));
    % Defining the left side covariance matrix for ordinary kriging
    CovSamUnsamo(n+1) = 1;
    
    % Estimating the gross thickness by Ordinary Kriging Method
    lamdao(:,ii) = CovSamPointo\ CovSamUnsamo;
    Lamdao=lamdao(:,ii);
    Lamdao(numel(h))=[];
	H_Krigo(ii,1) = sum(Lamdao'*Ho);
    
	% Reseting the X, Y coordinates and h to be perform the cross
	% validation on other points as well
    xo = xcord;
    yo = ycord;
    hho=h;
       
end


%% Visualization

% Simple Kriging
figure;
Error=H_Krig-h;
hist(Error,20);
title('\bfHistogram of Simple Kriging Cross-Validation for Flow Unit 5 gross thickness data','FontSize',10)
xlabel('\bfGross Thickness Error, Estimate Minus Data, ft','FontSize',10);
ylabel('\bfFrequency','FontSize',10);

figure;
R2simple = corrcoef(h,H_Krig);
R2simple = R2simple(2);
subplot(1,2,1)
plot(h,H_Krig,'ob');
hold on;
p = polyfit(h,H_Krig,1);
HKrigfit = polyval(p,h);
plot(h,HKrigfit,'k-','LineWidth',1.5);
xlim([0 30]);
ylim([0 30]);
hold off;
text(16,16,['\bfR^{2} = ' num2str(R2simple)]);
title('\bfSimple Kriging Cross-Validation','FontSize',10)
xlabel('\bfGross Thickness Data, ft','FontSize',10);
ylabel('\bfEstimated Gross Thickness, ft','FontSize',10);

subplot(1,2,2)
R2simpleerror = corrcoef(H_Krig,Error);
R2simpleerror = R2simpleerror(2);
plot(H_Krigo,Error,'ob');
hold on;
p1 = polyfit(H_Krig,Error,1);
Errorfit = polyval(p1,H_Krig);
plot(H_Krig,Errorfit,'k-','LineWidth',1.5);
xlim([0 30]);
ylim([-15 15]);
hold off;
text(15,0,['\bfR^{2} = ' num2str(R2simpleerror)]);
title('\bfSimple Kriging Cross-Validation Errors','FontSize',10)
xlabel('\bfSimple Kriging Estimate, ft','FontSize',10);
ylabel('\bfError, Estimate Minus Data','FontSize',10);

% Ordinary Kriging
figure;
Erroro=H_Krigo-h;
hist(Erroro,20);
title('\bfHistogram of Ordinary Kriging Cross-Validation for Flow Unit 5 gross thickness data','FontSize',10)
xlabel('\bfGross Thickness Error, Estimate Minus Data, ft','FontSize',10);
ylabel('\bfFrequency','FontSize',10);

figure;
R2ord = corrcoef(h,H_Krigo);
R2ord = R2ord(2);
subplot(1,2,1)
plot(h,H_Krigo,'ob');
hold on;
po = polyfit(h,H_Krigo,1);
HKrigfito = polyval(po,h);
plot(h,HKrigfito,'k-','LineWidth',1.5);
xlim([0 30]);
ylim([0 30]);
hold off;
text(16,16,['\bfR^{2} = ' num2str(R2ord)]);
title('\bfOrdinary Kriging Cross-Validation','FontSize',10)
xlabel('\bfGross Thickness Data, ft','FontSize',10);
ylabel('\bfEstimated Gross Thickness, ft','FontSize',10);

subplot(1,2,2)
R2orderror = corrcoef(H_Krigo,Erroro);
R2orderror = R2orderror(2);
plot(H_Krigo,Erroro,'ob');
hold on;
po1 = polyfit(H_Krigo,Erroro,1);
Errorofit = polyval(po1,H_Krigo);
plot(H_Krigo,Errorofit,'k-','LineWidth',1.5);
xlim([0 30]);
ylim([-15 15]);
hold off;
text(15,0,['\bfR^{2} = ' num2str(R2orderror)]);
title('\bfOrdinary Kriging Cross-Validation Errors','FontSize',10)
xlabel('\bfOrdinary Kriging Estimate, ft','FontSize',10);
ylabel('\bfError, Estimate Minus Data','FontSize',10);

% Ordinary vs. Simple
figure;
plot(H_Krig,H_Krigo,'ob');
title('\bfComparison of simple and ordinary kriging for Flow Unit 5 gross thickness data','FontSize',10)
xlabel('\bfSimple Kriging Estimate, ft','FontSize',10);
ylabel('\bfOrdinary Kriging Estimate, ft','FontSize',10);



























