clear; close all; clc;
addpath 'D:\git-repository\SkiPR\functions'
addpath 'E:\MCS'
%% Load Ancillary Data Files from LiDAR_GPR.mat and snowMachine.mat
dataDir = 'E:\MCS\MCS031722';
% Coordinate Data
load([dataDir,'\LiDAR\','MCS031722-Coords.mat'])
% LiDAR Data
load([dataDir,'\LiDAR\','MCS031722-LiDAR.mat'])
% LiDAR - GPR data
load([dataDir,'\GPR\','MCS031722-LiDAR-GPR.mat'])
% kdTree
load('MCS03172022kdtree.mat')
% Machine Learning
load('MCS03172022ML.mat')
%% Calulate SWE
% MLR
ML.MLR.spatialSWE = ML.MLR.spatialDensity.*LiDAR.A.A; % mm w.e.
% ANN
ML.ANN.ensemble.spatialSWE = ML.ANN.ensemble.avgSpatialDensity.*LiDAR.A.A; % mm w.e.
% RF
ML.RF.spatialSWE = ML.RF.spatialDensity.*LiDAR.A.A;
%% Comparisons with In-situ Observations
% Brought in in-situ data
% E:\MCS\MCS031722\inSitu
% Snow Pit Data
inSitu.snowPit.X = 604399; inSitu.snowPit.Y = 4868503.5;
inSitu.snowPit.Density = [180	167; 235	246; 315	339;427	423;351	361;352	372;...
379	384;376	391;381	406;401	429;425	416;440	438;423	432;425	390;391	378];
inSitu.snowPit.avgDensity = mean(mean(inSitu.snowPit.Density,2));
inSitu.snowPit.rmseDensity = sqrt(mean((diff(inSitu.snowPit.Density,1,2)).^2));
inSitu.snowPit.reDensity = (inSitu.snowPit.rmseDensity./inSitu.snowPit.avgDensity).*100;
inSitu.snowPit.Height = [160;150;140;130;120;110;100;90;80;70;60;50;40;30;20];
inSitu.snowPit.Depth = max(inSitu.snowPit.Height);
inSitu.snowPit.intervalHeight = 10.*ones(length(inSitu.snowPit.Density),1);
% inSitu.snowPit.intervalSWE = (mean(inSitu.snowPit.Density,2).*inSitu.snowPit.intervalHeight)./10; % mm swe
inSitu.snowPit.avgSWE = (inSitu.snowPit.avgDensity.*inSitu.snowPit.Depth)./100;

% Depth Probe Data
dataDir = 'E:\MCS\MCS031722\inSitu';
filename = 'MCS031722ProbeDepthsDensityT8corrected.csv';
inSitu.probe.data = readtable([dataDir,'\',filename]);
[inSitu.probe.data.X, inSitu.probe.data.Y] = deg2utm(inSitu.probe.data.Latitude,inSitu.probe.data.Longitude);
inSitu.probe.data.SWE = (inSitu.probe.data.Depth.*inSitu.probe.data.Density)./100;
%% kdtree searcher for Depths
% KD-tree Searcher
isKDtree = 0;
if isKDtree
tic
% winsize = .25;
winsize = mean([Coords.R.CellExtentInWorldX,Coords.R.CellExtentInWorldY]);
% mykdtree=KDTreeSearcher([Coords.Xi, Coords.Yi]); % searcher for GPR
[IDX,D]=rangesearch(mykdtree,[inSitu.probe.data.X inSitu.probe.data.Y],winsize); % search the points in the LiDAR grid
% Remove Empty Cells
ix =  find(~cellfun(@isempty,D));
D = D(ix); IDX = IDX(ix);
kdProbe.D = D; kdProbe.IDX = IDX;kdProbe.ix = ix;kdProbe.winsize = winsize;
toc
% Save the Output
save('MCS03172022kdtreeProbe.mat','kdProbe','-v7.3')
clear('D','IDX','winsize','ix')
else
% Load Previous KD-Tree
load('MCS03172022kdtreeProbe.mat')
end

%% Take the median
inSitu.probe.lidar.Depth = zeros(length(kdProbe.IDX),1);
inSitu.probe.lidar.X = zeros(length(kdProbe.IDX),1); 
inSitu.probe.lidar.Y = zeros(length(kdProbe.IDX),1);
for jj = 1:length(kdProbe.D)
    % Extract Binned Indicies from kdTree
    tmpIx = kdProbe.IDX{jj};tmpIx = tmpIx(:);
    % Median
    tmpDepth = median(LiDAR.A.A(tmpIx));
    [~,tmpDepthIx] = min(abs(tmpDepth-LiDAR.A.A(tmpIx)));
    inSitu.probe.lidar.Depth(jj) = tmpDepth;
    % Coordinates
    inSitu.probe.lidar.X(jj) = Coords.Xi(kdProbe.IDX{jj}(tmpDepthIx));
    inSitu.probe.lidar.Y(jj) = Coords.Yi(kdProbe.IDX{jj}(tmpDepthIx));
end

%% Locate Snow Pit
[~ , pitIx] = mink(sqrt((inSitu.snowPit.X-Coords.Xi).^2+(inSitu.snowPit.Y-Coords.Yi).^2),5);
inSitu.snowPit.valDepth = median(LiDAR.A.A(pitIx));
%% Depth Validation
depthVals = [inSitu.probe.lidar.Depth;inSitu.snowPit.valDepth].*100;
depthObs = [inSitu.probe.data.Depth;inSitu.snowPit.Depth];
inSitu.validation.probe.insituDepth = depthObs;
inSitu.validation.Depth = depthObs;
inSitu.validation.probe.estDepth = depthVals;
inSitu.validation.probe.R2 = corr(depthObs,depthVals).^2;
inSitu.validation.probe.ME = mean(depthObs-depthVals);
inSitu.validation.probe.MAE = mean(abs(depthObs-depthVals));
inSitu.validation.probe.RMSE = sqrt(mean(((depthObs-depthVals).^2)));
% Figure
figure(); plot(depthObs,depthVals,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([80,240],[80,240],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (cm)'); ylabel('LiDAR (cm)'); title('Snow Depth')
xlim([80,240]);ylim([80,240]); grid on ; grid minor; axis square
annotation('textbox',[0.3,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.probe.R2,2)),'   RMSE = ',num2str(round(inSitu.validation.probe.RMSE)),'   Bias = ',num2str(round(inSitu.validation.probe.ME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

%% Density & SWE Validation
densityIx = find(~isnan(inSitu.probe.data.Density));
inSitu.validation.X = [inSitu.probe.data.X(densityIx);inSitu.snowPit.X];
inSitu.validation.Y = [inSitu.probe.data.Y(densityIx);inSitu.snowPit.Y];
inSitu.validation.Density = [inSitu.probe.data.Density(densityIx);inSitu.snowPit.avgDensity];
% inSitu.validation.SWE = [inSitu.probe.data.Density(densityIx);inSitu.snowPit.avgDensity].*...
%     [inSitu.probe.data.Depth(densityIx);inSitu.snowPit.Depth]./100;
inSitu.validation.SWE = [inSitu.probe.data.SWE(densityIx);inSitu.snowPit.avgSWE];
% GPR Validation
inSitu.validation.GPR.estDensity = zeros(numel(inSitu.validation.Density),1);
inSitu.validation.GPR.estSWE = zeros(numel(inSitu.validation.SWE),1);
for kk = 1:length(inSitu.validation.Density)
    [~, tmpValIx] = mink(sqrt((lidarGPR.X-inSitu.validation.X(kk)).^2+(lidarGPR.Y - inSitu.validation.Y(kk)).^2),10);
    inSitu.validation.GPR.estDensity(kk) = median(lidarGPR.Density(tmpValIx));
    inSitu.validation.GPR.estSWE(kk) = median(lidarGPR.SWE(tmpValIx));
%     inSitu.validation.GPR.estSWE(kk) = median(lidarGPR.Density(tmpValIx).*lidarGPR.depthGPR(tmpValIx))./100;

end
% Density
inSitu.validation.GPR.densityR2 = corr(inSitu.validation.GPR.estDensity,inSitu.validation.Density).^2;
inSitu.validation.GPR.densityME = mean(inSitu.validation.GPR.estDensity-inSitu.validation.Density);
inSitu.validation.GPR.densityMAE = mean(abs(inSitu.validation.GPR.estDensity-inSitu.validation.Density));
inSitu.validation.GPR.densityRMSE = sqrt(mean(((inSitu.validation.GPR.estDensity-inSitu.validation.Density).^2)));
% SWE
inSitu.validation.GPR.sweR2 = corr(inSitu.validation.GPR.estSWE,inSitu.validation.SWE).^2;
inSitu.validation.GPR.sweME = mean(inSitu.validation.GPR.estSWE-inSitu.validation.SWE);
inSitu.validation.GPR.sweMAE = mean(abs(inSitu.validation.GPR.estSWE-inSitu.validation.SWE));
inSitu.validation.GPR.sweRMSE = sqrt(mean((inSitu.validation.GPR.estSWE-inSitu.validation.SWE).^2));
% Figure
% Density
figure(); plot(inSitu.validation.Density,inSitu.validation.GPR.estDensity,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([300,500],[300,500],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (kg/m^3)'); ylabel('LiDAR -- GPR (kg/m^3)'); title('Density')
xlim([300,500]);ylim([300,500]); grid on ; grid minor; axis square
% annotation('textbox',[0.3,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.GPR.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.GPR.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.GPR.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')
annotation('textbox',[0.235,0.7,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.GPR.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.GPR.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.GPR.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

% SWE
figure(); plot(inSitu.validation.SWE,inSitu.validation.GPR.estSWE,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([350,750],[350,750],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (mm)'); ylabel('LiDAR -- GPR (mm)'); title('Snow Water Equivalent')
xlim([350,750]);ylim([350,750]); grid on ; grid minor; axis square
annotation('textbox',[0.35,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.GPR.sweR2,2)),'   RMSE = ',num2str(round(inSitu.validation.GPR.sweRMSE)),'   Bias = ',num2str(round(inSitu.validation.GPR.sweME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

% Machine Learning Validation
inSitu.validation.ML.MLR.estDensity = zeros(numel(inSitu.validation.Density),1);
inSitu.validation.ML.MLR.estSWE = zeros(numel(inSitu.validation.SWE),1);
inSitu.validation.ML.ANN.estDensity = zeros(numel(inSitu.validation.Density),1);
inSitu.validation.ML.ANN.estSWE = zeros(numel(inSitu.validation.SWE),1);
inSitu.validation.ML.RF.estDensity = zeros(numel(inSitu.validation.Density),1);
inSitu.validation.ML.RF.estSWE = zeros(numel(inSitu.validation.SWE),1);

for kk = 1:length(inSitu.validation.Density)
        [~, tmpValIx] = mink(sqrt((Coords.Xi-inSitu.validation.X(kk)).^2+(Coords.Yi - inSitu.validation.Y(kk)).^2),10);
        inSitu.validation.ML.MLR.estDensity(kk) = median(ML.MLR.spatialDensity(tmpValIx));
        inSitu.validation.ML.MLR.estSWE(kk) = median(ML.MLR.spatialSWE(tmpValIx));
        inSitu.validation.ML.ANN.estDensity(kk) = median(ML.ANN.ensemble.avgSpatialDensity(tmpValIx));
        inSitu.validation.ML.ANN.estSWE(kk) = median(ML.ANN.ensemble.spatialSWE(tmpValIx));
        inSitu.validation.ML.RF.estDensity(kk) = median(ML.RF.spatialDensity(tmpValIx));
        inSitu.validation.ML.RF.estSWE(kk) = median(ML.RF.spatialSWE(tmpValIx));
end
% MLR
% Density
inSitu.validation.ML.MLR.densityR2 = corr(inSitu.validation.ML.MLR.estDensity,inSitu.validation.Density).^2;
inSitu.validation.ML.MLR.densityME = mean(inSitu.validation.ML.MLR.estDensity-inSitu.validation.Density);
inSitu.validation.ML.MLR.densityMAE = mean(abs(inSitu.validation.ML.MLR.estDensity-inSitu.validation.Density));
inSitu.validation.ML.MLR.densityRMSE = sqrt(mean(((inSitu.validation.ML.MLR.estDensity-inSitu.validation.Density).^2)));
% SWE
inSitu.validation.ML.MLR.sweR2 = corr(inSitu.validation.ML.MLR.estSWE,inSitu.validation.SWE).^2;
inSitu.validation.ML.MLR.sweME = mean(inSitu.validation.ML.MLR.estSWE-inSitu.validation.SWE);
inSitu.validation.ML.MLR.sweMAE = mean(abs(inSitu.validation.ML.MLR.estSWE-inSitu.validation.SWE));
inSitu.validation.ML.MLR.sweRMSE = sqrt(mean((inSitu.validation.ML.MLR.estSWE-inSitu.validation.SWE).^2));
% Figure
% Density
figure(); plot(inSitu.validation.Density,inSitu.validation.ML.MLR.estDensity,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([300,500],[300,500],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (kg/m^3)'); ylabel('MLR (kg/m^3)'); title('Density')
xlim([300,500]);ylim([300,500]); grid on ; grid minor; axis square
% annotation('textbox',[0.235,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.MLR.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.MLR.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.MLR.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')
annotation('textbox',[0.235,0.7,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.MLR.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.MLR.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.MLR.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

% SWE
figure(); plot(inSitu.validation.SWE,inSitu.validation.GPR.estSWE,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([350,750],[350,750],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (mm)'); ylabel('MLR (mm)'); title('Snow Water Equivalent')
xlim([350,750]);ylim([350,750]); grid on ; grid minor; axis square
annotation('textbox',[0.235,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.MLR.sweR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.MLR.sweRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.MLR.sweME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

% ANN
% Density
inSitu.validation.ML.ANN.densityR2 = corr(inSitu.validation.ML.ANN.estDensity,inSitu.validation.Density).^2;
inSitu.validation.ML.ANN.densityME = mean(inSitu.validation.ML.ANN.estDensity-inSitu.validation.Density);
inSitu.validation.ML.ANN.densityMAE = mean(abs(inSitu.validation.ML.ANN.estDensity-inSitu.validation.Density));
inSitu.validation.ML.ANN.densityRMSE = sqrt(mean(((inSitu.validation.ML.ANN.estDensity-inSitu.validation.Density).^2)));
% SWE
inSitu.validation.ML.ANN.sweR2 = corr(inSitu.validation.ML.ANN.estSWE,inSitu.validation.SWE).^2;
inSitu.validation.ML.ANN.sweME = mean(inSitu.validation.ML.ANN.estSWE-inSitu.validation.SWE);
inSitu.validation.ML.ANN.sweMAE = mean(abs(inSitu.validation.ML.ANN.estSWE-inSitu.validation.SWE));
inSitu.validation.ML.ANN.sweRMSE = sqrt(mean((inSitu.validation.ML.ANN.estSWE-inSitu.validation.SWE).^2));
% Figure
% Density
figure(); plot(inSitu.validation.Density,inSitu.validation.ML.ANN.estDensity,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([300,500],[300,500],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (kg/m^3)'); ylabel('ANN (kg/m^3)'); title('Density')
xlim([300,500]);ylim([300,500]); grid on ; grid minor; axis square
% annotation('textbox',[0.235,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.ANN.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.ANN.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.ANN.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')
annotation('textbox',[0.235,0.7,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.ANN.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.ANN.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.ANN.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

% SWE
figure(); plot(inSitu.validation.SWE,inSitu.validation.GPR.estSWE,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([350,750],[350,750],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (mm)'); ylabel('ANN (mm)'); title('Snow Water Equivalent')
xlim([350,750]);ylim([350,750]); grid on ; grid minor; axis square
annotation('textbox',[0.235,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.ANN.sweR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.ANN.sweRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.ANN.sweME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

% RF
% Density
inSitu.validation.ML.RF.densityR2 = corr(inSitu.validation.ML.RF.estDensity,inSitu.validation.Density).^2;
inSitu.validation.ML.RF.densityME = mean(inSitu.validation.ML.RF.estDensity-inSitu.validation.Density);
inSitu.validation.ML.RF.densityMAE = mean(abs(inSitu.validation.ML.RF.estDensity-inSitu.validation.Density));
inSitu.validation.ML.RF.densityRMSE = sqrt(mean(((inSitu.validation.ML.RF.estDensity-inSitu.validation.Density).^2)));
% SWE
inSitu.validation.ML.RF.sweR2 = corr(inSitu.validation.ML.RF.estSWE,inSitu.validation.SWE).^2;
inSitu.validation.ML.RF.sweME = mean(inSitu.validation.ML.RF.estSWE-inSitu.validation.SWE);
inSitu.validation.ML.RF.sweMAE = mean(abs(inSitu.validation.ML.RF.estSWE-inSitu.validation.SWE));
inSitu.validation.ML.RF.sweRMSE = sqrt(mean((inSitu.validation.ML.RF.estSWE-inSitu.validation.SWE).^2));
% Figure
% Density
figure(); plot(inSitu.validation.Density,inSitu.validation.ML.RF.estDensity,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([300,500],[300,500],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (kg/m^3)'); ylabel('RF (kg/m^3)'); title('Density')
xlim([300,500]);ylim([300,500]); grid on ; grid minor; axis square
% annotation('textbox',[0.235,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.RF.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.RF.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.RF.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')
annotation('textbox',[0.235,0.7,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.RF.densityR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.RF.densityRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.RF.densityME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

% SWE
figure(); plot(inSitu.validation.SWE,inSitu.validation.GPR.estSWE,'ok','MarkerSize',5,'MarkerFaceColor','k');
hold on; plot([350,750],[350,750],'--k','LineWidth',1.5)
l = lsline; l.LineWidth = 2;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlabel('In-Situ (mm)'); ylabel('RF (mm)'); title('Snow Water Equivalent')
xlim([350,750]);ylim([350,750]); grid on ; grid minor; axis square
annotation('textbox',[0.235,0.05,.2,.2],'string',['R^2 = ',num2str(round(inSitu.validation.ML.RF.sweR2,2)),'   RMSE = ',num2str(round(inSitu.validation.ML.RF.sweRMSE)),'   Bias = ',num2str(round(inSitu.validation.ML.RF.sweME))],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',14,'fontweight','bold')

%% Total SWE
stupidSWE = mean(inSitu.validation.Density).*LiDAR.A.A;
randIx = datasample([1:numel(ML.MLR.spatialDensity)]',round(0.01.*numel(ML.MLR.spatialDensity)),'Replace',false);
figure();
boxplot([ML.MLR.spatialSWE(randIx),ML.ANN.ensemble.spatialSWE(randIx),ML.RF.spatialSWE(randIx),stupidSWE(randIx)],'Notch','off')
% Change the boxplot color
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
% a1 = a(1:3);  
% set(a1, 'marker', '.','markeredgeColor',[0,0,0],'markersize',3); 
% a2 = a(3:6);% %median
% set(a2,'linewidth',3,'color',[0.5,0.5,0.5])
% a3 = a(7:21);%box
a1 = a(1:4);  
set(a1, 'marker', '.','markeredgeColor',[0,0,0],'markersize',3); 
a2 = a(5:8);% %median
set(a2,'linewidth',3,'color',[0.5,0.5,0.5])
a3 = a(9:28);%box
set(a3,'Color','k','linewidth',2)
grid on; grid minor; boxLabels = {'MLR','ANN','RF','Mean \rho'};
set(gca,'fontname','serif','fontweight','bold','fontsize',12,'XTickLabel',boxLabels,'TickLabelInterpreter','tex')
ylabel('Snow Water Equivalent (mm)'); title('Model Comparison')
ylim([0,1100])

bar([ML.RF.OOBPermutedPredictorRelativeImportance],'facecolor',[0.5,0.5,0.5],'edgecolor','k');
ylabel(['Relative Importance (%)'],'fontweight','bold');xlabel('Predictor','fontweight','bold');
paramNames3 = {' ','Hs','aspctHs','sxHs','dyHs','dxHs',...
    'Zs','aspctZs','sxZs','dyZs','dxZs','Zg','aspctZg','sxZg','dyZg','dxZg','Hveg','Sveg','NBR',' '};
set(gca,'xtick',[0:size(coms,2)+1],'xticklabel',paramNames3,'xticklabelrotation',45,'fontname','serif')
xlim([0,numel(ML.ANN.R)+1]); 
grid minor;grid on;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
%% Transect and Raster Validation Figures
figure();
for kk = 1%:3
figure();
imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.C.northness);colormap(bone)
freezeColors;
hold on;
if kk == 1
% MLR
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLR.spatialDensity,'AlphaData',0.65); 
end
if kk == 2
% ANN
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN.ensemble.avgSpatialDensity); 
end
if kk == 3
% RF
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.RF.spatialDensity); 
end
hold on; scatter(inSitu.validation.X./1000,inSitu.validation.Y./1000,100,inSitu.validation.Density,'filled','markeredgecolor','k')
daspect([1,1,1]);set(gca,'YDir','normal');
cb = colorbar;cb.FontSize = 14;
% cbfreeze(cb);
colormap(cmap); 
caxis([300,450]);
%caxis([quantile(spatialDensity(:),[0.025,0.975])]);
% xlabel('Easting (km)'); ylabel('Northing (km)');
% cb.Label.String = 'Average Snow Density (kg/m^3)'; 
if kk == 1
cb.Label.String = 'MLR Average Snow Density (kg/m^3)'; 
elseif kk == 2
cb.Label.String = 'ANN Average Snow Density (kg/m^3)'; 
elseif kk == 3
cb.Label.String = 'RF Average Snow Density (kg/m^3)';
end
%title('MLR Modeled Average Snow Density (kg/m^3)')
% xlabel('Easting (m)'); ylabel('Northing (m)');title('MLR Modeled Average Snow Density (kg/m^3)')

set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
% hold on; plot(Xs(snowPitIx),Ys(snowPitIx),'ok','markerfacecolor','k')
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'position',pos)
xlim([604,605.5])   
ylim([4868.0,4868.72])
xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
% hold on; plot(lidarGPR.X./1000,lidarGPR.Y./1000,'ok','markersize',5,'markerfacecolor','k')
end
%% GPR Transects Overlaid on Northness with Validation Density
figure();
imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(LiDAR.B.aspect+45)+sind(LiDAR.B.aspect+45))./2);colormap(bone)
freezeColors;
hold on;
% imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.D.D,"AlphaData",0.5,[0,10]);colormap(bone)
% freezeColors;
% imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.F.F,"AlphaData",1,[-500,500]);colormap(bone)
% freezeColors;
% hold on;
% imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45))./2,'AlphaData',1);colormap(bone)
% freezeColors;
scatter(lidarGPR.X./1000,lidarGPR.Y./1000,25,lidarGPR.Density,'filled')
hold on; scatter(inSitu.validation.X./1000,inSitu.validation.Y./1000,75,inSitu.validation.Density,'filled','markeredgecolor','k')
daspect([1,1,1]);set(gca,'YDir','normal');
cb = colorbar;cb.FontSize = 14;
% cbfreeze(cb);
colormap(cmap); 
caxis([300,450]);
%caxis([quantile(spatialDensity(:),[0.025,0.975])]);
% xlabel('Easting (km)'); ylabel('Northing (km)');
cb.Label.String = 'Average Snow Density (kg/m^3)'; 
%title('MLR Modeled Average Snow Density (kg/m^3)')
% xlabel('Easting (m)'); ylabel('Northing (m)');title('MLR Modeled Average Snow Density (kg/m^3)')

set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
% hold on; plot(Xs(snowPitIx),Ys(snowPitIx),'ok','markerfacecolor','k')
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'position',pos)
xlim([604.25,605.5])   
ylim([4868.0,4868.725])
xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
% hold on; plot(lidarGPR.X./1000,lidarGPR.Y./1000,'ok','markersize',5,'markerfacecolor','k')