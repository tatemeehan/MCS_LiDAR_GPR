%% Magneprobe
clear; close all; clc
addpath(genpath('D:\git-repository\SkiPR'))
addpath 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\UAVSAR-ASO-compare';
addpath 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\LIDAR'
cmap = csvread('D:\git-repository\GreenTrACS_MxRadar\colorMaps\RdYlBu.csv');
cmap = flipud(cmap);
%% Read data
dataDir = 'E:\MCS\MCS122123\Depth';
f1 = 'MCS122123_DEPTH1.csv';
f2 = 'MCS122123_DEPTH2.csv';
% Options
% opts = detectImportOptions([dataDir,'\',f1]);
% opts.VariableNamesLine = 2;
% opts.DataLines = [5 inf];

D1 = readtable([dataDir,'\',f1]);
D2 = readtable([dataDir,'\',f2]);
% Format LLZ
isS = strcmp(D1.ggan_s_ind,'S');
isW = strcmp(D1.ggae_w_ind,'W');
D1.ggailatitude = D1.ggailatitude./100;
D1.ggailatitude(isS) = D1.ggailatitude(isS).*-1;
D1.ggalongitude = D1.ggalongitude./100;
D1.ggalongitude(isW) = D1.ggalongitude(isW).*-1;
% Create LLZdepth array
Depth = [D1.ggalongitude,D1.ggailatitude,D1.ggaaltitude,D1.DepthCm./100;...
    D2.longitude_deg_,D2.latitude_deg_,D2.altitudeB,D2.DepthCm./100];
% Convert 2 utm
[Xp,Yp,zoneP] = deg2utm(Depth(:,2),Depth(:,1));
Depth(:,1) = Xp; Depth(:,2) = Yp;

% Mores Creek Snotel
snotelLat =  43.932067;
snotelLon = -115.665991;
[snotelX, snotelY] = deg2utm(snotelLat,snotelLon);
% snotelDepthChange = 0.05; % (m)
snotelDepth = 17.*2.54./100;
Depth = [Depth;repmat([snotelX,snotelY,0,snotelDepth],100,1)];
%% Read Gapfilled LiDAR Depths
% Data Directory 
dataDir = 'E:\MCS\MCS122123\LiDAR';
% Write Directory
writeDir = dataDir;
% A Snow Depth
filename = 'MCS_20231228_SNOWDEPTH_RFgapfilled.tif';
outfnA = filename(1:end-4);
fullfilename = fullfile(dataDir,filename);
[A,RA,~,~,lon,lat,utmX,utmY] = readLidarTif(fullfilename);
% Get UTM Coordinates as Vector
X = utmX(1,:);
Y = utmY(:,1);
Xi = utmX(:); 
Yi = utmY(:);
%% Co-Locate LiDAR and Depth Probe
% KD-tree Searcher
isKDtree = 0;
if isKDtree
tic
winsize = mean([RA.CellExtentInWorldX,RA.CellExtentInWorldY]);
mykdtree=KDTreeSearcher([Depth(:,1) Depth(:,2)]); % searcher for Probe
[IDX,D]=rangesearch(mykdtree,[Xi, Yi],winsize); % search the points in the LiDAR grid
% Remove Empty Cells
ix =  find(~cellfun(@isempty,D));
D = D(ix); IDX = IDX(ix);
kd.D = D; kd.IDX = IDX;kd.ix = ix;kd.winsize = winsize;
toc
% Save the Output
save([dataDir,outfnA(1:23),'kdtree2.mat'],'kd','-v7.3')

else
% Load Previous KD-Tree
load([dataDir,outfnA(1:23),'kdtree2.mat'])
end
%% Calculate the Error
dp = zeros(length(kd.IDX),1);
for kk = 1:length(kd.IDX)
    dp(kk) = mean([Depth(kd.IDX{kk},4)]);
end
% Error
depthError = dp-A(kd.ix);
%% Load PostPro LiDAR
    % Coordinate Data
    load([dataDir,'\',outfnA(1:13),'Coords.mat'])
    % LiDAR Data
    load([dataDir,'\',outfnA(1:13),'LiDAR.mat'])
%% Regress the Error
% Standardize Data
ixOG = kd.ix;
predictIx = find(~isnan(LiDAR.A.MLR));
% ixOG = find(~mask);
% NaN Mask
load("MCSborder.mat")
mask = isnan(LiDAR.A.MLR)+border;
mask = mask>0;
predictIx = find(~mask);
p = .0001;
nMC = 2500;
predictorMean = zeros(nMC,13);
predictorStd = predictorMean;
for kk = 1:nMC
    ix = datasample(predictIx,round(p.*numel(predictIx)),'Replace',false);
predictorMean(kk,:) = mean([LiDAR.A.RF(ix),LiDAR.C.C(ix),LiDAR.C.aspect(ix),LiDAR.C.slope(ix),...
    LiDAR.C.gradN(ix),LiDAR.C.gradE(ix),LiDAR.C.aspectN(ix),LiDAR.C.aspectE(ix),...
    LiDAR.C.northness(ix),LiDAR.C.eastness(ix),LiDAR.D.D(ix),LiDAR.E.E(ix),LiDAR.F.F(ix)],1);
predictorStd(kk,:) = std([LiDAR.A.RF(ix),LiDAR.C.C(ix),LiDAR.C.aspect(ix),LiDAR.C.slope(ix),...
    LiDAR.C.gradN(ix),LiDAR.C.gradE(ix),LiDAR.C.aspectN(ix),LiDAR.C.aspectE(ix),...
    LiDAR.C.northness(ix),LiDAR.C.eastness(ix),LiDAR.D.D(ix),LiDAR.E.E(ix),LiDAR.F.F(ix)],[],1);
end
predictorMean = mean(predictorMean);
predictorStd = mean(predictorStd);

nEnsemble = 5;
snowDepthError = zeros(numel(predictIx),nEnsemble);
% Machine Learning Infill Data Set
MLpred = [LiDAR.A.RF(predictIx),LiDAR.C.C(predictIx),LiDAR.C.aspect(predictIx),LiDAR.C.slope(predictIx),...
    LiDAR.C.gradN(predictIx),LiDAR.C.gradE(predictIx),LiDAR.C.aspectN(predictIx),LiDAR.C.aspectE(predictIx),...
    LiDAR.C.northness(predictIx),LiDAR.C.eastness(predictIx),LiDAR.D.D(predictIx),LiDAR.E.E(predictIx),LiDAR.F.F(predictIx)];
MLpred = (MLpred - predictorMean)./predictorStd;

rng(0) % For Reproducibility
probeIx = 1:numel(depthError);
for kk = 1:nEnsemble
    tic
    display(num2str(kk))
ix = datasample(probeIx,round(.1.*numel(probeIx)),'Replace',false);
% Machine Learning Training Data Set
ML = [depthError(ix),LiDAR.A.RF(ixOG(ix)),LiDAR.C.C(ixOG(ix)),LiDAR.C.aspect(ixOG(ix)),LiDAR.C.slope(ixOG(ix)),...
    LiDAR.C.gradN(ixOG(ix)),LiDAR.C.gradE(ixOG(ix)),LiDAR.C.aspectN(ixOG(ix)),LiDAR.C.aspectE(ixOG(ix)),...
    LiDAR.C.northness(ixOG(ix)),LiDAR.C.eastness(ixOG(ix)),LiDAR.D.D(ixOG(ix)),LiDAR.E.E(ixOG(ix)),LiDAR.F.F(ixOG(ix))];
ML (:,2:end) = (ML(:,2:end) - predictorMean)./predictorStd;
[rows, ~] = find(isnan(ML));
rmvIx = unique(rows);
ML(rmvIx,:) = [];

    [trainedModel, validationRMSE] = trainMCSsnowdepthErrorMLR(ML);

snowDepthError(:,kk) = trainedModel.predictFcn(MLpred);
toc

end

avgError = mean(snowDepthError,2);
distributedError = nan(size(LiDAR.A.A));
distributedError(predictIx) = avgError;

%% MAke Pretty Figure
figure();
    hI = imagesc(Coords.X./1000,Coords.Y./1000,distributedError);
    h = colorbar; ylabel(h,'Snow Depth Error (m)','FontWeight','bold','FontName','serif','FontSize',14)
    daspect([1,1,1]);colormap(cmap);caxis([-2 2])%caxis(round([quantile(distributedError(:),[0.01,0.99])./0.05]).*0.05);%caxis([0,5]);
%     title({'Mores Creek Summit: 04/05/2023','MLR Gap Filled Snow Depth'})
    title({strjoin(["Mores Creek Summit: ",num2str(outfnA(9:10)),"/",num2str(outfnA(11:12)),"/",num2str(outfnA(5:8))]),'MLR Snow Depth Error'})
    xlabel('Easting (km)');ylabel('Northing (km)');
    set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
    set(hI,'AlphaData',~mask)
    set(gca,'color',0.75*[1 1 1]);
    %% Write Snow Depth Error Geotiff
    % Need to Somehow Automate EPSG Code Lookup..
    epsgCode = 32611;
    geotiffwrite([writeDir,'\',outfnA(1:23),'_error.tif'],distributedError,Coords.R,'CoordRefSysCode',epsgCode);