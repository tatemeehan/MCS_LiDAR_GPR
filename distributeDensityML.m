%% snowMachine - Supervised Machine Learning to Distribute Snow Density 
clear; close all; clc;
cmap = csvread('RdYlBu.csv');
cmap = flipud(cmap);
%% Options
isWriteLiDARmat = 0;
isWriteGeoTiff = 1;
%% Load Ancillary Data Files from LiDAR_GPR.mat
dataDir = '/bsushare/hpmarshall-shared/LiDAR-GPR/20240315/';
writeDir = dataDir;
% Coordinate Data
load([dataDir,'20240315_MCS-Coords5m.mat'])
% LiDAR Data
load([dataDir,'20240315_MCS-LiDAR5m.mat'])
% LiDAR - GPR data
load([dataDir,'20240315_MCS-LiDAR-GPR_5m.mat'])
% kdTree
load([dataDir,'20240315_MCS-kdtree-5m.mat'])
%% Extract ML Distributors & Predictors

% LiDAR Distributor Data
% Apply the Model to Distribute Density within the Lidar Box
% Spatial Model of Density
Coords.ix = [1:numel(Coords.Xi)]';
Coords.ixMask = Coords.ix(~LiDAR.mask);
% Apply NaN Mask to Boundary Buffer
% MV = [ones(numel(Coords.ixMask),1),LiDAR.A.RF(Coords.ixMask),LiDAR.A.RFnorthness(Coords.ixMask),LiDAR.A.RFeastness(Coords.ixMask),LiDAR.A.RFslope(Coords.ixMask),LiDAR.A.RFgradN(Coords.ixMask),LiDAR.A.RFgradE(Coords.ixMask),LiDAR.A.RFaspect(Coords.ixMask),LiDAR.A.RFaspectN(Coords.ixMask),LiDAR.A.RFaspectE(Coords.ixMask),...
%     LiDAR.B.RF(Coords.ixMask),LiDAR.B.RFnorthness(Coords.ixMask),LiDAR.B.RFeastness(Coords.ixMask),LiDAR.B.RFslope(Coords.ixMask),LiDAR.B.RFgradN(Coords.ixMask),LiDAR.B.RFgradE(Coords.ixMask),LiDAR.B.RFaspect(Coords.ixMask),LiDAR.B.RFaspectN(Coords.ixMask),LiDAR.B.RFaspectE(Coords.ixMask),...
%     LiDAR.C.C(Coords.ixMask),LiDAR.C.northness(Coords.ixMask),LiDAR.C.eastness(Coords.ixMask),LiDAR.C.slope(Coords.ixMask),LiDAR.C.gradN(Coords.ixMask),LiDAR.C.gradE(Coords.ixMask),LiDAR.C.aspect(Coords.ixMask),LiDAR.C.aspectN(Coords.ixMask),LiDAR.C.aspectE(Coords.ixMask),...
%     LiDAR.D.D(Coords.ixMask),LiDAR.E.E(Coords.ixMask),LiDAR.F.F(Coords.ixMask)];

MV = [LiDAR.A.A(Coords.ixMask),LiDAR.A.northness(Coords.ixMask),LiDAR.A.eastness(Coords.ixMask),LiDAR.A.slope(Coords.ixMask),LiDAR.A.gradN(Coords.ixMask),LiDAR.A.gradE(Coords.ixMask),LiDAR.A.aspectN(Coords.ixMask),LiDAR.A.aspectE(Coords.ixMask),...
    LiDAR.B.B(Coords.ixMask),LiDAR.B.northness(Coords.ixMask),LiDAR.B.eastness(Coords.ixMask),LiDAR.B.slope(Coords.ixMask),LiDAR.B.gradN(Coords.ixMask),LiDAR.B.gradE(Coords.ixMask),LiDAR.B.aspectN(Coords.ixMask),LiDAR.B.aspectE(Coords.ixMask),...
    LiDAR.C.C(Coords.ixMask),LiDAR.C.northness(Coords.ixMask),LiDAR.C.eastness(Coords.ixMask),LiDAR.C.slope(Coords.ixMask),LiDAR.C.gradN(Coords.ixMask),LiDAR.C.gradE(Coords.ixMask),LiDAR.C.aspectN(Coords.ixMask),LiDAR.C.aspectE(Coords.ixMask),...
    LiDAR.D.D(Coords.ixMask),LiDAR.E.E(Coords.ixMask),LiDAR.F.F(Coords.ixMask)];

% Standardize Data
ixOG = find(~LiDAR.mask);
p = .001;
nMC = 1000;
predictorMean = zeros(nMC,size(MV,2));
predictorStd = predictorMean;
% Monte Carlo Simulation for Huge Data
for kk = 1:nMC
    ix = datasample(1:numel(ixOG),round(p.*numel(ixOG)),'Replace',false);
predictorMean(kk,:) = nanmean([MV(ix,:)],1);
predictorStd(kk,:) = nanstd([MV(ix,:)],[],1);
end
predictorMean = mean(predictorMean);
predictorStd = mean(predictorStd);
MV = (MV - predictorMean)./predictorStd;

% % Normalize the Data to the IQR
% for kk = 2:size(MV,2)
%     if kk == 29 % was 17
%     else
%         MV(:,kk) = (MV(:,kk)-median(MV(:,kk)))./(iqr(MV(:,kk))+eps);
%     end
% end
ML.distributors = MV;
clear MV

% Random Forest GapFilled Predictor Data
% MV = [ones(length(lidarGPR.Density(:)),1),LiDAR.A.RF(kd.ix),LiDAR.A.RFnorthness(kd.ix),LiDAR.A.RFeastness(kd.ix),LiDAR.A.RFslope(kd.ix),LiDAR.A.RFgradN(kd.ix),LiDAR.A.RFgradE(kd.ix),LiDAR.A.RFaspect(kd.ix),LiDAR.A.RFaspectN(kd.ix),LiDAR.A.RFaspectE(kd.ix),...
%     LiDAR.B.RF(kd.ix),LiDAR.B.RFnorthness(kd.ix),LiDAR.B.RFeastness(kd.ix),LiDAR.B.RFslope(kd.ix),LiDAR.B.RFgradN(kd.ix),LiDAR.B.RFgradE(kd.ix),LiDAR.B.RFaspect(kd.ix),LiDAR.B.RFaspectN(kd.ix),LiDAR.B.RFaspectE(kd.ix),...
%     LiDAR.C.C(kd.ix),LiDAR.C.northness(kd.ix),LiDAR.C.eastness(kd.ix),LiDAR.C.slope(kd.ix),LiDAR.C.gradN(kd.ix),LiDAR.C.gradE(kd.ix),LiDAR.C.aspect(kd.ix),LiDAR.C.aspectN(kd.ix),LiDAR.C.aspectE(kd.ix),...
%     LiDAR.D.D(kd.ix),LiDAR.E.E(kd.ix),LiDAR.F.F(kd.ix)];
MV = [LiDAR.A.A(kd.ix),LiDAR.A.northness(kd.ix),LiDAR.A.eastness(kd.ix),LiDAR.A.slope(kd.ix),LiDAR.A.gradN(kd.ix),LiDAR.A.gradE(kd.ix),LiDAR.A.aspectN(kd.ix),LiDAR.A.aspectE(kd.ix),...
    LiDAR.B.B(kd.ix),LiDAR.B.northness(kd.ix),LiDAR.B.eastness(kd.ix),LiDAR.B.slope(kd.ix),LiDAR.B.gradN(kd.ix),LiDAR.B.gradE(kd.ix),LiDAR.B.aspectN(kd.ix),LiDAR.B.aspectE(kd.ix),...
    LiDAR.C.C(kd.ix),LiDAR.C.northness(kd.ix),LiDAR.C.eastness(kd.ix),LiDAR.C.slope(kd.ix),LiDAR.C.gradN(kd.ix),LiDAR.C.gradE(kd.ix),LiDAR.C.aspectN(kd.ix),LiDAR.C.aspectE(kd.ix),...
    LiDAR.D.D(kd.ix),LiDAR.E.E(kd.ix),LiDAR.F.F(kd.ix)];

% Standardize Predictors
MV = (MV - predictorMean)./predictorStd;

% % Normalize Predictors 
% % Normalize the Data to the IQR
% for kk = 2:size(MV,2)
%     if kk == 29
%     else
%         MV(:,kk) = (MV(:,kk)-median(MV(:,kk)))./(iqr(MV(:,kk))+eps);
%     end
% end
MV = [lidarGPR.Density,MV];
% Add to Machine Learning Structure
ML.predictors = MV;
ML.paramNames = {'\rho_{obs}', 'constant',...
    'Hs','northnessHs','eastnessHS','slopeHs','dyHs','dxHs','aspctHs','aspctNHs','aspctEHs',...
    'Zs','northnessZs','eastnessZs','slopeZs','dyZs','dxZs','aspctZs','aspctNZs','aspctEZs',...
    'Zg','northnessZg','eastnessZg','slopeZg','dyZg','dxZg','aspctZg','aspctNZg','aspctEZg',...
    'Hveg','Sveg','NBR'};
clear MV


% Extract pertainent LiDAR data
A = LiDAR.A.A;
% HillShade
Caspect = LiDAR.C.aspect;
Cnorthness = LiDAR.C.northness;
Ceastness = LiDAR.C.eastness;
% SnowDepth
% depth = LiDAR.A.RF;
depth = LiDAR.A.A;
% Mask
alph=~LiDAR.mask;
% out of memory.
clear LiDAR
%% Machine Learning Regression Ensemble Modeling
modelChoice = [1,4];
for ii = modelChoice %1:3
% Algortihm
if ii == 1
    isMLR = 1;
    isRF = 0;
    isANN = 0;
    isGPR = 0;
end
if ii == 2
    isMLR = 0;
    isRF = 1;
    isANN = 0;
    isGPR = 0;
end
if ii == 3
    isMLR = 0;
    isRF = 0;
    isANN = 1;
    isGPR = 0;
end
if ii == 4
    isMLR = 0;
    isRF = 0;
    isANN = 0;
    isGPR = 1;
end
% if ii == 4
%     nEnsemble = 5;
% else
%     nEnsemble = 5;
% end
p = 0.75;
snowDensityDist = zeros(numel(Coords.ixMask),nEnsemble);
rng(0) % For Reproducibility
trainIx = datasample(1:numel(kd.ix),round(0.95*numel(kd.ix)),'Replace',false);
valIx = find(~ismember(1:numel(kd.ix),trainIx));
for kk = 1:nEnsemble
    tic
    display(num2str(kk))
    ix = datasample(trainIx,round(p*numel(trainIx)),'Replace',false);

% ix = datasample(1:numel(kd.ix),round(p*numel(kd.ix)),'Replace',false);
% Machine Learning Training Data Set
trainML = ML.predictors(ix,:);
[rows, ~] = find(isnan(trainML));
rmvIx = unique(rows);
trainML(rmvIx,:) = [];
if isANN
    [trainedModel, validationRMSE] = trainMCSsnowdensityANN(trainML);
elseif isRF
    [trainedModel, validationRMSE] = trainMCSsnowdensityRF(trainML);
    % [trainedModel, validationRMSE] = trainMCSsnowdensityRF1025(trainML);
elseif isMLR
    % [trainedModel, validationRMSE] = trainMCSsnowdensityMLR(trainML);
    [trainedModel, validationRMSE] = trainMCSsnowdensityMLRiPCA_5m(trainML);
elseif isGPR
    [trainedModel, validationRMSE] = trainMCSsnowdensityGPR_5m(trainML);
end
% Predict in Chunks
nChunks = 1000;
chunks = floor(numel(Coords.ixMask)./nChunks)+1;
for cc = 1:nChunks
%     display(num2str(cc))
    if cc == 1
        predIx = 1:chunks;
    elseif cc == nChunks
        predIx = max(predIx)+1:numel(Coords.ixMask);
    else
        predIx = max(predIx)+1:max(predIx)+chunks;
    end
    snowDensityDist(predIx,kk) = trainedModel.predictFcn(ML.distributors(predIx,:));
end
toc

end

snowDensityDistMean = mean(snowDensityDist,2);
snowDensityDistStd = std(snowDensityDist,[],2);

% Store Distributed Density in ML Structure 
if isMLR
    ML.MLRi.Density = nan(size(A));
    ML.MLRi.Density(Coords.ixMask) = snowDensityDistMean;
    ML.MLRi.stdDensity = nan(size(A));
    ML.MLRi.stdDensity(Coords.ixMask) = snowDensityDistStd;
elseif isRF
%     ML.RF.Density = nan(size(A));
%     ML.RF.Density(Coords.ixMask) = snowDensityDistMean;
%     ML.RF.stdDensity = nan(size(A));
%     ML.RF.stdDensity(Coords.ixMask) = snowDensityDistStd;
    ML.RF1025.Density = nan(size(A));
    ML.RF1025.Density(Coords.ixMask) = snowDensityDistMean;
    ML.RF1025.stdDensity = nan(size(A));
    ML.RF1025.stdDensity(Coords.ixMask) = snowDensityDistStd;
elseif isANN
    ML.ANN.Density = nan(size(A));
    ML.ANN.Density(Coords.ixMask) = snowDensityDistMean;
    ML.ANN.stdDensity = nan(size(A));
    ML.ANN.stdDensity(Coords.ixMask) = snowDensityDistStd;
elseif isGPR
    ML.GPR.Density = nan(size(A));
    ML.GPR.Density(Coords.ixMask) = snowDensityDistMean;
    ML.GPR.stdDensity = nan(size(A));
    ML.GPR.stdDensity(Coords.ixMask) = snowDensityDistStd;
end
%% Create Snow Denisity Figures
if isMLR
    FigMLR = figure('units','normalized','outerposition',[0 0 1 1]);
    % Hillshade
        imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(Caspect+45)+sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
%     imagesc((Coords.X)./1000,(Coords.Y)./1000,(-cosd(Caspect+45)-sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
    % imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
    freezeColors; hold on;
    % MLR Density
    hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLRi.Density,'AlphaData',0.625.*alph);
    daspect([1,1,1]);set(gca,'YDir','normal');
    colormap(cmap);
    % caxis([quantile(ML.MLRi.Density(Coords.ixMask),[0.005,0.995])]);
    clim([300 400])
    cb = colorbar;cb.FontSize = 14;
    cb.Label.String = 'Snow Density (kg/m^3)';
    title(['Mores Creek Summit: 03/15/2024']);%, dataDir(end-5:end)])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
    xtickformat('%.1f')
    ax.YAxis.Exponent = 0;
    ytickformat('%.1f')
    xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
    % Save Figure
    saveas(FigMLR, [writeDir,'\MLRsnowDensity_5m.png'],'png');
elseif isRF
    FigRF = figure('units','normalized','outerposition',[0 0 1 1]);
    % Hillshade
        imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(Caspect+45)+sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
%     imagesc((Coords.X)./1000,(Coords.Y)./1000,(-cosd(Caspect+45)-sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
    % imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
    freezeColors; hold on;
    % RF Density
    hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.RF1025.Density,'AlphaData',0.625.*alph);
    daspect([1,1,1]);set(gca,'YDir','normal');
    colormap(cmap);
    caxis([quantile(ML.MLR.Density(Coords.ixMask),[0.005,0.995])]);
    cb = colorbar;cb.FontSize = 14;
    cb.Label.String = 'RF Snow Density (kg/m^3)';
    title(['Mores Creek Summit Average Snow Density: ', dataDir(end-5:end)])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
    xtickformat('%.1f')
    ax.YAxis.Exponent = 0;
    ytickformat('%.1f')
    xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
    % Save Figure
    saveas(FigRF, [writeDir,'\RF1025snowDensity.png'],'png');
    elseif isANN
    FigANN = figure('units','normalized','outerposition',[0 0 1 1]);
    % Hillshade
        imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(Caspect+45)+sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
%     imagesc((Coords.X)./1000,(Coords.Y)./1000,(-cosd(Caspect+45)-sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
    % imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
    freezeColors; hold on;
    % RF Density
    hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN.Density,'AlphaData',0.625.*alph);
    daspect([1,1,1]);set(gca,'YDir','normal');
    colormap(cmap);
    caxis([quantile(ML.MLR.Density(Coords.ixMask),[0.005,0.995])]);
    cb = colorbar;cb.FontSize = 14;
    cb.Label.String = 'ANN Snow Density (kg/m^3)';
    title(['Mores Creek Summit Average Snow Density: ', dataDir(end-5:end)])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
    xtickformat('%.1f')
    ax.YAxis.Exponent = 0;
    ytickformat('%.1f')
    xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
    % Save Figure
    saveas(FigANN, [writeDir,'\ANNsnowDensityNE.png'],'png');
    elseif isGPR
    FigGPR = figure('units','normalized','outerposition',[0 0 1 1]);
    % Hillshade
        imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(Caspect+45)+sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
%     imagesc((Coords.X)./1000,(Coords.Y)./1000,(-cosd(Caspect+45)-sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
    % imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
    freezeColors; hold on;
    % RF Density
    hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.GPR.Density,'AlphaData',0.625.*alph);
    daspect([1,1,1]);set(gca,'YDir','normal');
    colormap(bone);
    % caxis([quantile(ML.GPR.Density(Coords.ixMask),[0.005,0.995])]);
    clim([320 360])
    cb = colorbar;cb.FontSize = 14;
    cb.Label.String = 'Snow Density (kg/m^3)';
    title(['Mores Creek Summit: 03/15/2024'])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
%     xtickformat('%.1f')
    ax.YAxis.Exponent = 0;
%     ytickformat('%.1f')
    xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
    % Save Figure
    saveas(FigGPR, [writeDir,'\GPRsnowDensity_5m.png'],'png');
end
end

% Write Processed Data!
if isWriteLiDARmat
    save([writeDir,'\',outfnA(1:13),'LiDAR.mat'],'LiDAR','-v7.3')
    save([writeDir,'\',outfnA(1:13),'Coords.mat'],'Coords','-v7.3')
end
if isWriteGeoTiff
    epsgCode = Coords.epsgCode;
    % % Need to Somehow Automate EPSG Code Lookup..
    % auxDir = '/bsuhome/tatemeehan/git-repo/auxData/';
    % tmpTif = [auxDir,'MCS_REFDEM_WGS84.tif'];
    % info = geotiffinfo(tmpTif);
    % epsgCode = info.GeoTIFFCodes.PCS;
    % epsgCode = 32611; % UTM Zone 11
    modelIx = find(ismember(1:4,modelChoice));
    for kk = 1:numel(modelChoice)
        if modelIx(kk) == 1
            geotiffwrite([writeDir,'20240315_MCS','_MLRdensity_5m.tif'],ML.MLRi.Density,Coords.R,'CoordRefSysCode',epsgCode);
        elseif modelIx(kk) == 2
            geotiffwrite([writeDir,'20240315_MCS','_RFdensity_5m.tif'],ML.RF.Density,Coords.R,'CoordRefSysCode',epsgCode);
        elseif modelIx(kk) == 3
            geotiffwrite([writeDir,'20240315_MCS','_ANNdensity_5m.tif'],ML.ANN.Density,Coords.R,'CoordRefSysCode',epsgCode);
        elseif modelIx(kk) == 4
            geotiffwrite([writeDir,'20240315_MCS','_GPRdensity_5m.tif'],ML.GPR.Density,Coords.R,'CoordRefSysCode',epsgCode);

        end
    end
end
keyboard
%% Multivariate Linear Regression
% Cross Validation at p percent
p = .1;
nMC = 250;
betaTrain = zeros(size(ML.MV,2)-1,nMC);
% betaTrain = zeros(length(coms{optOut(optIx,1)}(optOut(optIx,2),:)),nMC);

xvalRMSE = zeros(nMC,1); xvalMAE = xvalRMSE; xvalME = xvalRMSE;xvalR = xvalRMSE;
mcIx = int64(1:numel(kd.ix));
modIx = kd.ix;
for kk = 1:nMC
xvalIx = mcIx;
trainIx = int64(datasample(mcIx,round(p.*numel(modIx)),'replace',false));
xvalIx(trainIx) = [];
% Create the Training Model
betaTrain(:,kk) = mvregress(ML.MV(trainIx,2:end),ML.MV(trainIx,1));
% betaTrain(:,kk) = mvregress(ML.MV(trainIx,coms{optOut(optIx,1)}(optOut(optIx,2),:)+1),ML.MV(trainIx,1));

% Evaluate the Model
xvalDensity = ML.MV(xvalIx,2:end)*betaTrain(:,kk);
% xvalDensity = ML.MV(xvalIx,coms{optOut(optIx,1)}(optOut(optIx,2),:)+1)*betaTrain(:,kk);

% Predition Error
xvalError = ML.MV(xvalIx,1) - xvalDensity;
xvalRMSE(kk) = sqrt(mean(xvalError.^2));
xvalMAE(kk) = mean(abs(xvalError));
xvalME(kk) = mean(xvalError);
xvalR(kk) = corr(ML.MV(xvalIx,1), xvalDensity);
end
% Create the Full Model
betas = mean(betaTrain,2);
% Uncertainty Analysis
stdBeta = std(betaTrain,[],2);
ML.MLR.Density = ML.MV(:,2:end)*betas; ML.MLR.Density = inpaint_nans(ML.MLR.Density,5);
% ML.Density = ML.MV(:,coms{optOut(optIx,1)}(optOut(optIx,2),:)+1)*betas; ML.Density = inpaint_nans(ML.Density,5);


% Apply the Model to Distribute Density within the Lidar Box
% Spatial Model of Density
Coords.ix = [1:numel(Coords.Xi)]';
Coords.ixMask = Coords.ix(~LiDAR.mask);
% MV = [ones(numel(Coords.ix),1),LiDAR.A.A(Coords.ix),LiDAR.A.northness(Coords.ix),LiDAR.A.slope(Coords.ix),LiDAR.A.gradN(Coords.ix),LiDAR.A.gradE(Coords.ix),...
%     LiDAR.B.B(Coords.ix),LiDAR.B.northness(Coords.ix),LiDAR.B.slope(Coords.ix),LiDAR.B.gradN(Coords.ix),LiDAR.B.gradE(Coords.ix),...
%     LiDAR.C.C(Coords.ix),LiDAR.C.northness(Coords.ix),LiDAR.C.slope(Coords.ix),LiDAR.C.gradN(Coords.ix),LiDAR.C.gradE(Coords.ix),...
%     LiDAR.D.D(Coords.ix),LiDAR.E.E(Coords.ix),LiDAR.F.F(Coords.ix)];
% Apply NaN Mask to Boundary Buffer
MV = [ones(numel(Coords.ixMask),1),LiDAR.A.RF(Coords.ixMask),LiDAR.A.RFnorthness(Coords.ixMask),LiDAR.A.RFeastness(Coords.ixMask),LiDAR.A.RFslope(Coords.ixMask),LiDAR.A.RFgradN(Coords.ixMask),LiDAR.A.RFgradE(Coords.ixMask),LiDAR.A.RFaspect(Coords.ixMask),LiDAR.A.RFaspectN(Coords.ixMask),LiDAR.A.RFaspectE(Coords.ixMask),...
    LiDAR.B.RF(Coords.ixMask),LiDAR.B.RFnorthness(Coords.ixMask),LiDAR.B.RFeastness(Coords.ixMask),LiDAR.B.RFslope(Coords.ixMask),LiDAR.B.RFgradN(Coords.ixMask),LiDAR.B.RFgradE(Coords.ixMask),LiDAR.B.RFaspect(Coords.ixMask),LiDAR.B.RFaspectN(Coords.ixMask),LiDAR.B.RFaspectE(Coords.ixMask),...
    LiDAR.C.C(Coords.ixMask),LiDAR.C.northness(Coords.ixMask),LiDAR.C.eastness(Coords.ixMask),LiDAR.C.slope(Coords.ixMask),LiDAR.C.gradN(Coords.ixMask),LiDAR.C.gradE(Coords.ixMask),LiDAR.C.aspect(Coords.ixMask),LiDAR.C.aspectN(Coords.ixMask),LiDAR.C.aspectE(Coords.ixMask),...
    LiDAR.D.D(Coords.ixMask),LiDAR.E.E(Coords.ixMask),LiDAR.F.F(Coords.ixMask)];
% Normalize the Data to the IQR
for kk = 2:size(MV,2)
if kk == 29 % was 17
    else
        MV(:,kk) = (MV(:,kk)-median(MV(:,kk)))./(iqr(MV(:,kk))+eps);
    end
end
ML.distributors = MV;
tmpSpatialDensity = MV*betas;
ML.MLR.spatialDensity = nan(size(LiDAR.A.A));
ML.MLR.spatialDensity(Coords.ixMask) = tmpSpatialDensity;
ML.MLR.coeff = betas;
ML.MLR.coeffStd = stdBeta;
clear tmpSpatialDensity
% % MLR Interactions
% [ML.MLRi.Mdl, ML.MLRi.RMSE] = trainMLRinteractions(ML.MV);
% ML.MLRi.spatialDensity = zeros(length(MV),1);
% tic
% chunks = 1000;
% for kk = 1:chunks
% ML.MLRi.spatialDensity(kk:chunks:end)= ML.MLRi.Mdl.predictFcn(MV(kk:chunks:end,:));
% end
% toc
% ML.MLRi.spatialDensity = reshape(ML.MLRi.spatialDensity,size(LiDAR.A.A));
%% Feature Importance
isOptimizeMVR = 0;
if isOptimizeMVR
paramNames = {'\rho_{o}','Hs','aspctHs','sxHs','dyHs','dxHs',...
    'Zs','aspctZs','sxZs','dyZs','dxZs','Zg','aspctZg','sxZg','dyZg','dxZg','Hveg','Sveg','NBR'};
if exist('optOut')
else
    disp('This will Take a While..')
    mvregX = ML.MV(:,2:end);
    mvregY = ML.MV(:,1);
% Get All combinations of predictors
coms = cell(1,size(ML.MV(:,2:end),2));
for kk = 1:size(ML.MV(:,2:end),2)
% coms{kk} = combntns(1:size(ML.MV(:,2:end),2),kk);
coms{kk} = nchoosek(1:size(ML.MV(:,2:end),2),kk);
end
optOut = [];
%         % Random Sampling
%             p = 0.5;
%             % Randomly Sample .5 Percent
%             randIx = datasample(1:numel(lidarGPR.Density),round(p.*numel(lidarGPR.Density)));
%             tmpDensity = lidarGPR.Density(randIx);
%             evalIx = find(~ismember(1:numel(lidarGPR.Density),randIx));
%             testDensity = lidarGPR.Density(evalIx);
for kk = 1:size(coms,2)
    for jj = 1:size(coms{kk},1)
        % Random Sampling
        p = 0.5;
        % Randomly Sample 50 Percent
        randIx = datasample(1:numel(lidarGPR.Density),round(p.*numel(lidarGPR.Density)));
        tmpDensity = lidarGPR.Density(randIx);
        evalIx = find(~ismember(1:numel(lidarGPR.Density),randIx));
        testDensity = lidarGPR.Density(evalIx);
        testIx = coms{kk}(jj,:);
        % Test All combinations to maximize correlation of Density
        betas = mvregress(mvregX(randIx,testIx),mvregY(randIx));
        % Cross-validation
        tmpMod = mvregX(evalIx,testIx)*betas;
        tmpR = corrcoef(tmpMod,testDensity); tmpR = tmpR(1,2);
        tmpRMSE = sqrt(mean((tmpMod-testDensity).^2));
        optOut = [optOut;[kk,jj,tmpR,tmpRMSE]];
    end
    kk
end
end

% optOut
nparams = [linspace(0.3,1,size(paramNames,2))]';
% optIx = find(optOut(:,3)>0); % Take only Positive Correlations
% optIx = find(optOut(:,3)>0.5); % Take Correlations > 0.5
optIx = find(optOut(:,3)>0.5&optOut(:,4)<30); % R > 0.5 & RMSE < 30
% optIx = find(optOut(:,3)>0.7&optOut(:,4)<30); % R > 0.7 & RMSE < 30

% Detmine Optimal Parameters by statistics
% This is the way
numparams = zeros(numel(optIx),size(coms,2));
whichparams = numparams;
for kk = 1:length(optIx)
    [comIx] = optOut(optIx(kk),1:2);
    numparams(kk,:) = [ismember(1:size(coms,2),comIx(1))];
    whichparams(kk,:) = [ismember(1:size(coms,2),coms{comIx(1)}(comIx(2),:))];
end
nmods = cellfun(@numel,coms);
sumnumparams = sum(numparams);
averagenumparams = sumnumparams./nmods;
sumwhichparams = sum(whichparams);
averagewhichparams = sumwhichparams./length(optIx);
% Find Best model for each number of Params
numk = 100;
% Threshold
qthresh = quantile(optOut(:,3),0.99);
whichoptparams = [];
nummods = 0;
for kk = 1:size(coms,2)
    ixparams{kk} = find(optOut(:,1)==kk);
objparam{kk} = optOut((ixparams{kk}),3)-((optOut((ixparams{kk}),4)-min(optOut((ixparams{kk}),4)))./range(optOut((ixparams{kk}),4)));
% maximum correlation
bestmodix = find(objparam{kk}>qthresh);
if ~isempty(bestmodix)
bestmod{kk} = optOut(ixparams{kk}(bestmodix),(1:2));
    [comIx] = bestmod{kk}(:,1:2);
    numoptparams(kk) = kk;
tmpoptparams = zeros(length(bestmodix),size(coms,2));
for oo = 1:length(bestmodix)
    tmpoptparams(oo,:) = [ismember(1:size(coms,2),coms{comIx(1)}(comIx(oo,2),:))];
end
whichoptparams = [whichoptparams;tmpoptparams];
nummods = nummods+length(bestmodix);
end
end
% Stats
sumwhichoptparams = sum(whichoptparams);
averagewhichoptparams = sumwhichoptparams./(nummods);
ML.MLR.relativeImportance = averagewhichoptparams;
ML.MLR.paramNames = paramNames;

% Plots
% The Objective Function
figure();binscatter(optOut(:,3),optOut(:,4),250)
% All models > 0.7
figure();bar(averagewhichparams,'facecolor','k');
title('Feature Importance in Models with R > 0.5 and RMSE < 30 (kg/m^3)');
ylabel(['Percentage of Models (N = ',num2str(numel(optIx)),')']);xlabel('Feature');
set(gca,'xtick',[1:size(coms,2)],'xticklabel',paramNames,'xticklabelrotation',45,'fontname','serif')
xlim([1,size(coms,2)]); grid on;
% Optimal Models
figure();bar([averagewhichoptparams],'facecolor',[0.5,0.5,0.5],'edgecolor','k');
ylabel(['Relative Importance'],'fontweight','bold');xlabel('Predictor','fontweight','bold');
% ylabel(['Relative Feature Importance (N = ',num2str(nummods),')'],'fontweight','bold');xlabel('Feature','fontweight','bold');
% set(gca,'xtick',[1:size(coms,2)],'xticklabel',paramNames,'xticklabelrotation',45,'fontname','serif')
paramNames3 = {' ','\rho_{o}','Hs','aspctHs','sxHs','dyHs','dxHs',...
    'Zs','aspctZs','sxZs','dyZs','dxZs','Zg','aspctZg','sxZg','dyZg','dxZg','Hveg','Sveg','NBR',' '};
set(gca,'xtick',[0:size(coms,2)+1],'xticklabel',paramNames3,'xticklabelrotation',45,'fontname','serif')
xlim([0,size(coms,2)+1]); 
grid minor;grid on;
% Number of Parameters Versus Which Parameters
figure();
for kk = 1:size(coms,2)
    tmpIx = find(whichoptparams(kk,:));
    hold on;
    plot(tmpIx,kk.*ones(kk,1),'ok','markersize',10,'markerfacecolor','k')
end
ylabel('Number of Features in Optimal Model'); xlabel('Feature')
title('Feature Selection of Optimal Models')
set(gca,'ytick',[1:size(coms,2)],'xtick',[1:size(coms,2)],'xticklabel',paramNames,'xticklabelrotation',45,'fontname','serif')
grid on

end
%% Artificial Neural Network
Nh = round(size(ML.MV,2)./2);
for kk = 1:10
% ML.ANN.net = feedforwardnet(Nh,'trainscg');
ML.ANN.nets{kk}.net = feedforwardnet(Nh,'trainscg');
% ML.ANN.net.trainParam.max_fail = 5;
ML.ANN.nets{kk}.net.trainParam.max_fail = 5;
% ML.ANN.net = train(ML.ANN.net,ML.MV(:,3:end)',ML.MV(:,1)','useGPU','yes','showResources','yes');
ML.ANN.nets{kk}.net = train(ML.ANN.nets{kk}.net,ML.MV(:,3:end)',ML.MV(:,1)','useGPU','yes','showResources','yes');
% ML.ANN.Density = ML.ANN.net(ML.MV(:,3:end)');
% ML.ANN.spatialDensity = ML.ANN.net(MV(:,2:end)');
% ML.ANN.spatialDensity = reshape(ML.ANN.spatialDensity,size(LiDAR.A.A));
ML.ANN.Density{kk} = ML.ANN.nets{kk}.net(ML.MV(:,3:end)');
ML.ANN.spatialDensity{kk} = ML.ANN.nets{kk}.net(MV(:,2:end)');
ML.ANN.spatialDensity{kk} = reshape(ML.ANN.spatialDensity{kk},size(LiDAR.A.A));
end

annDensity = cat(3,ML.ANN.spatialDensity{:});
ML.ANN.ensemble.avgSpatialDensity = mean(annDensity,3);
ML.ANN.ensemble.stdSpatialDensity = std(annDensity,[],3);
ML.ANN.ensemble.avgDensity = ML.ANN.ensemble.avgSpatialDensity(kd.ix);
ML.ANN.ensemble.stdDensity = ML.ANN.ensemble.stdSpatialDensity(kd.ix);
ML.ANN.ensemble.spatialDensity = ML.ANN.spatialDensity;
ML.ANN.ensemble.Density = ML.ANN.Density;

% [ML.ANN10.Mdl, ML.ANN10.RMSE] = trainANN10(ML.MV);
% ML.ANN10.spatialDensity = ML.ANN10.Mdl.predictFcn(MV);
% ML.ANN10.spatialDensity = reshape(ML.ANN10.spatialDensity,size(LiDAR.A.A));
% 
% [ML.ANN25.Mdl, ML.ANN25.RMSE] = trainANN25(ML.MV);
% ML.ANN25.spatialDensity = ML.ANN25.Mdl.predictFcn(MV);
% ML.ANN25.spatialDensity = reshape(ML.ANN25.spatialDensity,size(LiDAR.A.A));

% Relative Importance 
% https://www.mathworks.com/matlabcentral/answers/299646-how-to-obtain-the-relative-importance-of-each-input-variable-for-a-neural-network
W1= ML.ANN.net.IW{1,1}; % Weights between input and hidden layers
W2= ML.ANN.net.LW{2,1}'; % Weights between hidden and output layers
k= size(W1,2); % number of input varables
g=1; % Select which output you are going to evaluate the relative importance of 
% input variables on (No need to change if you have only one output)
for n=1:k
for i=1:k
    a(:,i)=(abs(W1(:,n))./sum(abs(W1(:,1:i)),2)).*abs(W2(:,g));
end
b=sum(sum(a,2));
I(n,1)=sum(abs(W1(:,n))./sum(abs(W1),2))/b;
end
for i=1:k
    ML.ANN.R(i,1)=(I(i)/sum(I))*100; % Percentage
end
% Ensemble importance
for kk = 1:10
W1= ML.ANN.nets{kk}.net.IW{1,1}; % Weights between input and hidden layers
W2= ML.ANN.nets{kk}.net.LW{2,1}'; % Weights between hidden and output layers
weights1{kk} = W1;
weights2{kk} = W2;
bias{kk} = ML.ANN.nets{kk}.net.b;
k= size(W1,2); % number of input varables
g=1; % Select which output you are going to evaluate the relative importance of 
% input variables on (No need to change if you have only one output)
for n=1:k
for i=1:k
    a(:,i)=(abs(W1(:,n))./sum(abs(W1(:,1:i)),2)).*abs(W2(:,g));
end
b=sum(sum(a,2));
I(n,1)=sum(abs(W1(:,n))./sum(abs(W1),2))/b;
end
for i=1:k
    ML.ANN.nets{kk}.R(i,1)=(I(i)/sum(I))*100; % Percentage
    netsR(i,kk) = (I(i)/sum(I))*100;
end
end
ML.ANN.ensemble.relativeImportance = mean(netsR,2);
ML.ANN.ensemble.relativeImportanceStd = std(netsR,[],2);
paramNames = {'Hs','aspctHs','sxHs','dyHs','dxHs',...
    'Zs','aspctZs','sxZs','dyZs','dxZs','Zg','aspctZg','sxZg','dyZg','dxZg','Hveg','Sveg','NBR'};
ML.ANN.ensemble.paramNames = paramNames;
% % Reproduce Average from Average of Weights and Biases
% for kk = 1:10
%     bias1(:,kk) = bias{kk}{1};
%     bias2(kk) = bias{kk}{2};
%     weight1(:,:,kk) = weights1{kk};
%     weight2(:,kk) = weights2{kk};
% end
% bias1 = mean(bias1,2);
% bias2 = mean(bias2);
% weight1 = mean(weight1,3);
% weight2 = mean(weight2,2)';
% ML.ANN.avgNet = ML.ANN.nets{1}.net;
% ML.ANN.avgNet.b = {bias1;bias2};
% ML.ANN.avgNet.IW = {weight1;[]};
% ML.ANN.avgNet.LW = {[],[];weight2,[]};
% % ML.ANN.avgNet.spatialDensity = ML.ANN.avgNet(MV(:,2:end)');
% test = ML.ANN.avgNet(MV(:,2:end)');
% test = reshape(test,size(LiDAR.A.A));
% Optimal Models
% figure();bar([ML.ANN.R],'facecolor',[0.5,0.5,0.5],'edgecolor','k');
figure();bar([mean(netsR,2)],'facecolor',[0.5,0.5,0.5],'edgecolor','k');
hold on;
% errorbar(1:18,mean(netsR,2),mean(netsR,2)-std(netsR,[],2),mean(netsR,2)+std(netsR,[],2),'k','linewidth',2,'LineStyle','none');   
ylabel(['Relative Importance (%)'],'fontweight','bold');xlabel('Predictor','fontweight','bold');
paramNames3 = {' ','Hs','aspctHs','sxHs','dyHs','dxHs',...
    'Zs','aspctZs','sxZs','dyZs','dxZs','Zg','aspctZg','sxZg','dyZg','dxZg','Hveg','Sveg','NBR',' '};
set(gca,'xtick',[0:numel(ML.ANN.R)],'xticklabel',paramNames3,'xticklabelrotation',45,'fontname','serif')
xlim([0,numel(ML.ANN.R)+1]); 
grid minor;grid on;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
%% Random Forest Regression
rng(1); % For reproducibility
ML.RF.NumTrees = 30;
ML.RF.Mdl = TreeBagger(ML.RF.NumTrees,ML.MV(:,3:end),ML.MV(:,1),Method="regression",Surrogate="on",OOBPredictorImportance="on",PredictorSelection="interaction-curvature",MinLeafSize=8,NumPredictorsToSample=7) ;
ML.RF.Density = predict(ML.RF.Mdl,ML.MV(:,3:end));
ML.RF.spatialDensity = predict(ML.RF.Mdl,MV(:,2:end));
ML.RF.spatialDensity = reshape(ML.RF.spatialDensity,size(LiDAR.A.A));

% [ML.BT, ML.BT.RMSE] = trainBaggedTrees(ML.MV);
% ML.BT.spatialDensity = ML.BT.predictFcn(MV);
% ML.BT.spatialDensity = reshape(ML.BT.spatialDensity,size(LiDAR.A.A));

% Predictor Importance
ML.RF.OOBPermutedPredictorRelativeImportance = (ML.RF.Mdl.OOBPermutedPredictorDeltaError./sum(ML.RF.Mdl.OOBPermutedPredictorDeltaError)).*100;
paramNames = {'Hs','aspctHs','sxHs','dyHs','dxHs',...
    'Zs','aspctZs','sxZs','dyZs','dxZs','Zg','aspctZg','sxZg','dyZg','dxZg','Hveg','Sveg','NBR'};
ML.RF.paramNames = paramNames;
figure();bar([ML.RF.OOBPermutedPredictorRelativeImportance],'facecolor',[0.5,0.5,0.5],'edgecolor','k');
ylabel(['Relative Importance (%)'],'fontweight','bold');xlabel('Predictor','fontweight','bold');
paramNames3 = {' ','Hs','aspctHs','sxHs','dyHs','dxHs',...
    'Zs','aspctZs','sxZs','dyZs','dxZs','Zg','aspctZg','sxZg','dyZg','dxZg','Hveg','Sveg','NBR',' '};
set(gca,'xtick',[0:size(coms,2)+1],'xticklabel',paramNames3,'xticklabelrotation',45,'fontname','serif')
xlim([0,numel(ML.ANN.R)+1]); 
grid minor;grid on;
set(gca,'fontname','serif','fontweight','bold','fontsize',12)

%% Plot Density Figure
alph=~LiDAR.mask;

figure();
% Hillshade
% imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.C.northness);colormap(bone)
imagesc((Coords.X)./1000,(Coords.Y)./1000,(+cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
% imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
freezeColors; hold on;

% MLR
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLR.spatialDensity,'AlphaData',0.625.*alph); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLRi.spatialDensity); 
% ANN
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN.ensemble.avgSpatialDensity,'AlphaData',0.625); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN.spatialDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN10.spatialDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN25.spatialDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,meanAnnDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,stdAnnDensity); 



% RF
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.RF.spatialDensity,'AlphaData',0.625); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.BT.spatialDensity); 

% Add Snow Pits
% hold on; scatter(Xs(allSnowPitIx)./1000,Ys(allSnowPitIx)./1000,35,snowPitDensityFout,'filled','markeredgecolor','k')
% set(hI,'AlphaData',alph)

% NAN Mask
% imagesc(Data_Array,'AlphaData',imAlpha);
% set(gca,'color',0*[1 1 1]);


daspect([1,1,1]);set(gca,'YDir','normal');
colormap(cmap); 
% caxis([300,450]);
caxis([quantile(ML.MLR.spatialDensity(Coords.ixMask),[0.005,0.995])]);
% caxis([quantile(ML.ANN.ensemble.avgSpatialDensity(:),[0.005,0.995])]);
% imagesc((Coords.X)./1000,(Coords.Y)./1000,~alph,'AlphaData',alph)
% set(gca,'color',0.*[0.75,0.75,0.75])
cb = colorbar;cb.FontSize = 14;
% xlabel('Easting (km)'); ylabel('Northing (km)');
cb.Label.String = 'Average Snow Density (kg/m^3)'; %title('MLR Modeled Average Snow Density (kg/m^3)')
% xlabel('Easting (m)'); ylabel('Northing (m)');title('MLR Modeled Average Snow Density (kg/m^3)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
% hold on; plot(Xs(snowPitIx),Ys(snowPitIx),'ok','markerfacecolor','k')
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'position',pos)
% xlim([741.5,746.0])   
% ylim([4321.5,4325.0])
xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
%% SWE Figure
figure();
% Hillshade
% imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.C.northness);colormap(bone)
imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45))./2);colormap(bone)
freezeColors; hold on;
% MLR
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLR.spatialSWE,'AlphaData',0.625);
% ANN
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN.ensemble.spatialSWE,'AlphaData',0.625); 
% RF
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.RF.spatialSWE,'AlphaData',0.625); 
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLR.spatialSWE - ML.ANN.ensemble.spatialSWE ,'AlphaData',0.625); 
daspect([1,1,1]);set(gca,'YDir','normal');
colormap(cmap); 
caxis([200,800]);
%caxis([quantile(spatialDensity(:),[0.025,0.975])]);
cb = colorbar;cb.FontSize = 14;
% xlabel('Easting (km)'); ylabel('Northing (km)');
cb.Label.String = 'Snow Water Equivalent (mm)'; %title('MLR Modeled Average Snow Density (kg/m^3)')
% xlabel('Easting (m)'); ylabel('Northing (m)');title('MLR Modeled Average Snow Density (kg/m^3)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
% hold on; plot(Xs(snowPitIx),Ys(snowPitIx),'ok','markerfacecolor','k')
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'position',pos)
% xlim([741.5,746.0])   
% ylim([4321.5,4325.0])
xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
%% Snow Depth Figure
figure();
% Hillshade
% imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.C.northness);colormap(bone)
imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45))./2);colormap(bone)
freezeColors; hold on;
% LiDAR Snow Depth
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.A.A.*100,'AlphaData',0.625); 
daspect([1,1,1]);set(gca,'YDir','normal');
colormap(cmap); 
caxis([50,250]);
%caxis([quantile(spatialDensity(:),[0.025,0.975])]);
cb = colorbar;cb.FontSize = 14;
% xlabel('Easting (km)'); ylabel('Northing (km)');
cb.Label.String = 'Depth (cm)'; %title('MLR Modeled Average Snow Density (kg/m^3)')
% xlabel('Easting (m)'); ylabel('Northing (m)');title('MLR Modeled Average Snow Density (kg/m^3)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
% hold on; plot(Xs(snowPitIx),Ys(snowPitIx),'ok','markerfacecolor','k')
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'position',pos)
% xlim([741.5,746.0])   
% ylim([4321.5,4325.0])
xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
%% Plot Density Difference Figure
figure();
% Hillshade
% imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.C.northness);colormap(bone)
imagesc((Coords.X)./1000,(Coords.Y)./1000,(cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45))./2);colormap(bone)
freezeColors; hold on;

% MLR
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLR.spatialDensity-ML.ANN.ensemble.avgSpatialDensity,'AlphaData',0.625); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.MLRi.spatialDensity); 
% ANN
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN.ensemble.avgSpatialDensity,'AlphaData',0.625); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN.spatialDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN10.spatialDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.ANN25.spatialDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,meanAnnDensity); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,stdAnnDensity); 



% RF
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.RF.spatialDensity,'AlphaData',0.625); 
% hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,ML.BT.spatialDensity); 

% hold on; scatter(Xs(allSnowPitIx)./1000,Ys(allSnowPitIx)./1000,35,snowPitDensityFout,'filled','markeredgecolor','k')
daspect([1,1,1]);set(gca,'YDir','normal');
colormap(cmap); 
caxis([300,450]);
%caxis([quantile(spatialDensity(:),[0.025,0.975])]);
cb = colorbar;cb.FontSize = 14;
% xlabel('Easting (km)'); ylabel('Northing (km)');
cb.Label.String = 'Average Snow Density (kg/m^3)'; %title('MLR Modeled Average Snow Density (kg/m^3)')
% xlabel('Easting (m)'); ylabel('Northing (m)');title('MLR Modeled Average Snow Density (kg/m^3)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
% hold on; plot(Xs(snowPitIx),Ys(snowPitIx),'ok','markerfacecolor','k')
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'position',pos)
% xlim([741.5,746.0])   
% ylim([4321.5,4325.0])
xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
%% NBR 
figure();
% hI = imagesc(X,Y,spatialDensity); 
hI = imagesc((Coords.X)./1000,(Coords.Y)./1000,LiDAR.F.F); 
% hI = imagesc(X./1000,Y./1000,densityBaggedTrees-spatialDensity); 
% hold on; plot(GPR.Easting./1000,GPR.Northing./1000,'.m')
% hold on; scatter(Xs(allSnowPitIx)./1000,Ys(allSnowPitIx)./1000,35,snowPitDensityFout,'filled','markeredgecolor','k')
daspect([1,1,1]);set(gca,'YDir','normal');
colormap(bone); 
% caxis([300,450]);
%caxis([quantile(spatialDensity(:),[0.025,0.975])]);
cb = colorbar;cb.FontSize = 14;
% xlabel('Easting (km)'); ylabel('Northing (km)');
cb.Label.String = 'Normalized Burn Ratio'; %title('MLR Modeled Average Snow Density (kg/m^3)')
% xlabel('Easting (m)'); ylabel('Northing (m)');title('MLR Modeled Average Snow Density (kg/m^3)')
set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
% hold on; plot(Xs(snowPitIx),Ys(snowPitIx),'ok','markerfacecolor','k')
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'position',pos)
% xlim([741.5,746.0])
% ylim([4321.5,4325.0])
xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);