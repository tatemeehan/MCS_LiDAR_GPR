%% Post-Process LiDAR
clear;close all; clc;
cmap = csvread('./colormaps/RdYlBu.csv');
cmap = flipud(cmap);
%% Processing Steps
% Write Output as GeoTiff
isWriteGeoTiff = 1;
% Write LiDAR as .mat
isWriteLiDARmat = 1;
% Data Directory 
dataDir = '/bsushare/hpmarshall-shared/LiDAR-GPR/20240418/';
% Write Directory
writeDir = dataDir;
%% Load Raster Layers
% A Snow Depth
filename = '20220317_MCS-snowdepth.tif';
outfnA = filename(1:end-4);
fullfilename = fullfile(dataDir,filename);
[A,RA,~,~,lon,lat,utmX,utmY] = readLidarTif(fullfilename);
tmpix = ~isnan(A);
threshold = quantile(A(tmpix),0.999);% MAX snowdepth
% threshold = 5.25; % MAX snowdepth
A(A>threshold) = NaN;
A(A<0) = NaN;
xq = RA.XWorldLimits(1):RA.CellExtentInWorldX:RA.XWorldLimits(2);
yq = RA.YWorldLimits(1):RA.CellExtentInWorldY:RA.YWorldLimits(2);
[Xq,Yq] = meshgrid(xq,yq);
Yq = flipud(Yq);
% Get UTM Coordinates as Vector
X = utmX(1,:);
Y = utmY(:,1);
Xi = utmX(:); 
Yi = utmY(:);
% Create Coordinate Structure
Coords.R = RA;
Coords.lon = lon;
Coords.lat = lat;
Coords.utmX = utmX;
Coords.utmY = utmY;
Coords.X = X;
Coords.Y = Y;
Coords.Xi = Xi;
Coords.Yi = Yi;
clear('X','Y','Xi','Yi','lon','lat','utmX','utmY')

% B Snow Elevation
% filename = '20220317_MCS-snow-T8.tif';
filename = '20220317_MCS-snow.tif';
outfnB = filename(1:end-4);
fullfilename = fullfile(dataDir,filename);
[B,RB,~,~,~,~,~,~] = readLidarTif(fullfilename);
B = mapinterp(B,RB,Xq,Yq);
% Resize Snow On
B = imresize(B,size(A));
nanix = find(B<0 | B > 3000);
B(nanix) = nan; 
clear('RB')

% C Reference DEM
dataDir2 = '/bsuhome/tatemeehan/git-repo/auxData';
% filename = '20220317_MCS-refDEM-T8.tif';
filename = 'MCS_REFDEM_WGS84.tif';
fullfilename = fullfile(dataDir2,filename);
[C,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
nanix = find(C<0);
C(nanix) = nan;
nanMask = double(~isnan(C));
nanMask(nanMask ~= 1) = NaN;

% D Canopy Height
% filename = '20220317_MCS-canopyheight-T8.tif';
% filename = 'hanover_reprocess-canopy.tif';
filename = 'chm_mcs_1m.tif';
fullfilename = fullfile(dataDir2,filename);
[D,RD,~,~,~,~,~,~] = readLidarTif(fullfilename);
D = mapinterp(D,RD,Xq,Yq);
D = imresize(D,size(A));
zedix = find(D<0 & D>=-1e1);
D(zedix) = 0;
nanix = find(D<-1e1 | D > 1e3);
D(nanix) = nan;
clear('RD','zedix')

% E Canopy Proximity
filename = 'proximity_map_mcs_1m.tif';
fullfilename = fullfile(dataDir2,filename);
[E,RE,~,~,~,~,~,~] = readLidarTif(fullfilename);
E = mapinterp(E,RE,Xq,Yq);
E = imresize(E,size(A));
clear('RE')

% F Normalized Burn Ratio
filename =  'dnbr_mcs_utm_11n.tif';
fullfilename = fullfile(dataDir2,filename);
[F,RF,~,~,~,~,~,~] = readLidarTif(fullfilename);
% Trim NBR
[F,RF] = mapcrop(F,RF,RA.XWorldLimits,RA.YWorldLimits);
% Resize/Scale NBR to 0.5m
scale = round(RF.CellExtentInWorldX./RA.CellExtentInWorldX);
[F,RF] = mapresize(F,RF,scale,'cubic');
[F,~] = mapcrop(F,RF,RA.XWorldLimits,RA.YWorldLimits);
F = imresize(F,size(A));
clear('RF','RA','Xq','Yq','xq','yq')

% NaN Mask
border = csvread([dataDir2,'/MCSborder.csv']);
mask = isnan(C)+border;
mask = mask>0;
clear('border')
%% Intelligent Gap Filling
% Inpaint Nans on MCS Snow Depth Predictors
D = inp8nt_nans(D,mask);
% Vegheight Threshold
Vix = D>0.5;
E = inp8nt_nans(E,mask);
F = inp8nt_nans(F,mask);

% Create LiDAR Structure
LiDAR.A.Feature = 'SnowDepth';
LiDAR.B.Feature = 'SnowElevation';
LiDAR.C.Feature = 'ReferenceDEM';
LiDAR.D.Feature = 'CanopyHeight';
LiDAR.E.Feature = 'CanopyProximity';
LiDAR.F.Feature = 'NormalizedBurnRatio';
% Input Rasters
LiDAR.A.A = A;
LiDAR.B.B = B;
LiDAR.C.C = C;
LiDAR.D.D = D;
LiDAR.D.Vix = Vix;
LiDAR.E.E = E;
LiDAR.F.F = F;
LiDAR.mask = mask;
clear('A', "B", 'C', 'D', "E", "F","Vix")

% Median Smoothing
% LiDAR.A.A = medfilt2(LiDAR.A.A,[5,5]);
% LiDAR.B.B = medfilt2(LiDAR.B.B,[5,5]);
% LiDAR.C.C = medfilt2(LiDAR.C.C,[5,5]);
tmpMask = isnan(LiDAR.C.C);
tmpC = medfilt2(LiDAR.C.C,[5,5]);
tmpMask2 = isnan(tmpC);
tmpix = find(abs(tmpMask-tmpMask2)>0);
tmpC(tmpix) = LiDAR.C.C(tmpix);
LiDAR.C.C = tmpC;
clear("tmpMask","tmpMask2","tmpC")


% Create Georeference
latlim = [min(Coords.lat(:)),max(Coords.lat(:))];
lonlim = [min(Coords.lon(:)) max(Coords.lon(:))];
sizeLidar = size(LiDAR.A.A);
georef = georefpostings(latlim,lonlim,sizeLidar,'RowsStartFrom','west','ColumnsStartFrom','north');
clear("latlim","lonlim","sizeLidar",'shadIx')

% % Snow Depth Derivatives
% [LiDAR.A.aspect,LiDAR.A.slope,LiDAR.A.gradN,LiDAR.A.gradE] = gradientm(LiDAR.A.A,georef);
% LiDAR.A.aspect = inpaint_nans(LiDAR.A.aspect,5);
% 
% % Snow-On Elevation Derivatives
% [LiDAR.B.aspect,LiDAR.B.slope,LiDAR.B.gradN,LiDAR.B.gradE] = gradientm(LiDAR.B.B,georef);
% LiDAR.B.aspect = inpaint_nans(LiDAR.B.aspect,5);

% Bare Ground Elevation Derivatives
[LiDAR.C.aspect,LiDAR.C.slope,LiDAR.C.gradN,LiDAR.C.gradE] = gradientm(LiDAR.C.C,georef);
% LiDAR.C.aspect = inpaint_nans(LiDAR.C.aspect,5);
% LiDAR.C.aspect = gapfillMCS(LiDAR.C.aspect);
% LiDAR.C.slope = gapfillMCS(LiDAR.C.slope);
% LiDAR.C.gradN = gapfillMCS(LiDAR.C.gradN);
% LiDAR.C.gradE = gapfillMCS(LiDAR.C.gradE);
LiDAR.C.aspect = inp8nt_nans(LiDAR.C.aspect,mask);
LiDAR.C.slope = inp8nt_nans(LiDAR.C.slope,mask);
LiDAR.C.gradN = inp8nt_nans(LiDAR.C.gradN,mask);
LiDAR.C.gradE = inp8nt_nans(LiDAR.C.gradE,mask);
% tmpCaspect = LiDAR.C.aspect(tmpix);
% tmpCslope = LiDAR.C.slope(tmpix);
% tmpCgradn = LiDAR.C.gradN(tmpix);
% tmpCgrade = LiDAR.C.gradE(tmpix);

    % Predictors
    if exist('smoothFlag')
    else
        smoothFlag = 1;
%         LiDAR.A.aspect = medfilt2(LiDAR.A.aspect,[25,25]);
%         LiDAR.A.slope = medfilt2(LiDAR.A.slope,[25,25]);
%         LiDAR.A.gradN = medfilt2(LiDAR.A.gradN,[25,25]);
%         LiDAR.A.gradE = medfilt2(LiDAR.A.gradE,[25,25]);
%         LiDAR.B.aspect = medfilt2(LiDAR.B.aspect,[25,25]);
%         LiDAR.B.slope = medfilt2(LiDAR.B.slope,[25,25]);
%         LiDAR.B.gradN = medfilt2(LiDAR.B.gradN,[25,25]);
%         LiDAR.B.gradE = medfilt2(LiDAR.B.gradE,[25,25]);
        LiDAR.C.aspect = medfilt2(LiDAR.C.aspect,[5,5]);
%         LiDAR.C.aspect(tmpix) = tmpCaspect;
        LiDAR.C.slope = medfilt2(LiDAR.C.slope,[5,5]);
%         LiDAR.C.slope(tmpix) = tmpCslope;
        LiDAR.C.gradN = medfilt2(LiDAR.C.gradN,[5,5]);
%         LiDAR.C.gradN(tmpix) = tmpCgradn;
        LiDAR.C.gradE = medfilt2(LiDAR.C.gradE,[5,5]);
%         LiDAR.C.gradE(tmpix) = tmpCgrade;
        LiDAR.C.aspect = inp8nt_nans(LiDAR.C.aspect,mask);
        LiDAR.C.slope = inp8nt_nans(LiDAR.C.slope,mask);
        LiDAR.C.gradN = inp8nt_nans(LiDAR.C.gradN,mask);
        LiDAR.C.gradE = inp8nt_nans(LiDAR.C.gradE,mask);
        clear('tmpCaspect',"tmpCslope","tmpCgradn","tmpCgrade",'tmpix')
        % Aspect North Normalized
%         LiDAR.A.aspectN = cosd(LiDAR.A.aspect);
%         LiDAR.B.aspectN = cosd(LiDAR.B.aspect);
        LiDAR.C.aspectN = cosd(LiDAR.C.aspect);
        % Aspect East Normalized
%         LiDAR.A.aspectE = sind(LiDAR.A.aspect);
%         LiDAR.B.aspectE = sind(LiDAR.B.aspect);
        LiDAR.C.aspectE = sind(LiDAR.C.aspect);
        % Northness
%         LiDAR.A.northness = cosd(LiDAR.A.aspect).*sind(LiDAR.A.slope);
%         LiDAR.B.northness = cosd(LiDAR.B.aspect).*sind(LiDAR.B.slope);
        LiDAR.C.northness = cosd(LiDAR.C.aspect).*sind(LiDAR.C.slope);
        % Eastness
%         LiDAR.A.eastness = sind(LiDAR.A.aspect).*sind(LiDAR.A.slope);
%         LiDAR.B.eastness = sind(LiDAR.B.aspect).*sind(LiDAR.B.slope);
        LiDAR.C.eastness = sind(LiDAR.C.aspect).*sind(LiDAR.C.slope);
        % Normalize aspect between 0 and 1;
%         LiDAR.A.northness = abs(LiDAR.A.aspect-180)./180;
%         LiDAR.B.northness = abs(LiDAR.B.aspect-180)./180;
%         LiDAR.C.northness = abs(LiDAR.C.aspect-180)./180;
%         LiDAR.A.northness = cosd(LiDAR.A.aspect+45)+sind(LiDAR.A.aspect+45);
%         LiDAR.B.northness = cosd(LiDAR.B.aspect+45)+sind(LiDAR.B.aspect+45);
%         LiDAR.C.northness = cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45);
    end
%% Machine Learining for Gap Filling Snow Depths
for ii = 1:3
% Algortihm
if ii == 1
    isMLR = 1;
    isRF = 0;
    isANN = 0;
end
if ii == 2
    isMLR = 0;
    isRF = 1;
    isANN = 0;
end
if ii == 3
    isMLR = 0;
    isRF = 0;
    isANN = 1;
end
% Standardize Data
% ix = find(~isnan(LiDAR.C.C));
ixOG = find(~mask);
p = .0001;
nMC = 2500;
predictorMean = zeros(nMC,12);
predictorStd = predictorMean;
for kk = 1:nMC
    ix = datasample(ixOG,round(p.*numel(ixOG)),'Replace',false);
predictorMean(kk,:) = mean([LiDAR.C.C(ix),LiDAR.C.aspect(ix),LiDAR.C.slope(ix),...
    LiDAR.C.gradN(ix),LiDAR.C.gradE(ix),LiDAR.C.aspectN(ix),LiDAR.C.aspectE(ix),...
    LiDAR.C.northness(ix),LiDAR.C.eastness(ix),LiDAR.D.D(ix),LiDAR.E.E(ix),LiDAR.F.F(ix)],1);
predictorStd(kk,:) = std([LiDAR.C.C(ix),LiDAR.C.aspect(ix),LiDAR.C.slope(ix),...
    LiDAR.C.gradN(ix),LiDAR.C.gradE(ix),LiDAR.C.aspectN(ix),LiDAR.C.aspectE(ix),...
    LiDAR.C.northness(ix),LiDAR.C.eastness(ix),LiDAR.D.D(ix),LiDAR.E.E(ix),LiDAR.F.F(ix)],[],1);
end
predictorMean = mean(predictorMean);
predictorStd = mean(predictorStd);

% ixOG = find(~isnan(LiDAR.A.A));
ixOG = find(~isnan(LiDAR.A.A)&~isnan(LiDAR.B.B));
gapIxA = find(isnan(LiDAR.A.A));
gapIxB = find(isnan(LiDAR.B.B));
gapIx = unique([gapIxA(:);gapIxB(:)]);
maskix = find(mask);
gapIx = setdiff(gapIx,maskix);
% gapIxA = find(ismember(gapIxA,gapIx));
% gapIxB = find(ismember(gapIxB,gapIx));
% gapIxA = setdiff(gapIxA,maskix);
% gapIxB = intersect(gapIx,gapIxB);
% wtfIxA = find(ismember(gapIxA,gapIx));
% wtfIxB = find(ismember(gapIxB,gapIx));

nEnsemble = 5;
snowDepthGap = zeros(numel(gapIx),nEnsemble);
% Machine Learning Infill Data Set
MLgap = [LiDAR.C.C(gapIx),LiDAR.C.aspect(gapIx),LiDAR.C.slope(gapIx),...
    LiDAR.C.gradN(gapIx),LiDAR.C.gradE(gapIx),LiDAR.C.aspectN(gapIx),LiDAR.C.aspectE(gapIx),...
    LiDAR.C.northness(gapIx),LiDAR.C.eastness(gapIx),LiDAR.D.D(gapIx),LiDAR.E.E(gapIx),LiDAR.F.F(gapIx)];
MLgap = (MLgap - predictorMean)./predictorStd;

rng(0) % For Reproducibility
for kk = 1:nEnsemble
    tic
    display(num2str(kk))
ix = datasample(ixOG,round(.0005.*numel(ixOG)),'Replace',false);
% Machine Learning Training Data Set
ML = [LiDAR.A.A(ix),LiDAR.C.C(ix),LiDAR.C.aspect(ix),LiDAR.C.slope(ix),...
    LiDAR.C.gradN(ix),LiDAR.C.gradE(ix),LiDAR.C.aspectN(ix),LiDAR.C.aspectE(ix),...
    LiDAR.C.northness(ix),LiDAR.C.eastness(ix),LiDAR.D.D(ix),LiDAR.E.E(ix),LiDAR.F.F(ix)];
ML (:,2:end) = (ML(:,2:end) - predictorMean)./predictorStd;
[rows, ~] = find(isnan(ML));
rmvIx = unique(rows);
ML(rmvIx,:) = [];
if isANN
    [trainedModel, validationRMSE] = trainMCSsnowdepthANN(ML);
elseif isRF
    [trainedModel, validationRMSE] = trainRandomForestMCSsnowDepth(ML);
elseif isMLR
    [trainedModel, validationRMSE] = trainMCSsnowdepthMLR(ML);
end

snowDepthGap(:,kk) = trainedModel.predictFcn(MLgap);
toc

end

% Gap Fill Snow Depth and Snow Surface Elevation.
snowDepthGapMean = mean(snowDepthGap,2);
snowDepthGapStd = std(snowDepthGap,[],2);

tmpA = LiDAR.A.A;
if isMLR
    LiDAR.A.MLR = LiDAR.A.A;
    LiDAR.A.MLR(gapIx) = snowDepthGapMean;
    LiDAR.B.MLR = LiDAR.B.B;
    LiDAR.B.MLR(gapIx) = LiDAR.C.C(gapIx)+snowDepthGapMean;
elseif isRF
    LiDAR.A.RF = LiDAR.A.A;
    LiDAR.A.RF(gapIx) = snowDepthGapMean;
    LiDAR.B.RF = LiDAR.B.B;
    LiDAR.B.RF(gapIx) = LiDAR.C.C(gapIx)+snowDepthGapMean;
elseif isANN
    LiDAR.A.ANN = LiDAR.A.A;
    LiDAR.A.ANN(gapIx) = snowDepthGapMean;
    LiDAR.B.ANN = LiDAR.B.B;
    LiDAR.B.ANN(gapIx) = LiDAR.C.C(gapIx)+snowDepthGapMean;
end

%% Create Snow Depth Figures
alph=isnan(tmpA);
if isMLR
%     FigMLR = figure('Position', get(0, 'Screensize'));
    FigMLR = figure('units','normalized','outerposition',[0 0 1 1]);
elseif isRF
%     FigRF = figure('Position', get(0, 'Screensize'));
    FigRF = figure('units','normalized','outerposition',[0 0 1 1]);
elseif isANN
%     FigANN = figure('Position', get(0, 'Screensize'));
    FigANN = figure('units','normalized','outerposition',[0 0 1 1]);

end
subplot(1,2,1)
hI = imagesc(Coords.X./1000,Coords.Y./1000,tmpA);
h = colorbar; ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif','FontSize',14)
daspect([1,1,1]);colormap(cmap);caxis(round([quantile(tmpA(:),[0.1,0.99])./0.05]).*0.05);%caxis([0,5]);
% title({'Mores Creek Summit: 04/05/2023','Hanover Reprocessed Snow Depth'})
title({strjoin(["Mores Creek Summit: ",num2str(outfnA(9:10)),"/",num2str(outfnA(11:12)),"/",num2str(outfnA(5:8))]),'LiDAR Snow Depth'})
xlabel('Easting (km)');ylabel('Northing (km)');
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
set(hI,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

subplot(1,2,2)
if isANN
    hI = imagesc(Coords.X./1000,Coords.Y./1000,LiDAR.A.ANN);
    h = colorbar; ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif','FontSize',14)
    daspect([1,1,1]);colormap(cmap);caxis(round([quantile(tmpA(:),[0.1,0.99])./0.05]).*0.05);%caxis([0,5]);
%     title({'Mores Creek Summit: 04/05/2023','ANN Gap Filled Snow Depth'})
    title({strjoin(["Mores Creek Summit: ",num2str(outfnA(9:10)),"/",num2str(outfnA(11:12)),"/",num2str(outfnA(5:8))]),'ANN Gap Filled Snow Depth'})
    xlabel('Easting (km)');ylabel('Northing (km)');
    set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
    set(hI,'AlphaData',~mask)
    set(gca,'color',0.75*[1 1 1]);
    saveas(FigANN, [writeDir,'\ANNsnowDepth.png'],'png');
elseif isRF
    hI = imagesc(Coords.X./1000,Coords.Y./1000,LiDAR.A.RF);
    h = colorbar; ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif','FontSize',14)
    daspect([1,1,1]);colormap(cmap);caxis(round([quantile(tmpA(:),[0.1,0.99])./0.05]).*0.05);%caxis([0,5]);
%     title({'Mores Creek Summit: 04/05/2023','RF Gap Filled Snow Depth'})
    title({strjoin(["Mores Creek Summit: ",num2str(outfnA(9:10)),"/",num2str(outfnA(11:12)),"/",num2str(outfnA(5:8))]),'RF Gap Filled Snow Depth'})
    xlabel('Easting (km)');ylabel('Northing (km)');
    set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
    set(hI,'AlphaData',~mask)
    set(gca,'color',0.75*[1 1 1]);
    saveas(FigRF, [writeDir,'\RFsnowDepth.png'],'png');
elseif isMLR
    hI = imagesc(Coords.X./1000,Coords.Y./1000,LiDAR.A.MLR);
    h = colorbar; ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif','FontSize',14)
    daspect([1,1,1]);colormap(cmap);caxis(round([quantile(tmpA(:),[0.1,0.99])./0.05]).*0.05);%caxis([0,5]);
%     title({'Mores Creek Summit: 04/05/2023','MLR Gap Filled Snow Depth'})
    title({strjoin(["Mores Creek Summit: ",num2str(outfnA(9:10)),"/",num2str(outfnA(11:12)),"/",num2str(outfnA(5:8))]),'MLR Gap Filled Snow Depth'})
    xlabel('Easting (km)');ylabel('Northing (km)');
    set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
    set(hI,'AlphaData',~mask)
    set(gca,'color',0.75*[1 1 1]);
    saveas(FigMLR, [writeDir,'\MLRsnowDepth.png'],'png');
end


end
%% Write Processed Data!
if isWriteLiDARmat
    save([writeDir,'\',outfnA(1:13),'LiDAR.mat'],'LiDAR','-v7.3')
    save([writeDir,'\',outfnA(1:13),'Coords.mat'],'Coords','-v7.3')
end
if isWriteGeoTiff
    % Need to Somehow Automate EPSG Code Lookup..
%         epsgCode = 32611;
    % Automated EPSG Code Lookup
    info = geotiffinfoT8(fullfilename);
    epsgCode = info.GeoTIFFCodes.PCS;
    geotiffwrite([writeDir,'\',outfnA,'_MLRgapfilled.tif'],LiDAR.A.MLR,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'\',outfnB,'_MLRgapfilled.tif'],LiDAR.B.MLR,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'\',outfnA,'_RFgapfilled.tif'],LiDAR.A.RF,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'\',outfnB,'_RFgapfilled.tif'],LiDAR.B.RF,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'\',outfnA,'_ANNgapfilled.tif'],LiDAR.A.ANN,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'\',outfnB,'_ANNgapfilled.tif'],LiDAR.B.ANN,Coords.R,'CoordRefSysCode',epsgCode);
end
%% Shad's Gap Filling Excercise
keyboard % Stop the Code here.. Proceed if Patient/Curious
for ii = 1:3
% Algortihm
if ii == 1
    isMLR = 1;
    isRF = 0;
    isANN = 0;
end
if ii == 2
    isMLR = 0;
    isRF = 1;
    isANN = 0;
end
if ii == 3
    isMLR = 0;
    isRF = 0;
    isANN = 1;
end
ShadA = tmpA;
% Create a Grid of NaN Values
% nanGrid = ones(size(ShadA));
% Create 100 m gaps every 1000 m
shadIx = [950.*2:1050.*2,1950.*2:2050.*2,2950.*2:2050.*2,2950.*2:3050.*2,...
    3950.*2:4050.*2,4950.*2:5050.*2,5950.*2:6050.*2,6950.*2:7050.*2];
ShadA(shadIx,shadIx) = NaN; % Checker Board Gaps
ShadA(shadIx,:) = NaN;ShadA(:,shadIx) = NaN; % Stripes
tmpix = find(isnan(tmpA));
tmpix2 = find(isnan(ShadA));
shadIx = setdiff(tmpix2,tmpix);
clear("tmpix","tmpix2")

% Machine Learning
ixOG = find(~isnan(ShadA));
% gapIx = find(isnan(ShadA));
% maskix = find(mask);
% gapIx = setdiff(gapIx,maskix);
gapIx = shadIx;
nEnsemble = 5;
% Machine Learning Infill Data Set
MLgap = [LiDAR.C.C(gapIx),LiDAR.C.aspect(gapIx),LiDAR.C.slope(gapIx),...
    LiDAR.C.gradN(gapIx),LiDAR.C.gradE(gapIx),LiDAR.C.aspectN(gapIx),LiDAR.C.aspectE(gapIx),...
    LiDAR.C.northness(gapIx),LiDAR.C.eastness(gapIx),LiDAR.D.D(gapIx),LiDAR.E.E(gapIx),LiDAR.F.F(gapIx)];
MLgap = (MLgap - predictorMean)./predictorStd;
[rows, ~] = find(isnan(MLgap));
rmvIx = unique(rows);
MLgap(rmvIx,:) = [];
shadIx(rmvIx) = [];
snowDepthGap = zeros(numel(shadIx),nEnsemble);
gapIx = shadIx;


rng(0) % For Reproducibility
for kk = 1:nEnsemble
    tic
    display(num2str(kk))
ix = datasample(ixOG,round(.0005.*numel(ixOG)),'Replace',false);
% Machine Learning Training Data Set
ML = [ShadA(ix),LiDAR.C.C(ix),LiDAR.C.aspect(ix),LiDAR.C.slope(ix),...
    LiDAR.C.gradN(ix),LiDAR.C.gradE(ix),LiDAR.C.aspectN(ix),LiDAR.C.aspectE(ix),...
    LiDAR.C.northness(ix),LiDAR.C.eastness(ix),LiDAR.D.D(ix),LiDAR.E.E(ix),LiDAR.F.F(ix)];
ML (:,2:end) = (ML(:,2:end) - predictorMean)./predictorStd;
[rows, ~] = find(isnan(ML));
rmvIx = unique(rows);
ML(rmvIx,:) = [];
if isANN
    [trainedModel, validationRMSE] = trainMCSsnowdepthANN(ML);
elseif isRF
    [trainedModel, validationRMSE] = trainRandomForestMCSsnowDepth(ML);
elseif isMLR
    [trainedModel, validationRMSE] = trainMCSsnowdepthMLR(ML);
end

snowDepthGap(:,kk) = trainedModel.predictFcn(MLgap);
toc

end

% Gap Fill Snow Depth and Snow Surface Elevation.
snowDepthGapMean = mean(snowDepthGap,2);
% snowDepthGapStd = std(snowDepthGap,[],2);
tmpShadA = ShadA;
ShadA(gapIx) = snowDepthGapMean;

% Validation
randIx = datasample(1:numel(shadIx),round(0.01.*numel(shadIx)),Replace=false);
% figure(); histogram(ShadA(shadIx));hold on; histogram(tmpA(shadIx));
figure();scatter(tmpA(shadIx(randIx)),ShadA(shadIx(randIx)),5,[tmpA(shadIx(randIx))-ShadA(shadIx(randIx))],'filled');colormap(cmap)
RMSE(ii) = sqrt(mean((ShadA(shadIx)-tmpA(shadIx)).^2));%./mean(tmpA(shadIx))
ME(ii) = mean(ShadA(shadIx)-tmpA(shadIx));
tmpR = corrcoef(tmpA(shadIx),ShadA(shadIx)).^2;
R2(ii) = tmpR(2,1);


%% Shad Gaps Figure
alph=isnan(tmpShadA);
figure();
% subplot(1,2,1)
hI = imagesc(Coords.X./1000,Coords.Y./1000,tmpShadA);
h = colorbar; ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif')
daspect([1,1,1]);colormap(cmap);caxis([0,5]);
title({'Mores Creek Summit: 04/05/2023','Hanover Reprocessed Snow Depth'})
xlabel('Easting (km)');ylabel('Northing (km)');
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
set(hI,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

figure()
subplot(1,2,1)
alph=isnan(ShadA);
hI = imagesc(Coords.X./1000,Coords.Y./1000,ShadA);
h = colorbar; ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif')
daspect([1,1,1]);colormap(cmap);caxis([0,5]);
if isANN
    title({'Mores Creek Summit: 04/05/2023','ANN Gap Filled Snow Depth'})
elseif isRF
    title({'Mores Creek Summit: 04/05/2023','RF Gap Filled Snow Depth'})
elseif isMLR
    title({'Mores Creek Summit: 04/05/2023','MLR Gap Filled Snow Depth'})
end
xlabel('Easting (km)');ylabel('Northing (km)');
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
set(hI,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

% Difference Figure
% figure();
subplot(1,2,2)
alph=isnan(tmpA);
hI = imagesc(Coords.X./1000,Coords.Y./1000,tmpA-ShadA);
h = colorbar; ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif')
daspect([1,1,1]);colormap(cmap);caxis([-1,1]);
if isANN
    title({'Mores Creek Summit: 04/05/2023','ANN Gap Filled Snow Depth'})
elseif isRF
    title({'Mores Creek Summit: 04/05/2023','RF Gap Filled Snow Depth'})
elseif isMLR
    title({'Mores Creek Summit: 04/05/2023','MLR Gap Filled Snow Depth'})
end
xlabel('Easting (km)');ylabel('Northing (km)');
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
set(hI,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
end