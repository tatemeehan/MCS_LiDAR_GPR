%% upscale LiDAR
clear;close all; clc;
%% Load Colormap
cmap = csvread('./colormaps/RdYlBu.csv');
cmap = flipud(cmap);
%% Processing Steps
% Write Output as GeoTiff
isWriteGeoTiff = 1;
% Write LiDAR as .mat
isWriteLiDARmat = 1;
% Data Directory 
dataDir = '/bsushare/hpmarshall-shared/LiDAR-GPR/20240315/';
% aux DataDir
auxDir = '/bsuhome/tatemeehan/git-repo/auxData/';
% Write Directory
writeDir = dataDir;

%% Load Raster Layers
% A Snow Depth
filename = '20240315_MCS-snowdepth_RFgapfilled.tif';
outfnA = filename(1:end-4);
fullfilename = fullfile(auxDir,filename);
[A,RA,~,~,~,~,~,~,epsgCode] = readLidarTif(fullfilename);
% tmpix = ~isnan(A);
% threshold = quantile(A(tmpix),0.999);% MAX snowdepth
% % threshold = 5.25; % MAX snowdepth
% A(A>threshold) = NaN;
% A(A<0) = NaN;
scale = 0.1;
[A,RA] = mapresize(A,RA,scale);
[xq,yq] = pixcenters(RA,size(A));
% xq = RA.XWorldLimits(1):RA.CellExtentInWorldX:RA.XWorldLimits(2);
% yq = RA.YWorldLimits(1):RA.CellExtentInWorldY:RA.YWorldLimits(2);
[Xq,Yq] = meshgrid(xq,yq);
% Yq = flipud(Yq);
% Get UTM Coordinates as Vector
% Get X,Y MeshGrid like Matrix
X = ones(RA.RasterSize(1),1)*linspace(RA.XWorldLimits(1),RA.XWorldLimits(2),RA.RasterSize(2));
Y = linspace(RA.YWorldLimits(1),RA.YWorldLimits(2),RA.RasterSize(1))'*ones(1,RA.RasterSize(2));
[lat,lon] = projinv(RA.ProjectedCRS,X,Y); 
utmprojection = projcrs(epsgCode);
[utmX,utmY] = projfwd(utmprojection,lat,lon);
utmY = flipud(utmY);
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

% Create Georeference
latlim = [min(Coords.lat(:)),max(Coords.lat(:))];
lonlim = [min(Coords.lon(:)) max(Coords.lon(:))];
sizeLidar = size(A);
georef = georefpostings(latlim,lonlim,sizeLidar,'RowsStartFrom','west','ColumnsStartFrom','north');
clear("latlim","lonlim","sizeLidar")
Coords.georef = georef;
Coords.epsgCode = epsgCode;


% B Snow Elevation
filename = '20240315_MCS-snow_RFgapfilled.tif';
outfnB = filename(1:end-4);
fullfilename = fullfile(auxDir,filename);
[B,RB,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
B = mapinterp(B,RB,Xq,Yq);
% Resize Snow On
% B = imresize(B,size(A));
% nanix = find(B<0 | B > 3000);
% B(nanix) = nan; 
clear('RB')

% C Reference DEM
dataDir2 = '/bsuhome/tatemeehan/git-repo/auxData';
% filename = '20220317_MCS-refDEM-T8.tif';
filename = 'MCS_REFDEM_WGS84.tif';
fullfilename = fullfile(auxDir,filename);
[C,RC,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
C(C<0) = NaN;
% NaN Mask
% border = csvread([auxDir,'MCSborder.csv']);
load([auxDir,'MCSborder.mat']);
mask = isnan(C)+border;
mask = mask>0;
mask = imresize(mask,size(A),Method="nearest");
clear('border')
C = mapinterp(C,RC,Xq,Yq);

% D Canopy Height
% filename = '20220317_MCS-canopyheight-T8.tif';
% filename = 'hanover_reprocess-canopy.tif';
filename = 'chm_mcs_1m.tif';
fullfilename = fullfile(auxDir,filename);
[D,RD,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
D = mapinterp(D,RD,Xq,Yq);
zedix = find(D<0 & D>=-1e1);
D(zedix) = 0;
nanix = find(D<-1e1 | D > 1e3);
D(nanix) = nan;
clear('RD','zedix')

% E Canopy Proximity
filename = 'proximity_map_mcs_1m.tif';
fullfilename = fullfile(auxDir,filename);
[E,RE,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
E = mapinterp(E,RE,Xq,Yq);
clear('RE')

% F Normalized Burn Ratio
filename =  'dnbr_mcs_utm_11n.tif';
fullfilename = fullfile(auxDir,filename);
[F,RF,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
% Trim NBR
[F,RF] = mapcrop(F,RF,RA.XWorldLimits,RA.YWorldLimits);
% Resize/Scale NBR to 0.5m
scale = round(RF.CellExtentInWorldX./RA.CellExtentInWorldX);
[F,RF] = mapresize(F,RF,scale,'cubic');
[F,~] = mapcrop(F,RF,RA.XWorldLimits,RA.YWorldLimits);
F = imresize(F,size(A));
clear('RF','RA','Xq','Yq','xq','yq')
%% Intelligent Gap Filling
% Inpaint Nans on MCS Snow Depth Predictors
D = inp8nt_nans(D,mask,5);
% Vegheight Threshold
Vix = D>0.5;
E = inp8nt_nans(E,mask,5);
F = inp8nt_nans(F,mask,5);

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

%% Compute Derivatives
% % Snow Depth Derivatives
[LiDAR.A.aspect,LiDAR.A.slope,LiDAR.A.gradN,LiDAR.A.gradE] = gradientm(LiDAR.A.A,georef);
LiDAR.A.aspect = inp8nt_nans(LiDAR.A.aspect,LiDAR.mask,5);
LiDAR.A.slope = inp8nt_nans(LiDAR.A.slope,LiDAR.mask,5);
LiDAR.A.gradN = inp8nt_nans(LiDAR.A.gradN,LiDAR.mask,5);
LiDAR.A.gradE = inp8nt_nans(LiDAR.A.gradE,LiDAR.mask,5);

% Aspect North Normalized
LiDAR.A.aspectN = cosd(LiDAR.A.aspect);
% Aspect East Normalized
LiDAR.A.aspectE = sind(LiDAR.A.aspect);
% Northness
LiDAR.A.northness = cosd(LiDAR.A.aspect).*sind(LiDAR.A.slope);
% Eastness
LiDAR.A.eastness = sind(LiDAR.A.aspect).*sind(LiDAR.A.slope);
 
% % Snow-On Elevation Derivatives
[LiDAR.B.aspect,LiDAR.B.slope,LiDAR.B.gradN,LiDAR.B.gradE] = gradientm(LiDAR.B.B,georef);
LiDAR.B.aspect = inp8nt_nans(LiDAR.B.aspect,LiDAR.mask,5);
LiDAR.B.slope = inp8nt_nans(LiDAR.B.slope,LiDAR.mask,5);
LiDAR.B.gradN = inp8nt_nans(LiDAR.B.gradN,LiDAR.mask,5);
LiDAR.B.gradE = inp8nt_nans(LiDAR.B.gradE,LiDAR.mask,5);

% Aspect North Normalized
LiDAR.B.aspectN = cosd(LiDAR.B.aspect);
% Aspect East Normalized
LiDAR.B.aspectE = sind(LiDAR.B.aspect);
% Northness
LiDAR.B.northness = cosd(LiDAR.B.aspect).*sind(LiDAR.B.slope);
% Eastness
LiDAR.B.eastness = sind(LiDAR.B.aspect).*sind(LiDAR.B.slope);

% Bare Ground Elevation Derivatives
[LiDAR.C.aspect,LiDAR.C.slope,LiDAR.C.gradN,LiDAR.C.gradE] = gradientm(LiDAR.C.C,georef);
LiDAR.C.aspect = inp8nt_nans(LiDAR.C.aspect,LiDAR.mask,5);
LiDAR.C.slope = inp8nt_nans(LiDAR.C.slope,LiDAR.mask,5);
LiDAR.C.gradN = inp8nt_nans(LiDAR.C.gradN,LiDAR.mask,5);
LiDAR.C.gradE = inp8nt_nans(LiDAR.C.gradE,LiDAR.mask,5);
% Aspect North Normalized
LiDAR.C.aspectN = cosd(LiDAR.C.aspect);
% Aspect East Normalized
LiDAR.C.aspectE = sind(LiDAR.C.aspect);
% Northness
LiDAR.C.northness = cosd(LiDAR.C.aspect).*sind(LiDAR.C.slope);
% Eastness
LiDAR.C.eastness = sind(LiDAR.C.aspect).*sind(LiDAR.C.slope);

%% Write Pre-Processed Data!
if isWriteLiDARmat
    save([writeDir,outfnA(1:13),'LiDAR5m.mat'],'LiDAR','-v7.3')
    save([writeDir,outfnA(1:13),'Coords5m.mat'],'Coords','-v7.3')
end
if isWriteGeoTiff
    % Need to Somehow Automate EPSG Code Lookup..
%         epsgCode = 32611;
    % Automated EPSG Code Lookup
    % info = geotiffinfoT8([auxDir,'20240315_MCS-snowdepth_RFgapfilled.tif']);
    % epsgCode = info.GeoTIFFCodes.PCS;
    geotiffwrite([writeDir,'20240315_MCS-snowdepth_RF_5m.tif'],LiDAR.A.A,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'20240315_MCS-snow_RF_5m.tif'],LiDAR.B.B,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'MCS_REF_DEM_5m.tif'],LiDAR.C.C,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'MCS_REF_CanopyHeight_5m.tif'],LiDAR.D.D,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'MCS_REF_CanopyProx_5m.tif'],LiDAR.E.E,Coords.R,'CoordRefSysCode',epsgCode);
    geotiffwrite([writeDir,'MCS_REF_dNBR_5m.tif'],LiDAR.F.F,Coords.R,'CoordRefSysCode',epsgCode);
end
