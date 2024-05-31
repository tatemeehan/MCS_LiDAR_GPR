%% LIDAR -- GPR .m
clear;close all; clc;
cmap = csvread('./colormaps/RdYlBu.csv');
cmap = flipud(cmap);
%% Read LiDAR .tif
% Load LiDAR as .mat Generated From postProLiDAR
isLoadLiDARmat = 1;
% Write Output as GeoTiff
isWriteGeoTiff = 0;
% Write LiDAR as .mat
isWriteLiDARmat = 0;
% Write LiDAR - GPR .mat
isWriteLidarGPRmat = 1;
% Write GPR Transect Data
isGPRcsv = 1;
% Calculate kdTree
isKDtree = 1;
% Data Directory 
dataDir = '/bsushare/hpmarshall-shared/LiDAR-GPR/20240315/';
% Write Directory
writeDir = dataDir;
if isLoadLiDARmat
    %% Load Ancillary Data Files from LiDAR_GPR.mat
    % Coordinate Data
    load([dataDir,'20240315_MCS-Coords5m.mat'])
    % LiDAR Data
    load([dataDir,'20240315_MCS-LiDAR5m.mat'])

    % Create Georeference
    latlim = [min(Coords.lat(:)),max(Coords.lat(:))];
    lonlim = [min(Coords.lon(:)) max(Coords.lon(:))];
    sizeLidar = size(LiDAR.A.A);
    georef = georefpostings(latlim,lonlim,sizeLidar,'RowsStartFrom','west','ColumnsStartFrom','north');
    clear("latlim","lonlim","sizeLidar")
%     LiDAR.georef = georef;
isCalcDerivatives = all([~isfield(LiDAR.A,'RFaspect'), ~isfield(LiDAR.A,'northness')]);
if isCalcDerivatives
    % Compute Derivatives using Random Forest Gap Filled Layer
    % Snow Depth Derivatives
    [LiDAR.A.RFaspect,LiDAR.A.RFslope,LiDAR.A.RFgradN,LiDAR.A.RFgradE] = gradientm(LiDAR.A.RF,georef);
    % Cure NaNs
    LiDAR.A.RFaspect = inp8nt_nans(LiDAR.A.RFaspect,LiDAR.mask);
    LiDAR.A.RFslope = inp8nt_nans(LiDAR.A.RFslope,LiDAR.mask);
    LiDAR.A.RFgradN = inp8nt_nans(LiDAR.A.RFgradN,LiDAR.mask);
    LiDAR.A.RFgradE = inp8nt_nans(LiDAR.A.RFgradE,LiDAR.mask);
    % Median Filter
    LiDAR.A.RFaspect = medfilt2(LiDAR.A.RFaspect,[5,5]);
    LiDAR.A.RFslope = medfilt2(LiDAR.A.RFslope,[5,5]);
    LiDAR.A.RFgradN = medfilt2(LiDAR.A.RFgradN,[5,5]);
    LiDAR.A.gradE = medfilt2(LiDAR.A.RFgradE,[5,5]);
    % Cure NaNs
    LiDAR.A.RFaspect = inp8nt_nans(LiDAR.A.RFaspect,LiDAR.mask);
    LiDAR.A.RFslope = inp8nt_nans(LiDAR.A.RFslope,LiDAR.mask);
    LiDAR.A.RFgradN = inp8nt_nans(LiDAR.A.RFgradN,LiDAR.mask);
    LiDAR.A.RFgradE = inp8nt_nans(LiDAR.A.RFgradE,LiDAR.mask);
    % Aspect North Normalized
    LiDAR.A.RFaspectN = cosd(LiDAR.A.RFaspect);
    % Aspect East Normalized
    LiDAR.A.RFaspectE = sind(LiDAR.A.RFaspect);
    % Northness
    LiDAR.A.RFnorthness = cosd(LiDAR.A.RFaspect).*sind(LiDAR.A.RFslope);
    % Eastness
    LiDAR.A.RFeastness = sind(LiDAR.A.RFaspect).*sind(LiDAR.A.RFslope);


    % Snow-On Elevation Derivatives
    [LiDAR.B.RFaspect,LiDAR.B.RFslope,LiDAR.B.RFgradN,LiDAR.B.RFgradE] = gradientm(LiDAR.B.RF,georef);
    % Cure NaNs
    LiDAR.B.RFaspect = inp8nt_nans(LiDAR.B.RFaspect,LiDAR.mask);
    LiDAR.B.RFslope = inp8nt_nans(LiDAR.B.RFslope,LiDAR.mask);
    LiDAR.B.RFgradN = inp8nt_nans(LiDAR.B.RFgradN,LiDAR.mask);
    LiDAR.B.RFgradE = inp8nt_nans(LiDAR.B.RFgradE,LiDAR.mask);
    % Median Filter
    LiDAR.B.RFaspect = medfilt2(LiDAR.B.RFaspect,[5,5]);
    LiDAR.B.RFslope = medfilt2(LiDAR.B.RFslope,[5,5]);
    LiDAR.B.RFgradN = medfilt2(LiDAR.B.RFgradN,[5,5]);
    LiDAR.B.gradE = medfilt2(LiDAR.B.RFgradE,[5,5]);
    % Cure NaNs
    LiDAR.B.RFaspect = inp8nt_nans(LiDAR.B.RFaspect,LiDAR.mask);
    LiDAR.B.RFslope = inp8nt_nans(LiDAR.B.RFslope,LiDAR.mask);
    LiDAR.B.RFgradN = inp8nt_nans(LiDAR.B.RFgradN,LiDAR.mask);
    LiDAR.B.RFgradE = inp8nt_nans(LiDAR.B.RFgradE,LiDAR.mask);
    % Aspect North Normalized
    LiDAR.B.RFaspectN = cosd(LiDAR.B.RFaspect);
    % Aspect East Normalized
    LiDAR.B.RFaspectE = sind(LiDAR.B.RFaspect);
    % Northness
    LiDAR.B.RFnorthness = cosd(LiDAR.B.RFaspect).*sind(LiDAR.B.RFslope);
    % Eastness
    LiDAR.B.RFeastness = sind(LiDAR.B.RFaspect).*sind(LiDAR.B.RFslope);
end
else

% A Snow Depth
filename = 'hanover_reprocess-snowdepth.tif';
fullfilename = fullfile(dataDir,filename);
[A,tiffR,~,~,lon,lat,utmX,utmY] = readLidarTif(fullfilename);
% Get UTM Coordinates as Vector
X = utmX(1,:);
Y = utmY(:,1);
Xi = utmX(:); 
Yi = utmY(:);
% Create Coordinate Structure
Coords.R = tiffR;
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
filename = '20220317_MCS-snow-T8.tif';
fullfilename = fullfile(dataDir,filename);
[B,~,~,~,~,~,~,~] = readLidarTif(fullfilename);

% C Reference DEM
% filename = '20220317_MCS-refDEM-T8.tif';
filename = 'MCS_REFDEM_WGS84.tif';
fullfilename = fullfile(dataDir,filename);
[C,~,~,~,~,~,~,~] = readLidarTif(fullfilename);

% D Canopy Height
filename = '20220317_MCS-canopyheight-T8.tif';
fullfilename = fullfile(dataDir,filename);
[D,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
% Vegheight Threshold
Vix = D>0.5;

% E Canopy Proximity
filename = '20220317_MCS-canopyproximity-T8.tif';
fullfilename = fullfile(dataDir,filename);
[E,~,~,~,~,~,~,~] = readLidarTif(fullfilename);

% F Normalized Burn Ratio
filename =  'MCS-NBR-T8.tif';
fullfilename = fullfile(dataDir,filename);
[F,~,~,~,~,~,~,~] = readLidarTif(fullfilename);

%% Process LiDAR Rasters
% Gap Fill
% A
nanIx = find(A(:) == -9999);
A(nanIx) = NaN;
A = inpaint_nans(A,5);
% B
nanIx = find(B(:) == -9999);
B(nanIx) = NaN;
B = inpaint_nans(B,5);
% C
nanIx = find(C(:) == -9999);
C(nanIx) = NaN;
C = inpaint_nans(C,5);
% D
nanIx = find(D(:) == -9999);
D(nanIx) = NaN;
D = inpaint_nans(D,5);
% E
nanIx = find(E(:) == -9999);
E(nanIx) = NaN;
E = inpaint_nans(E,5);
% F
nanIx = find(F(:) == -9999);
F(nanIx) = NaN;
F = inpaint_nans(F,5);

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
clear('A', "B", 'C', 'D', "E", "F","Vix")

% Median Smoothing
LiDAR.A.A = medfilt2(LiDAR.A.A,[5,5]);
LiDAR.B.B = medfilt2(LiDAR.B.B,[5,5]);
LiDAR.C.C = medfilt2(LiDAR.C.C,[5,5]);

% Create Georeference
latlim = [min(Coords.lat(:)),max(Coords.lat(:))];
lonlim = [min(Coords.lon(:)) max(Coords.lon(:))];
sizeLidar = size(LiDAR.A.A);
georef = georefpostings(latlim,lonlim,sizeLidar,'RowsStartFrom','west','ColumnsStartFrom','north');
clear("latlim","lonlim","sizeLidar",'nanIx')

% Snow Depth Derivatives
[LiDAR.A.aspect,LiDAR.A.slope,LiDAR.A.gradN,LiDAR.A.gradE] = gradientm(LiDAR.A.A,georef);
LiDAR.A.aspect = inpaint_nans(LiDAR.A.aspect,5);

% Snow-On Elevation Derivatives
[LiDAR.B.aspect,LiDAR.B.slope,LiDAR.B.gradN,LiDAR.B.gradE] = gradientm(LiDAR.B.B,georef);
LiDAR.B.aspect = inpaint_nans(LiDAR.B.aspect,5);

% Bare Ground Elevation Derivatives
[LiDAR.C.aspect,LiDAR.C.slope,LiDAR.C.gradN,LiDAR.C.gradE] = gradientm(LiDAR.C.C,georef);
LiDAR.C.aspect = inpaint_nans(LiDAR.C.aspect,5);

    % Predictors
    if exist('smoothFlag')
    else
        smoothFlag = 1;
        LiDAR.A.aspect = medfilt2(LiDAR.A.aspect,[25,25]);
        LiDAR.A.slope = medfilt2(LiDAR.A.slope,[25,25]);
        LiDAR.A.gradN = medfilt2(LiDAR.A.gradN,[25,25]);
        LiDAR.A.gradE = medfilt2(LiDAR.A.gradE,[25,25]);
        LiDAR.B.aspect = medfilt2(LiDAR.B.aspect,[25,25]);
        LiDAR.B.slope = medfilt2(LiDAR.B.slope,[25,25]);
        LiDAR.B.gradN = medfilt2(LiDAR.B.gradN,[25,25]);
        LiDAR.B.gradE = medfilt2(LiDAR.B.gradE,[25,25]);
        LiDAR.C.aspect = medfilt2(LiDAR.C.aspect,[25,25]);
        LiDAR.C.slope = medfilt2(LiDAR.C.slope,[25,25]);
        LiDAR.C.gradN = medfilt2(LiDAR.C.gradN,[25,25]);
        LiDAR.C.gradE = medfilt2(LiDAR.C.gradE,[25,25]);
        
        % Normalize aspect between 0 and 1;
        LiDAR.A.northness = abs(LiDAR.A.aspect-180)./180;
        LiDAR.B.northness = abs(LiDAR.B.aspect-180)./180;
        LiDAR.C.northness = abs(LiDAR.C.aspect-180)./180;
%         LiDAR.A.northness = cosd(LiDAR.A.aspect+45)+sind(LiDAR.A.aspect+45);
%         LiDAR.B.northness = cosd(LiDAR.B.aspect+45)+sind(LiDAR.B.aspect+45);
%         LiDAR.C.northness = cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45);
    end
end
% Write LiDAR and Coods .mat
if isWriteLiDARmat
    save([writeDir,'20240315_MCS-LiDAR.mat'],'LiDAR','-v7.3')
    save([writeDir,'20240315_MCS-Coords.mat'],'Coords','-v7.3')
end
%% Read GPR .csv
filename = 'MCS031524-TWT.csv';
% filename = '/bsuhome/tatemeehan/git-repo/auxData/MCS2024-TWT.csv';
% Read GPR data
GPR = readtable(fullfile(dataDir,filename));
% GPR = readtable(filename);
gprX = GPR.Easting; gprY = GPR.Northing;
gprTWT = GPR.TWT; gprZ = GPR.ElevationWGS84;
% % For Standardized Data Set
% tmpdate = '20240213';
% tmpDateTime = datetime(tmpdate,'InputFormat','yyyyMMdd');
% [~,dateIx] = min(abs(GPR.DateTime - tmpDateTime));
% gprTWT = (GPR.zTWT.*GPR.stdTWT(dateIx))+GPR.meanTWT(dateIx);
% gprTWT = (GPR.zTWT.*GPR.stdTWT)+GPR.meanTWT(dateIx);


%% Co-Locate LiDAR and GPR
% KD-tree Searcher
if isKDtree
tic
winsize = mean([Coords.R.CellExtentInWorldX,Coords.R.CellExtentInWorldY]);
mykdtree=KDTreeSearcher([GPR.Easting GPR.Northing]); % searcher for GPR
[IDX,D]=rangesearch(mykdtree,[Coords.Xi, Coords.Yi],winsize); % search the points in the LiDAR grid
% Remove Empty Cells
ix =  find(~cellfun(@isempty,D));
D = D(ix); IDX = IDX(ix);
kd.D = D; kd.IDX = IDX;kd.ix = ix;kd.winsize = winsize;
toc
% Save the Output
save([dataDir,'20240315_MCS-kdtree-5m.mat'],'kd','-v7.3')
% save([dataDir,'2024_MCS-kdtree.mat'],'kd','-v7.3')


else
% Load Previous KD-Tree
% load('MCS03172022kdtree.mat')
load([dataDir,'20240315_MCS-kdtree-5m.mat'])
% load([dataDir,'20240213_MCS-kdtree.mat'])
end

%% Pair LiDAR and GPR Observations
GPRx = zeros(length(kd.ix),1);GPRy = GPRx; GPRtwt = GPRx;GPRz = GPRtwt;
GPRtime = GPRtwt; GPRdate = GPRtwt;
    for jj = 1:length(kd.ix)
        % Calculate Distance
        dist = kd.D{jj};dist = dist(:);
        % Find Points within radius r of grid point
        tmpIx = kd.IDX{jj};tmpIx = tmpIx(:);
        % Median
        GPRtwt(jj) = median(gprTWT(tmpIx));
        if isnan(GPRtwt(jj))
                keyboard 
        end
        GPRz(jj) = median(gprZ(tmpIx));
        % Coordinates
        GPRx(jj) = Coords.Xi(kd.ix(jj));
        GPRy(jj) = Coords.Yi(kd.ix(jj));
%         % Date and Time
%         GPRtime(jj) = min(gprTime(tmpIx));
%         GPRdate(jj) = min(gprDate(tmpIx));
    end
    % Correct GPR TWT for Slope angle
    GPRtwt = GPRtwt./cosd(LiDAR.B.slope(kd.ix));
    % GPRtwt = GPRtwt./cosd(LiDAR.B.RFslope(kd.ix));

% House Keeping
lidarGPR.kd = kd;
lidarGPR.X = GPRx; lidarGPR.Y = GPRy; lidarGPR.Z = GPRz;lidarGPR.TWT = GPRtwt;
clear ('gprTWT','gprX','gprY','gprZ','GPRtwt','GPRx','GPRy','GPRz','GPRdate','GPRtime','tmpIx','dist','IDX','ix','D','winsize')

%% Calculate LiDAR - GPR Density
% Load Bulk Density Data
isBiasCorrectTWT = 0;
isInsitu = 1;
if isBiasCorrectTWT
if isInsitu 
insitu = readtable([dataDir,'MCS20240213_SWE.csv']);
% Find GPR Loctions within radius of Density Obs
r = 100;
ix = [];
for kk = 1:numel(lidarGPR.X)
dist = sqrt((insitu.UTM_X-lidarGPR.X(kk)).^2+(insitu.UTM_Y-lidarGPR.Y(kk)).^2);
if min(dist) > r
    tmpdensity(kk) = nan;
else
            % Inverse Distance Weighting
            p = 1;
        w = 1./dist.^p; w(isinf(w))= 1;
        w=w(:);
    tmpdensity(kk) = sum(w.*insitu.Density_kgm3)./sum(w);
end
% [minIx] = find(dist<=r);
% ix = [ix;minIx];
end
ix = find(~isnan(tmpdensity));
% !!Bias Correction!!
% 1.
avgTWT = median(lidarGPR.TWT(ix));
avgDepth = median(LiDAR.A.RF(kd.ix(ix)));
calculatedV = mean(DryCrimVRMS(tmpdensity(ix)));
deltaTWT = (2.*avgDepth./calculatedV) - avgTWT;
% 2. 
% calculatedV = (DryCrimVRMS(tmpdensity(ix)));
% avgTWT = (lidarGPR.TWT(ix));
% avgDepth = (LiDAR.A.RF(kd.ix(ix)));
% deltaTWT = median((2.*avgDepth(:)./calculatedV(:)) - avgTWT(:));
lidarGPR.TWT = lidarGPR.TWT+deltaTWT;
else
% Input Average Density
avgRho = 300;%mean(insitu.Density_kgm3); % Prior Density % 02132024 Avg Density 300
avgV = DryCrimVRMS(avgRho);
avgTWT = median(lidarGPR.TWT(ix));
avgDepth = avgV.*avgTWT./2;
deltaDepth = median(LiDAR.A.RF(kd.ix(ix))) - avgDepth;
deltaTWT = 2.*deltaDepth./avgV;
lidarGPR.TWT = lidarGPR.TWT+deltaTWT;
end
end
% Calculate Density
isDrySnow = 1;
if isDrySnow
% Dry Snow
% velocity = 2.*(LiDAR.A.RF(kd.ix)./lidarGPR.TWT);
velocity = 2.*(LiDAR.A.A(kd.ix)./lidarGPR.TWT);
lidarGPR.Density = DryCrim(velocity).*1000;
else
% Wet Snow
% velocity = 2.*(LiDAR.A.RF(kd.ix)./lidarGPR.TWT);
velocity = 2.*(LiDAR.A.A(kd.ix)./lidarGPR.TWT);

lwc = 0.0025;
lidarGPR.Density = WetCrim(velocity,lwc ).*1000;
end
infIx = find(isinf(lidarGPR.Density));
lidarGPR.Density(infIx) = quantile(lidarGPR.Density,0.99);

% Raw Data Quantiles 10% and 90% thresholds
Qdensity = quantile(lidarGPR.Density,[0.25,0.5,0.75]);
% Outlier Threshold
% This is happenchance near the IQR
% IQR Threshold
threshHi = lidarGPR.Density>Qdensity(3);
threshLo = lidarGPR.Density<Qdensity(1);
threshIx = find(threshHi| threshLo); % For Threshold
% not used
notThreshIx = 1:numel(lidarGPR.Density);
notThreshIx(threshIx) = [];
% Allocate Co-Located Depths
lidarGPR.Depth = zeros(size(kd.ix));

%% KDtreesearcher Smooth GPR-Lidar-Density
% Simple removal of Outliers
tmpLidarGPRdensity = lidarGPR.Density;
% Median Interpolation of Outliers
winsize = 12.5;
mykdtree=KDTreeSearcher([lidarGPR.X, lidarGPR.Y]); % searcher for Colocation
[IDXp,Dp]=rangesearch(mykdtree,[lidarGPR.X lidarGPR.Y],winsize); % search the points in the Colocated Domain
tmpBin = find(~cellfun(@isempty,IDXp));
IDXp = IDXp(tmpBin);Dp=Dp(tmpBin);
nanix = [];
ixUsed = zeros(length(lidarGPR.Density),1);
% Weighted Average Smoothing
    for jj = 1:length(lidarGPR.Density)
                % Calculate Distance
        dist = Dp{jj};dist = dist(:);
        % Find Points within radius r of grid point
        tmpIx = IDXp{jj};tmpIx = tmpIx(:);
        badIx = ismember(tmpIx,threshIx);
        tmpIx = tmpIx(~badIx); dist = dist(~badIx);
        % if length(tmpIx)<5%isempty(tmpIx)
        if isempty(tmpIx)
            lidarGPR.Density(jj) = NaN;
            lidarGPR.TWT(jj) = NaN;
            lidarGPR.Depth(jj) = NaN;
            nanix = [nanix;jj];
        else
        % Inverse Distance Weighting
%         w = 1./dist.^p; w(isinf(w))= 1;
        % calculate weights - bisquare kernel
%         w=15/16*(1-((dist./winsize)+eps).^2).^2;
%         w=w(:);
%         smoothLidarDensity(jj) = sum(w.*LidarDensityNorm(tmpIx))./sum(w);
%         smoothLidarDensity(jj) = sum(w.*LidarGPRdensitysqz(tmpIx))./sum(w);

        % Median Smoothing
        tmpDensity = median(tmpLidarGPRdensity(tmpIx));
        [~,tmpDensityIx] = min(abs(tmpDensity-tmpLidarGPRdensity(tmpIx)));
        ixUsed(jj) = tmpIx(tmpDensityIx);
        lidarGPR.Density(jj) = tmpLidarGPRdensity(tmpIx(tmpDensityIx));
        lidarGPR.TWT(jj) = lidarGPR.TWT(tmpIx(tmpDensityIx));
        lidarGPR.Depth(jj) = LiDAR.A.A(kd.ix(tmpIx(tmpDensityIx)));
        % lidarGPR.Depth(jj) = LiDAR.A.RF(kd.ix(tmpIx(tmpDensityIx)));
        end
    end
    %inpaintnans
    lidarGPR.Density = inpaint_nans(lidarGPR.Density,5);
    lidarGPR.TWT = inpaint_nans(lidarGPR.TWT,5);
    lidarGPR.Depth = inpaint_nans(lidarGPR.Depth,5);
    % Calculate SWE
    lidarGPR.SWE = lidarGPR.Depth.*lidarGPR.Density;
    unUsedIx = find(ixUsed == 0);
    usedIx = find(ixUsed);
    tmpUsedIx = interp1(usedIx, ixUsed(usedIx), unUsedIx, 'nearest');
    ixUsed(unUsedIx) = tmpUsedIx;
    % CRIM wavespeed
    lidarGPR.Velocity = DryCrimVRMS(lidarGPR.Density);
    % Dielectric Permittivity
    lidarGPR.Permittivity = (0.3./lidarGPR.Velocity).^2;
    % GPR Depth
    lidarGPR.depthGPR = lidarGPR.TWT.*lidarGPR.Velocity./2.*100;
    %% Export .csv of GPR Transects
    if isGPRcsv
        GPRout = table(lidarGPR.X,lidarGPR.Y,lidarGPR.TWT,lidarGPR.Depth.*100,lidarGPR.Density,lidarGPR.SWE,lidarGPR.Velocity,lidarGPR.Permittivity);
        GPRout = renamevars(GPRout,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8"], ...
            ["X (WGS84 UTM 11N)","Y (WGS84 UTM 11N)","TWT (ns)","Depth (cm)","Density (kg/m3)","SWE (mm)","Velocity (m/ns)","Permittivity"]);
        writetable(GPRout,[writeDir,'20240315_MCS-LiDAR-GPR-density_5m.csv'])
        clear('GPRout')
    end
    % Export LiDAR_GPR .mat
    if isWriteLidarGPRmat
        save([writeDir,'20240315_MCS-LiDAR-GPR_5m.mat'],'lidarGPR','-v7.3')
    end

    % House Keeping
    clear('GPRout','IDXp','Dp','nanix','ixUsed','unUsedIx','usedIx','badIx','dist','threshHi', 'threshLo','threshIx','notThreshIx','tmpLidarGPRdensity',...
        'Qdensity','infIx','tmpIx','tmpDensityIx','tmpDensity','tmpBin','tmpUsedIx')
