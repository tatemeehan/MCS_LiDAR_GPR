%% carSAR 
clear; close all; clc
%% Launch ParPool
parpool('local',4);
%% Read carSAR Phase and Coherence Data
dataDir = 'E:\MCS\MCS031924\CarSAR\processed\';
phzFn = '2024_03_19am_2024_03_19pm_utm_HH_int.tif';
corFn = '2024_03_19am_2024_03_19pm_utm_HH_cor.tif';
% InSAR Data
% Phase
[phz,phzR,~,~,lonPhz,latPhz,utmXphz,utmYphz] = readLidarTif([dataDir,phzFn]);
phz(phz>pi | phz <-pi) = NaN;
phz = medfilt2(phz,[10,10]);
% Coherence
[cor,corR,~,~,lonCor,latCor,utmXcor,utmYCor] = readLidarTif([dataDir,corFn]);
cor(cor>1 | cor<0) = NaN;
cor = medfilt2(cor,[10,10]);
%% Load Snow Depth & Density
dataDir = 'E:\MCS\MCS031524\';
bareDataDir = 'E:\MCS\MCS040623\LiDAR\';
densityFn = '20240315_MCS_GPRdensity.tif';
depthFn = '20240315_MCS-snowdepth_RFgapfilled.tif';
bareFn = 'MCS_REFDEM_WGS84.tif';
[density,R,~,~,lon,lat,utmX,utmY] = readLidarTif([dataDir,densityFn]);
[depth,~,~,~,~,~,~,~] = readLidarTif([dataDir,depthFn]);
[elvis,~,~,~,~,~,~,~] = readLidarTif([bareDataDir,bareFn]);
% Calculate Slope Angle
% Create Georeference
latlim = [min(lat(:)),max(lat(:))];
lonlim = [min(lon(:)) max(lon(:))];
sizeLidar = size(elvis);
georef = georefpostings(latlim,lonlim,sizeLidar,'RowsStartFrom','west','ColumnsStartFrom','north');
[~,slope,~,~] = gradientm(elvis,georef);

%% Interpolate Phase and Coherence to LiDAR Grid
% Get Phase and Coherence on EXACT same Grid
% Interpolate to M3 Grid
% xq = phzR.XWorldLimits(1):phzR.CellExtentInWorldX:phzR.XWorldLimits(2);
% yq = phzR.YWorldLimits(1):phzR.CellExtentInWorldY:phzR.YWorldLimits(2);
[xq,yq] = pixcenters(phzR,size(phz));
[X,Y] = meshgrid(xq,yq);
[Xq,Yq] = meshgrid(xq,yq);
% Yq = flipud(Yq);
cor = mapinterp(cor,corR,Xq,Yq);
cor = imresize(cor,size(phz));

% Interpolate Depth and Density and Slope
density = mapinterp(density,R,Xq,Yq);
density = imresize(density,size(phz));
depth = mapinterp(depth,R,Xq,Yq);
depth = imresize(depth,size(phz));
slope = mapinterp(slope,R,Xq,Yq);
slope = imresize(slope,size(phz));
%% Calculate Phase Change Response
f = 1.5; % GHz
c = 0.3; % speedolight
lambda = c./f;
thetai = 90-slope; % Incidence Angle with Zero Degree Brackets (Approx)
deltaPerm = ((phz.*lambda)./(4.*pi.*(depth))+cosd(thetai)).^2+sind(thetai).^2;
lwc = 0:0.0001:.1;
amLWC = 0;
amPerm = (c./WetCrimVRMS(density,amLWC)).^2;
pmPerm = amPerm+deltaPerm;
deltaLWC = zeros(size(phz));
% Solve for LWC Change
tic
p = gcp('nocreate');
if isempty(p)
    % There is no parallel pool
    poolsize = 0;
else
    % There is a parallel pool of <p.NumWorkers> workers
    poolsize = p.NumWorkers;
end
if poolsize == 0
    for ii = 1:numel(phz)
        calPerm = (c./WetCrimVRMS(density(ii),lwc)).^2;
        [~,minIx] = min(abs(calPerm(:)-pmPerm(ii)));
        deltaLWC(ii) = lwc(minIx).*100;
    end
else
    parfor (ii = 1:numel(phz),poolsize)
        calPerm = (c./WetCrimVRMS(density(ii),lwc)).^2;
        [~,minIx] = min(abs(calPerm(:)-pmPerm(ii)));
        deltaLWC(ii) = lwc(minIx).*100;
    end
end
toc

cmap = csvread('.\colormaps\RdYlBu.csv');
cmap = flipud(cmap);
%% Create Figures
% Phase
figure();
imagesc(utmXphz(1,:)./1000,(utmYphz(:,1))./1000,phz);daspect([1,1,1]);colormap("bone");hc=colorbar;
ylabel(hc,'\Delta \phi (rad)','fontname','serif','fontweight','bold','fontsize',12)
xlabel('Easting (km)');ylabel('Northing (km)');
% clim([-pi pi])
caxis([-pi pi])

set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
% Coherence
figure();
imagesc(utmXphz(1,:)./1000,(utmYphz(:,1))./1000,cor);daspect([1,1,1]);colormap("bone");hc=colorbar;
ylabel(hc,'Coherence','fontname','serif','fontweight','bold','fontsize',12)
xlabel('Easting (km)');ylabel('Northing (km)');
% clim([0 1])
caxis([0 1])
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
% LWC
figure();
imagesc(utmXphz(1,:)./1000,(utmYphz(:,1))./1000,deltaLWC);daspect([1,1,1]);colormap("bone");hc=colorbar;
ylabel(hc,'\Delta LWC (%)','fontname','serif','fontweight','bold','fontsize',12)
xlabel('Easting (km)');ylabel('Northing (km)');
clim([quantile(deltaLWC(deltaLWC>0),[0.05,0.95])])
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)