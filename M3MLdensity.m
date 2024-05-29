%% Compare ML Density to M3 Density
clear; close all; clc;
%% Launch ParPool
parpool('local',4);
%% Load ML & M3 data
% M3 Density Model
% M3dataDir = 'E:\MCS\MCS040623\M3';
% M3filename = 'snow_density_2023-04-06T230000.tif';
M3dataDir = 'E:\MCS\MCS031524\M3';
M3filename = 'snow_density_2024-03-15T23 00 00.tif';
M3fullfilename = fullfile(M3dataDir,M3filename);
[M3,M3R,~,~,M3lon,M3lat,M3utmX,M3utmY] = readLidarTif(M3fullfilename);
% Interpolate to M3 Grid
xq = M3R.XWorldLimits(1):M3R.CellExtentInWorldX:M3R.XWorldLimits(2);
yq = M3R.YWorldLimits(1):M3R.CellExtentInWorldY:M3R.YWorldLimits(2);
[Xq,Yq] = meshgrid(xq,yq);
Yq = flipud(Yq);
% Load .tif
% Snow Density
filename = '20240315_MCS_MLRdensity.tif';
MLdataDir = 'E:\MCS\MCS031524';
fullfilename = fullfile(MLdataDir,filename);
[A,RA,~,~,lon,lat,utmX,utmY] = readLidarTif(fullfilename);
% tmpix = ~isnan(A);
% threshold = quantile(A(tmpix),0.999);% MAX snowdepth
% % threshold = 5.25; % MAX snowdepth
% A(A>threshold) = NaN;
% A(A<0) = NaN;
% xq = RA.XWorldLimits(1):RA.CellExtentInWorldX:RA.XWorldLimits(2);
% yq = RA.YWorldLimits(1):RA.CellExtentInWorldY:RA.YWorldLimits(2);
% [Xq,Yq] = meshgrid(xq,yq);
% Yq = flipud(Yq);
% Store Results
MLresultsMLR.MLRi.Density = A;
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

% GPR Density
filename = '20240315_MCS_GPRdensity.tif';
MLdataDir = 'E:\MCS\MCS031524';
fullfilename = fullfile(MLdataDir,filename);
[B,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
MLresultsGPR.GPR.Density = B;

clear('A','B')
% % Load Coords
% LiDARdataDir = 'E:\MCS\MCS040623\LiDAR';
% filename = 'MCS040523-Coords';
% fullfilename = fullfile(LiDARdataDir,filename);
% load(fullfilename);
% 
% % ML MLR Density Model
% MLdataDir = 'E:\MCS\MCS040623\GPR';
% MLfilename = 'MCS040623-MLresultsMLR.mat';
% MLfullfilename = fullfile(MLdataDir,MLfilename);
% load(MLfullfilename);
% 
% % ML GPR Density Model
% MLdataDir = 'E:\MCS\MCS040623\GPR';
% MLfilename = 'MCS040623-MLresultsGPR.mat';
% MLfullfilename = fullfile(MLdataDir,MLfilename);
% load(MLfullfilename);


%% Upscale to 50 m
% Get pixel centers !!!!
% [x,y] = pixcenters(R,size(geotiff));
% [X,Y] = meshgrid(x,y);
MLresultsGPR.GPR.Density50m = mapinterp(MLresultsGPR.GPR.Density,Coords.R,Xq,Yq);
MLresultsMLR.MLRi.Density50m =  mapinterp(MLresultsMLR.MLRi.Density,Coords.R,Xq,Yq);
MLresultsGPR.GPR.Density50m = imresize(MLresultsGPR.GPR.Density50m ,size(M3),'bilinear');
MLresultsMLR.MLRi.Density50m = imresize(MLresultsMLR.MLRi.Density50m,size(M3),'bilinear');
% Trim Rasters to Boundary
nanIx = find(~isnan(MLresultsGPR.GPR.Density50m(:)));
[row,col] = ind2sub(size(MLresultsGPR.GPR.Density50m),nanIx);
minRowIx = min(row); maxRowIx = max(row);
minColIx = min(col); maxColIx = max(col);
x = xq(minRowIx:maxRowIx); y = yq(minColIx:maxColIx);
MLRi = MLresultsMLR.MLRi.Density50m(minRowIx:maxRowIx,minColIx:maxColIx);
GPR = MLresultsGPR.GPR.Density50m(minRowIx:maxRowIx,minColIx:maxColIx);
M3d = M3(minRowIx:maxRowIx,minColIx:maxColIx);

% MLresultsGPR.GPR.Density50m = imresize(MLresultsGPR.GPR.Density,1./100,'bilinear');
% MLresultsMLR.MLRi.Density50m = imresize(MLresultsMLR.MLRi.Density,1./100,'bilinear');
% % Find NaN Padding
% nanIx = find(~isnan(MLresultsGPR.GPR.Density50m(:)));
% % Resize Coordinates
% Coords.utmX50 = imresize(Coords.utmX,1./100,"nearest");
% Coords.utmY50 = imresize(Coords.utmY,1./100,"nearest");

% Compute moving Stats
[R50mGPR,RMSE50mGPR,ME50mGPR] = movCorr2D(M3d,GPR,5);
[R50mMLRi,RMSE50mMLRi,ME50mMLRi] = movCorr2D(M3d,MLRi,5);

%% Upscale to 5 m
% Interpolate to 5m LiDAR Grid
scale = 10;
xq = Coords.R.XWorldLimits(1):Coords.R.CellExtentInWorldX.*scale:Coords.R.XWorldLimits(2);
yq = Coords.R.YWorldLimits(1):Coords.R.CellExtentInWorldY.*scale:Coords.R.YWorldLimits(2);
[Xq,Yq] = meshgrid(xq,yq);
Yq = flipud(Yq);
yq = flipud(yq);
[MLresultsGPR.GPR.Density5m] = mapinterp(MLresultsGPR.GPR.Density,Coords.R,Xq,Yq);
MLresultsMLR.MLRi.Density5m =  mapinterp(MLresultsMLR.MLRi.Density,Coords.R,Xq,Yq);
% Moving window Statistics
[R5m,RMSE5m,ME5m] = movCorr2D(MLresultsGPR.GPR.Density5m,MLresultsMLR.MLRi.Density5m,11);
% MLresultsGPR.GPR.Density5m = imresize(MLresultsGPR.GPR.Density50m ,size(M3),'bilinear');
% MLresultsMLR.MLRi.Density5m = imresize(MLresultsMLR.MLRi.Density50m,size(M3),'bilinear');
% Trim Rasters to Boundary
% nanIx = find(~isnan(MLresultsGPR.GPR.Density50m(:)));
% [row,col] = ind2sub(size(MLresultsGPR.GPR.Density50m),nanIx);
% minRowIx = min(row); maxRowIx = max(row);
% minColIx = min(col); maxColIx = max(col);
% x = xq(minRowIx:maxRowIx); y = yq(minColIx:maxColIx);
% MLRi = MLresultsMLR.MLRi.Density50m(minRowIx:maxRowIx,minColIx:maxColIx);
% GPR = MLresultsGPR.GPR.Density50m(minRowIx:maxRowIx,minColIx:maxColIx);
% M3d = M3(minRowIx:maxRowIx,minColIx:maxColIx);
%% Find M3 Density within 50 m LiDAR Domain
ix = zeros(numel(Coords.utmX50(nanIx)),1);
for kk = 1:numel(Coords.utmX50(nanIx))
    dist = sqrt(((Coords.utmX50(nanIx(kk))-M3utmX(:)).^2)+((Coords.utmY50(nanIx(kk))-M3utmY(:)).^2));
    [~,ix(kk)] = min(dist);
end

% Extract Tile
for kk = 1:4
    if kk == 1
    dist = sqrt(((min(Coords.utmX50(:))-M3utmX(:)).^2)+((min(Coords.utmY50(:))-M3utmY(:)).^2));
    elseif kk == 2
            dist = sqrt(((min(Coords.utmX50(:))-M3utmX(:)).^2)+((max(Coords.utmY50(:))-M3utmY(:)).^2));
    elseif kk == 3
            dist = sqrt(((max(Coords.utmX50(:))-M3utmX(:)).^2)+((max(Coords.utmY50(:))-M3utmY(:)).^2));
    elseif kk == 4
            dist = sqrt(((max(Coords.utmX50(:))-M3utmX(:)).^2)+((min(Coords.utmY50(:))-M3utmY(:)).^2));
    end
    [~,bbix(kk)] = min(dist);
end
[row,col] = ind2sub(size(M3),bbix);
M3Tile = M3(min(row):max(row),min(col):max(col));
M3Tile = imresize(M3Tile,size(MLresultsMLR.MLRi.Density50m),'bilinear');
mask = find(isnan(MLresultsGPR.GPR.Density50m(:)));
M3Tile(mask) = nan;
%% Compute some Stats
% Find NaN Padding
nanIx = find(~isnan(MLRi(:)));
rmseMLR = sqrt(mean((MLRi(nanIx(:))-M3d(nanIx(:))).^2))
rmseGPR = sqrt(mean((GPR(nanIx(:))-M3d(nanIx(:))).^2))
meMLR = mean(M3d(nanIx(:))-(MLRi(nanIx(:))))
meGPR = mean(M3d(nanIx(:))-(GPR(nanIx(:))))
corrcoef(MLRi(nanIx(:)),M3d(nanIx(:))).^2
corrcoef(GPR(nanIx(:)),M3d(nanIx(:))).^2
[~,MLRiSSIM]= ssim(MLRi,M3d,'Radius',1.5);
nanmean(MLRiSSIM(:))
[~,GPRSSIM]= ssim(GPR,M3d,'Radius',1.5);
nanmean(GPRSSIM(:))


% rmse = sqrt(mean((MLresultsGPR.GPR.Density50m(nanIx(:))-M3(ix)).^2));
% rmse = sqrt(mean((MLresultsMLR.MLRi.Density50m(nanIx(:))-M3(ix)).^2));

corrcoef(MLresultsGPR.GPR.Density50m(nanIx(:)),M3(ix))
corrcoef(MLresultsMLR.MLRi.Density50m(nanIx(:)),M3(ix))
% SSIM
[MLRiSSIM,MLRiSSIMs] = ssim(MLresultsMLR.MLRi.Density50m,M3Tile,'Radius',1.5);
[GPRSSIM,GPRSSIMs] = ssim(MLresultsGPR.GPR.Density50m,M3Tile,'Radius',1.5);
% Moving Correlation
GPRcorr = movcorr(MLresultsGPR.GPR.Density50m,M3,5,'omitnan');

%% Make Some figures
cmap = csvread('.\colormaps\RdYlBu.csv');
cmap = flipud(cmap);

alph=isnan(GPR);

figure();
subplot(1,2,1)
plot(M3d(:),MLRi(:),'.k','MarkerSize',1);lsline;hold on
binscatter(M3d(:),MLRi(:),250);
colormap(gca,'gray')
xlim([250,500]);ylim([250,500])
lsline
grid on; grid minor;
xlabel('MLRi Density');ylabel('M3 Density')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
axis square
subplot(1,2,2)
plot(M3d(:),GPR(:),'.k','MarkerSize',1);lsline;hold on;
binscatter(M3d(:),GPR(:),250);
colormap(gca,'gray')
xlim([250,500]);ylim([250,500])
lsline;
grid on; grid minor;
xlabel('GPR Density');ylabel('M3 Density')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
axis square

% Rasters
figure(); 
subplot(1,2,1)
hM3MLRi = imagesc(x./1000,flipud(y)./1000,M3d - MLRi);colorbar
daspect([1,1,1]);
colormap(cmapAdapt(M3d - MLRi,cmap))
title('M3 - MLRi')
caxis([-50 150])
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
set(hM3MLRi,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
subplot(1,2,2)
hM3GPR = imagesc(x./1000,flipud(y)./1000,M3d - GPR);colorbar
daspect([1,1,1]);
caxis([-50 150])
colormap(cmapAdapt(M3d - GPR,cmap))
title('M3 - GPR')
set(hM3GPR,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
set(gca,'fontname','serif','fontweight','bold','fontsize',12)

figure();
subplot(1,3,1)
hM3 = imagesc(x./1000,flipud(y)./1000,M3d);ch = colorbar; 
ylabel(ch,'Density (kg/m^3)','fontname','serif','fontweight','bold','fontsize',12)
daspect([1,1,1]);
colormap(cmap)
title('M3')
% clim([quantile(M3d(:),[0.025,0.975])])
% clim([275 425])
% caxis([250,450])
clim([300,400])
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
set(hM3,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

subplot(1,3,2)
hMLRi = imagesc(x./1000,flipud(y)./1000,MLRi);ch = colorbar;
ylabel(ch,'Density (kg/m^3)','fontname','serif','fontweight','bold','fontsize',12)
daspect([1,1,1]);
colormap(cmap)
title('MLRi')
% clim([quantile(MLRi(:),[0.025,0.975])])
% clim([275 425])
% caxis([250,450])
clim([300,400])
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
set(hMLRi,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

subplot(1,3,3)
hGPR = imagesc(x./1000,flipud(y)./1000,GPR);ch = colorbar;
ylabel(ch,'Density (kg/m^3)','fontname','serif','fontweight','bold','fontsize',12)
daspect([1,1,1]);
colormap(cmap)
title('GP')
% clim([quantile(GPR(:),[0.025,0.975])])
% clim([275 425])
% caxis([250,450])
clim([300,400])
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
set(hGPR,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

%% Moving Statistics Figures
alph2 = isnan(MLresultsGPR.GPR.Density5m);
figure();
tiledlayout(3,3,"TileSpacing","compact")
% M3 vs. MLRi
% R
nexttile
h = imagesc(Coords.X./1000,Coords.Y./1000,R50mMLRi);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
clim([-1,1])
title('R: M3 -- MLRi')
colorbar
%RMSE
nexttile
h = imagesc(Coords.X./1000,Coords.Y./1000,RMSE50mMLRi);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
title('RMSE: M3 -- MLRi')
% clim([quantile(RMSE50mMLRi(:),[0.05,0.95])])
clim([0 50])
colorbar

%ME
nexttile
h = imagesc(Coords.X./1000,Coords.Y./1000,ME50mMLRi);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
title('ME: M3 -- MLRi')
% clim([quantile(ME50mMLRi(:),[0.05,0.95])])
clim([-50,50])
colorbar

% M3 vs. GPR
%R
nexttile
h = imagesc(Coords.X./1000,Coords.Y./1000,R50mGPR);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
clim([-1 1])
colorbar
title('R: M3 -- GP')
%RMSE
nexttile
h= imagesc(Coords.X./1000,Coords.Y./1000,RMSE50mGPR);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
% clim([quantile(RMSE50mGPR(:),[0.05,0.95])])
clim([0 50])
colorbar
title('RMSE: M3 -- GP')

%ME
nexttile
h = imagesc(Coords.X./1000,Coords.Y./1000,ME50mGPR);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);
title("ME: M3 -- GP")
% clim([[quantile(ME50mGPR(:),[0.05,0.95])]])
clim([-50 50])
colorbar


% MLRi Vs. GPR
%R
nexttile
h=imagesc(Coords.X./1000,Coords.Y./1000,R5m);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph2)
set(gca,'color',0.75*[1 1 1]);
% clim([quantile(R5m(:),[0.05,0.95])])
clim([-1 1])
colorbar
title('R: MLRi -- GP')
%RMSE
nexttile
h = imagesc(Coords.X./1000,Coords.Y./1000,RMSE5m);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph2)
set(gca,'color',0.75*[1 1 1]);
% clim([quantile(RMSE5m(:),[0.05,0.95])])
clim([0 50])
colorbar
title('RMSD: MLRi -- GP')
%ME
nexttile
h = imagesc(Coords.X./1000,Coords.Y./1000,ME5m);daspect([1,1,1]);colormap(cmap)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',8)
set(h,'AlphaData',~alph2)
set(gca,'color',0.75*[1 1 1]);
% clim([quantile(ME5m(:),[0.05,0.95])])
clim([-50 50])
colorbar
title('ME: MLRi -- GP')