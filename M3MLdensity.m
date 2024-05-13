%% Compare ML Density to M3 Density
clear; close all; clc;
%% Load ML & M3 data
% M3 Density Model
M3dataDir = 'E:\MCS\MCS040623\M3';
M3filename = 'snow_density_2023-04-06T230000.tif';
M3fullfilename = fullfile(M3dataDir,M3filename);
[M3,M3R,~,~,M3lon,M3lat,M3utmX,M3utmY] = readLidarTif(M3fullfilename);
% Interpolate to M3 Grid
xq = M3R.XWorldLimits(1):M3R.CellExtentInWorldX:M3R.XWorldLimits(2);
yq = M3R.YWorldLimits(1):M3R.CellExtentInWorldY:M3R.YWorldLimits(2);
[Xq,Yq] = meshgrid(xq,yq);
Yq = flipud(Yq);

% Load Coords
LiDARdataDir = 'E:\MCS\MCS040623\LiDAR';
filename = 'MCS040523-Coords';
fullfilename = fullfile(LiDARdataDir,filename);
load(fullfilename);

% ML MLR Density Model
MLdataDir = 'E:\MCS\MCS040623\GPR';
MLfilename = 'MCS040623-MLresultsMLR.mat';
MLfullfilename = fullfile(MLdataDir,MLfilename);
load(MLfullfilename);

% ML GPR Density Model
MLdataDir = 'E:\MCS\MCS040623\GPR';
MLfilename = 'MCS040623-MLresultsGPR.mat';
MLfullfilename = fullfile(MLdataDir,MLfilename);
load(MLfullfilename);


%% Upscale to 50 m
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

%% Make Some figures
cmap = csvread('..\codeRepo\colormaps\RdYlBu.csv');
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
hM3 = imagesc(x./1000,flipud(y)./1000,M3d);colorbar
daspect([1,1,1]);
colormap(cmap)
title('M3')
caxis([250,450])
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
set(hM3,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

subplot(1,3,2)
hMLRi = imagesc(x./1000,flipud(y)./1000,MLRi);colorbar
daspect([1,1,1]);
colormap(cmap)
title('MLRi')
caxis([250,450])
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
set(hMLRi,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);

subplot(1,3,3)
hGPR = imagesc(x./1000,flipud(y)./1000,GPR);colorbar
daspect([1,1,1]);
colormap(cmap)
title('GPR')
caxis([250,450])
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
set(hGPR,'AlphaData',~alph)
set(gca,'color',0.75*[1 1 1]);