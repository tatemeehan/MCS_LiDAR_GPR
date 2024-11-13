%% Compute SWE
clear; close all; clc;
addpath(genpath('E:\MCS\MCS_LiDAR_GPR\'))
cmap = csvread('RdYlBu.csv');
cmap = flipud(cmap);

%
% dataDir = ['E:\MCS\MCS031524\'];
dataDir = ['E:\MCS\MCS031524\MCS03152024_5m\'];
writeDir = dataDir;
% Load Density .Tif
fn = '20240315_MCS_GPRdensity_5m.tif';
% fn = '20240315_MCS_GPRdensity.tif';
[Density,R,~,~,lon,lat,utmX,utmY,epsgCode] = readLidarTif([dataDir,fn]);
X = utmX(1,:);
Y = utmY(:,1);
isLoadMat = 1;
isLoadTif = 0;
if isLoadMat
% Load LiDAR.mat
    % Coordinate Data
    load([dataDir,'20240315_MCS-Coords5m.mat'])
    % LiDAR Data
    load([dataDir,'20240315_MCS-LiDAR5m.mat'])
    % Extract pertainent LiDAR data
    A = LiDAR.A.A;
    % HillShade
    Caspect = LiDAR.C.aspect;
%     Cnorthness = LiDAR.C.northness;
%     Ceastness = LiDAR.C.eastness;
    % SnowDepth
    % depth = LiDAR.A.RF;
    depth = LiDAR.A.A;
%     Mask
    alph=~isnan(A);
%     out of memory.
%     clear LiDAR
end
if isLoadTif
    % Load LiDAR.tif
[A,~,~,~,~,~,~,~,~] = readLidarTif([dataDir,'20240315_MCS-snowdepth_RFgapfilled.tif']);
% C Reference DEM
dataDir2 = 'E:\MCS\';
% filename = '20220317_MCS-refDEM-T8.tif';
% Create Georeference
latlim = [min(lat(:)),max(lat(:))];
lonlim = [min(lon(:)) max(lon(:))];
sizeLidar = size(A);
georef = georefpostings(latlim,lonlim,sizeLidar,'RowsStartFrom','west','ColumnsStartFrom','north');
filename = 'MCS_REFDEM_WGS84.tif';
[C,~,~,~,~,~,~,~,~] = readLidarTif([dataDir2,filename]);
[Caspect,~,~,~] = gradientm(C,georef);
alph=~isnan(A);

end
% SWE
SWE = A.*Density;

%% Write Geotiff
geotiffwrite([writeDir,'20240315_MCS','_GPRswe_5m.tif'],SWE,R,'CoordRefSysCode',epsgCode);
geotiffwrite([writeDir,'20240315_MCS','_snowdepth_RFgapfilled_5m.tif'],A,R,'CoordRefSysCode',epsgCode);


%% Make Figures
   FigDepth = figure();
    % Hillshade
        imagesc((X)./1000,(Y)./1000,(cosd(Caspect+45)+sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
%     imagesc((Coords.X)./1000,(Coords.Y)./1000,(-cosd(Caspect+45)-sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
    % imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
    freezeColors; hold on;
    % MLR Density
    hI = imagesc((X)./1000,(Y)./1000,A,'AlphaData',0.625.*alph);
    daspect([1,1,1]);set(gca,'YDir','normal');
    colormap(cmap);
    caxis([quantile(A(alph),[0.005,0.995])]);
    clim([0.5 3])
    cb = colorbar;cb.FontSize = 14;
    cb.Label.String = 'Snow Depth (m)';
    title(['Mores Creek Summit: 03/15/2024']);%, dataDir(end-5:end)])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
    xtickformat('%.0f')
    ax.YAxis.Exponent = 0;
    ytickformat('%.0f')
    xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
    % Save Figure
    saveas(FigDepth, [writeDir,'\MCS031524_snowDepth_5m.png'],'png');
% SWE figre
       FigSWE = figure();%'units','normalized','outerposition',[0 0 1 1]);
    % Hillshade
        imagesc((X)./1000,(Y)./1000,(cosd(Caspect+45)+sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
%     imagesc((Coords.X)./1000,(Coords.Y)./1000,(-cosd(Caspect+45)-sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
    % imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
    freezeColors; hold on;
    % MLR Density
    hI = imagesc((X)./1000,(Y)./1000,SWE,'AlphaData',0.625.*alph);
    daspect([1,1,1]);set(gca,'YDir','normal');
    colormap(cmap);
%     caxis([quantile(SWE(alph),[0.005,0.995])]);
    clim([100 1000])
    cb = colorbar;cb.FontSize = 14;
    cb.Label.String = 'Snow Water Equivalent (mm)';
    title(['Mores Creek Summit: 03/15/2024']);%, dataDir(end-5:end)])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
    xtickformat('%.0f')
    ax.YAxis.Exponent = 0;
    ytickformat('%.0f')
    xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
    % Save Figure
    saveas(FigSWE, [writeDir,'\MCS031524_swe_5m.png'],'png');


% Density figre
       FigDensity = figure();%'units','normalized','outerposition',[0 0 1 1]);
    % Hillshade
        imagesc((X)./1000,(Y)./1000,(cosd(Caspect+45)+sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
%     imagesc((Coords.X)./1000,(Coords.Y)./1000,(-cosd(Caspect+45)-sind(Caspect+45))./2,'AlphaData',0.0625.*~alph+alph);colormap(bone)
    % imagesc((Coords.X)./1000,(Coords.Y)./1000,(-LiDAR.C.northness-LiDAR.C.eastness)./2,'AlphaData',0.25.*~alph+alph);colormap(bone)
    freezeColors; hold on;
    % MLR Density
    hI = imagesc((X)./1000,(Y)./1000,Density,'AlphaData',0.625.*alph);
    daspect([1,1,1]);set(gca,'YDir','normal');
    colormap(cmap);
%     caxis([quantile(SWE(alph),[0.005,0.995])]);
    clim([320 360])
    cb = colorbar;cb.FontSize = 14;
    cb.Label.String = 'Snow Density (kg/m^3)';
    title(['Mores Creek Summit: 03/15/2024']);%, dataDir(end-5:end)])
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif');
    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
    xtickformat('%.0f')
    ax.YAxis.Exponent = 0;
    ytickformat('%.0f')
    xlabel('Easting (km)','fontsize',14); ylabel('Northing (km)','fontsize',14);
    % Save Figure
    saveas(FigDensity, [writeDir,'\MCS031524_density_5m.png'],'png');