%% Make A GIF!
clear; close all; clc;
%% Read LiDAR Data
dataDir = 'E:\MCS\GIF';
% Dec. Snow Depth
filename = 'MCS_20231228_SNOWDEPTH_RFgapfilled.tif';
fullfilename = fullfile(dataDir,filename);
[A,R,~,~,lon,lat,utmX,utmY] = readLidarTif(fullfilename);
% Resize Raster
A = imresize(A,.1,"nearest");
% ResizeCoords
lat = imresize(lat,.1,"nearest");
lon = imresize(lon,.1,"nearest");
utmX = imresize(utmX,.1,"nearest");
utmY = imresize(utmY,.1,"nearest");
Y = utmY(:,1);
X = utmX(1,:);




%% Read LiDAR Data
dataDir = 'E:\MCS\GIF';
% January Snow Depth
filename = 'MCS_20240115_SNOWDEPTH_RFgapfilled.tif';
fullfilename = fullfile(dataDir,filename);
[B,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
% Resize Raster
B = imresize(B,.1,"nearest");

%% Read LiDAR Data
dataDir = 'E:\MCS\GIF';
% Feb. Snow Depth
filename = '20240213_MCS-snowdepth_RFgapfilled.tif';
fullfilename = fullfile(dataDir,filename);
[C,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
% Resize Raster
C = imresize(C,.1,"nearest");

%% Msrch Snow Depth
filename = '20240315_MCS-snowdepth_RFgapfilled.tif';
fullfilename = fullfile(dataDir,filename);
[D,~,~,~,~,~,~,~] = readLidarTif(fullfilename);
% Resize Raster
D = imresize(D,.1,"nearest");
%% Image Stack
% Stack = cat(3,B,C,D,A);
Stack = cat(3,A,B,C,D);

%% Make The Gif
cmap = csvread('..\codeRepo\colormaps\RdYlBu.csv');
cmap = flipud(cmap);
alph=isnan(A);
figure();
for ii = 1:40
hI = imagesc(X,Y,Stack(:,:,ceil(ii./10)));
% h = colorbar; %ylabel(h,'Snow Depth (m)','FontWeight','bold','FontName','serif','FontSize',14)
daspect([1,1,1]);colormap(cmap);caxis([0.5,2.5])%caxis(round([quantile(Stack(:),[0.1,0.99])./0.05]).*0.05);%caxis([0,5]);
axis off
% title({'Mores Creek Summit: 04/05/2023','Hanover Reprocessed Snow Depth'})
% title({strjoin(["Mores Creek Summit: ",num2str(outfnA(9:10)),"/",num2str(outfnA(11:12)),"/",num2str(outfnA(5:8))]),'LiDAR Snow Depth'})
% xlabel('Easting (km)');ylabel('Northing (km)');
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
set(hI,'AlphaData',~alph)
set(gca,'color','none');
pause(0.1)
exportgraphics(gcf,[dataDir,'\','gifford.gif'],'Append',true);
end


