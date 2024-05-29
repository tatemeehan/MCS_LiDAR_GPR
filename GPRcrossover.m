%% GPR Cross Over Loctions
cmap = csvread('D:\git-repository\GreenTrACS_MxRadar\colorMaps\RdYlBu.csv');
cmap = flipud(cmap);
% Load Data
% GPR Data
GPR = readtable(['E:\MCS\MCS021324\GPR\processed\GPR-TWT.csv']);
gprX = GPR.Easting; gprY = GPR.Northing;
gprTWT = GPR.TWT; gprZ = GPR.ElevationWGS84;
% KDtree
load('E:\MCS\MCS021324\GPR\processed\20240213_MCS-kdtree.mat')
ix = kd.ix;
D = kd.D;
IDX = kd.IDX;
% Snow Depth Raster
dataDir = 'E:\MCS\MCS021324\LiDAR';
filename = '20240213_MCS-snowdepth_RFgapfilled.tif';
outfnA = filename(1:end-4);
fullfilename = fullfile(dataDir,filename);
[A,RA,~,~,lon,lat,utmX,utmY] = readLidarTif(fullfilename);
% Get UTM Coordinates as Vector
X = utmX(1,:);
Y = utmY(:,1);
Xi = utmX(:); 
Yi = utmY(:);

isGPRcrossOver = 1;
if isGPRcrossOver
% Get GPR into UnixTime
gprUnixTime = posixtime(GPR.DateTime); gprUnixTime = inpaint_nans(gprUnixTime);
timeThresh = 10;
iter = 1;
    for jj = 1:length(ix)
        % Calculate Distance
        dist = D{jj};dist = dist(:);
        % Find Points within radius r of grid point
        tmpIx = IDX{jj};tmpIx = tmpIx(:);
        % Find Points that CrossOver
        tmpTime = gprUnixTime(tmpIx);
        diffTime = abs(tmpTime-tmpTime(1));
        bin1ix = find(diffTime<timeThresh);
        bin2ix = find(diffTime>timeThresh);
        if ~isempty(bin2ix) %& ~ismember(jj,threshIx)
            % Variance
            GPRtwtvar(iter,1) = var(gprTWT(tmpIx));
            GPRtwtvar1(iter,1) = var(gprTWT(tmpIx(bin1ix)));
            GPRtwtvar2(iter,1) = var(gprTWT(tmpIx(bin2ix)));

        % Median
        GPRtwt1(iter,1) = median(gprTWT(tmpIx(bin1ix)));
        GPRtwt2(iter,1) = median(gprTWT(tmpIx(bin2ix)));
        % Coordinates
        GPRxcross(iter,1) = Xi(ix(jj));
        GPRycross(iter,1) = Yi(ix(jj));
        iter = iter + 1;
        end
    end
    % Repeat this process to Calculate the median of Nodes
    binRadius = 5; % 10 meter radius
    for kk = 1:length(GPRxcross)
        allDist(:,kk) = sqrt((GPRxcross(kk)-GPRxcross).^2+(GPRycross(kk)-GPRycross).^2);
        isInBin(:,kk) = allDist(:,kk)<binRadius;
    end
    isInBin = tril(isInBin);
    crossover = 1;
    iter = 1;
    ixBin = [];
    while crossover
        % Take the median within 10 meters
        ix1 = find(isInBin(:,iter));
        rowBin = max(ix1);
        ix2 = find(isInBin(rowBin,:));
        ixBin = [ixBin;ix1];
        
        nodeGPRtwt1(iter) = median(GPRtwt1(ix1));
        nodeGPRtwt2(iter) = median(GPRtwt2(ix1));
        nodeGPRxcross(iter) = median(GPRxcross(ix1));
        nodeGPRycross(iter) = median(GPRycross(ix1));
        iter = max(ixBin)+1;
        if iter > length(GPRxcross)
            crossover = 0;
        end
    end
    nodeGPRxcross(nodeGPRxcross == 0) = [];
    nodeGPRycross(nodeGPRycross == 0) = [];
    nodeGPRtwt1(nodeGPRtwt1==0) = [];
    nodeGPRtwt2(nodeGPRtwt2==0) = [];
    
    % Determine Nominal GPR Grid Spacing
    GPRdist = zeros(length(nodeGPRxcross),1);

    for kk = 1:length(nodeGPRxcross)
        neighborIx = 1:length(nodeGPRxcross);
        neighborIx(kk) = [];
        tmpGPRdist = (sqrt((nodeGPRxcross(kk) - nodeGPRxcross(neighborIx)).^2+(nodeGPRycross(kk) - nodeGPRycross(neighborIx)).^2));
        % Threshold Neighbors between 5 & 100 m (not nessarry but oh well) 
        tmpGPRdistIx = find(tmpGPRdist>5 & tmpGPRdist < 100);
        if isempty(tmpGPRdistIx)
            GPRdist(kk) = NaN;
        else
            GPRdist(kk) = min(tmpGPRdist(tmpGPRdistIx));
        end
    end
    % Summary Statistics
    avgGridH = nanmean(GPRdist);
    [Rtwtdiff,Ptwtdiff] = corrcoef(GPRtwt1,GPRtwt2);
    meanTWTdiff = round(mean(GPRtwt1-GPRtwt2),1);
    stdTWTdiff = round(std(GPRtwt1-GPRtwt2),1);
    meannodeTWTdiff = round(mean(nodeGPRtwt1-nodeGPRtwt2),1);
    stdnodeTWTdiff = round(std(nodeGPRtwt1-nodeGPRtwt2),1);
    
    % Summary Figures
    % Plot the Corss-over Locations
    figure();plot(GPRxcross,GPRycross,'.k')
    title('GPR Cross-over Locations')
    xlabel('Easting (m)')
    ylabel('Northing (m)')
    grid on
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif')
    annotation('textbox',[0.525,0.5525,.2,.2],'string',['Average Grid Spacing = ',num2str(round(avgGridH)),' m'],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',12)

    % Corrleation Plot of TWT differences
    figure();
    binscatter(GPRtwt1,GPRtwt2,100);colorbar off;
    title('Correlation of Two-way Travel-time Difference at Cross-over Locations')
    xlabel('TWT (ns)'); ylabel('TWT (ns)');
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif')
    grid on
    annotation('textbox',[0.625,0.2325,.2,.2],'string',['R = ',num2str(round(Rtwtdiff(1,2),2)),' RMSE = ', num2str(stdTWTdiff)],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',12)
    

    % Histogram of TWT differences
    smoothTWTdiff = ksdensity(GPRtwt1-GPRtwt2,-4:0.05:4);
    figure();
    histogram(GPRtwt1-GPRtwt2,'facecolor',[0.5,0.5,0.5],'edgecolor','k','normalization','pdf','facealpha',0.5,'NumBins',35);
    hold on; plot(-4:0.05:4,smoothTWTdiff,'k','linewidth',2)
    set(gca,'fontsize',12,'fontweight','bold','fontname','serif')
    title('Two-Way Travel-time Difference at Cross-over Locations')
    ylabel('PDF')
    xlabel('TWT Difference (ns)')
    xlim([-4,4])
    grid on
    annotation('textbox',[0.625,0.24,.2,.2],'string',['\mu = ',num2str(meanTWTdiff),' \sigma = ',num2str(stdTWTdiff)],'fitboxtotext','on','linestyle','none','fontname','serif','fontsize',12)
end