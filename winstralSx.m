function [Sx,Sb,SxSep] = winstralSx(X, Y, Z, windDir, dmax, dsep)
%winstralSx computes the maximum upwind slope parameter Sx 
%   Input: X & Y - MeshGrid Style Coordinates of Raster. If vectors are
%          supplied the Coordinates will be gridded. UTM Coordinates!
%          Z - Digital Elevation Model Raster of Surface of interest
%          windDir - Either a Scalar or Time Serires of Desired Wind
%          Direciton For upwind analysis. If input is a meterorlogical time
%          series, this code will determine the median direction of snow
%          transportable winds following the transport threshold of Li and
%          Pomeroy (1997)
%          dmax - The Search Radius distance (Default 100 m)
%          dsep - The Longer Length Scale for SxSep (Default 1000 m)
%
%   Output: Sx - The maximum upwind slope raster
%           Sb - The Slope Break Parameter (Sx - SxSep)
%           SxSep - The long range upwind slope parameter

% Check Inputs
[m,n] = size(Z);
[mx,nx] = size(X); 
[my,ny] = size(Y);
if (mx == 1 || nx == 1) && (my == 1 || ny == 1)
    X = X(:); Y = Y(:);
%     if length(X) ~= length(Y)
%         error('X and Y must be  Vectors of the same size!')
%     end
    [X,Y] = meshgrid(X,Y);
elseif (mx > 1 && nx > 1) && (my > 1 && ny > 1)
    if mx ~= my || nx ~= ny
        error('X and Y must be Matricies of the same size!')
    end
else
    error('X and Y must be Either Vectors or Matricies of the same size!')
end

if numel(Z) ~= numel(X)
    error('Z must be of compatible size as X & Y coordinates!')
end
% Check Wind Direction
if isscalar(windDir)
    % Change from North Orientation 0 degrees to Trigonometric Orientation
    WindSpeed = ones(numel(windDir),1);
    [u,v] = pol2cart(deg2rad(windDir),WindSpeed);
    [windDir,~]=cart2pol(v,u);
    % Convert to Degrees
    windDir = rad2deg(windDir);
    % Convert to 0 - 360 with 0 facing East using Meterological Convention
    windDir = mod(windDir + 360,360);
else
    threshold = 5; % Approximately Lower Theshold from Li & Pomeroy (1997).
    threshold = threshold.*(ones(numel(windDir),1));
    wIx = find(windDir > threshold);
    windDir = median(windDir(wIx));
    % Change from North Orientation 0 degrees to Trigonometric Orientation
    WindSpeed = ones(numel(windDir),1);
    [u,v] = pol2cart(deg2rad(windDir),WindSpeed);
    [windDir,~]=cart2pol(v,u);
    % Convert to Degrees
    windDir = rad2deg(windDir);
    % Convert to 0 - 360 with 0 facing East using Meterological Convention
    windDir = mod(windDir + 360,360);
end

% Check dmax
if exist("dmax")
else
    % Set Default at 100 m
    dmax = 100;
end
if exist("dsep")
else
    % Set default at 1000 m
    dsep = 1000;
end

% Calculate Distance to All Points
dx = abs(floor(X(1,2) - X(1,1)));
Zsep = Z; Xsep = X; Ysep = Y;
if dx < 10
    if dmax <= 50
        % Resize
        X = imresize(X,dx/5);
        Y = imresize(Y,dx/5);
        Z = imresize(Z,dx/5);
    else
        % Resize
        X = imresize(X,dx/10);
        Y = imresize(Y,dx/10);
        Z = imresize(Z,dx/10);
    end
    if dsep<=250
        Xsep = imresize(Xsep,dx/10);
        Ysep = imresize(Ysep,dx/10);
        Zsep = imresize(Zsep,dx/10);
    else
        Xsep = imresize(Xsep,dx/30);
        Ysep = imresize(Ysep,dx/30);
        Zsep = imresize(Zsep,dx/30);
    end
    [dIx,D] = rangesearch([X(:),Y(:)],[X(:),Y(:)],dmax);
    [dsepIx,Dsep] = rangesearch([Xsep(:),Ysep(:)],[Xsep(:),Ysep(:)],dsep);
else
    % Use Native Resolution
    [dIx,D] = rangesearch([X(:),Y(:)],[X(:),Y(:)],dmax);
    if dx < 30
        % Check Outlier Scale Break Resolution
        Xsep = imresize(Xsep,dx/30);
        Ysep = imresize(Ysep,dx/30);
        Zsep = imresize(Zsep,dx/30);
        [dsepIx,Dsep] = rangesearch([Xsep(:),Ysep(:)],[Xsep(:),Ysep(:)],dsep);
    else
        % Use Native Resolution
        [dsepIx,Dsep] = rangesearch([Xsep(:),Ysep(:)],[Xsep(:),Ysep(:)],dsep);
    end
end
% Loop For Upwind Maximum Slope
ix = 1:numel(Z);
Sx = zeros(size(Z)); SxSep = zeros(size(Zsep));
for ii = 1:numel(Z)
    x0 = X(ii); y0 = Y(ii); ix0 = ix(ii);
    % Theta for Sx
    [theta,~] = cart2pol(x0-X(dIx{ix0})',y0-Y(dIx{ix0})');
    theta = mod(rad2deg(theta)+360,360);
    % Wind Direciton Buffers
     windBin = - 15 : 5 : 10;
    for kk = 1:length(windBin)
         binIx = find(theta>=windDir+windBin(kk) & theta<windDir+(windBin(kk)+5));
         % Calculate Sx
         if isempty(binIx)
             tmpSx(kk) = nan;
         elseif (atand((Z(dIx{ix0}(binIx))-Z(ix0))./...
                 sqrt((X(dIx{ix0}(binIx))-X(ix0)).^2+(Y(dIx{ix0}(binIx))-Y(ix0)).^2)) < 0)
             tmpSx(kk) = min(atand((Z(dIx{ix0}(binIx))-Z(ix0))./...
                 sqrt((X(dIx{ix0}(binIx))-X(ix0)).^2+(Y(dIx{ix0}(binIx))-Y(ix0)).^2)));
         else
             tmpSx(kk) = max(atand((Z(dIx{ix0}(binIx))-Z(ix0))./...
                 sqrt((X(dIx{ix0}(binIx))-X(ix0)).^2+(Y(dIx{ix0}(binIx))-Y(ix0)).^2)));
         end
    end
         % Average the Wind Bins
         Sx(ix0) = nanmean(tmpSx);
end
% Wind Break Loop
ix = 1:numel(Zsep);
for ii = 1:numel(Zsep)
    x0 = Xsep(ii); y0 = Ysep(ii); ix0 = ix(ii);
    % Use dmax as dsep cutoff
%     rmvIx = find(Dsep{ix0}<dmax);
    tmpDsepIx = dsepIx{ix0};
    tmpDsepIx(Dsep{ix0}<dmax) = [];
    % Theta SxSep
%     [thetaSep,~] = cart2pol(x0-Xsep(dsepIx{ix0})',y0-Ysep(dsepIx{ix0})');
    [thetaSep,~] = cart2pol(x0-Xsep(tmpDsepIx)',y0-Ysep(tmpDsepIx)');
    thetaSep =  mod(rad2deg(thetaSep)+360,360);
    % Wind Direciton Buffers
    windBin = - 15 : 5 : 10;
    for kk = 1:length(windBin)
        % Calculate SxSep
        binIxSep = find(thetaSep>=windDir+windBin(kk) & thetaSep<windDir+(windBin(kk)+5));
        if isempty(binIxSep)
            tmpSxSep(kk) = nan;
        elseif atand((Zsep(tmpDsepIx(binIxSep))-Zsep(ix0))./...
                sqrt((Xsep(tmpDsepIx(binIxSep))-Xsep(ix0)).^2+(Ysep(tmpDsepIx(binIxSep))-Ysep(ix0)).^2)) < 0
            tmpSxSep(kk) = min(atand((Zsep(tmpDsepIx(binIxSep))-Zsep(ix0))./...
                sqrt((Xsep(tmpDsepIx(binIxSep))-Xsep(ix0)).^2+(Ysep(tmpDsepIx(binIxSep))-Ysep(ix0)).^2)));
        else
            tmpSxSep(kk) = max(atand((Zsep(tmpDsepIx(binIxSep))-Zsep(ix0))./...
                sqrt((Xsep(tmpDsepIx(binIxSep))-Xsep(ix0)).^2+(Ysep(tmpDsepIx(binIxSep))-Ysep(ix0)).^2)));
        end
%         elseif atand((Zsep(dsepIx{ix0}(binIxSep))-Zsep(ix0))./...
%                 sqrt((Xsep(dsepIx{ix0}(binIxSep))-Xsep(ix0)).^2+(Ysep(dsepIx{ix0}(binIxSep))-Ysep(ix0)).^2)) < 0
%             tmpSxSep(kk) = min(atand((Zsep(dsepIx{ix0}(binIxSep))-Zsep(ix0))./...
%                 sqrt((Xsep(dsepIx{ix0}(binIxSep))-Xsep(ix0)).^2+(Ysep(dsepIx{ix0}(binIxSep))-Ysep(ix0)).^2)));
%         else
%             tmpSxSep(kk) = max(atand((Zsep(dsepIx{ix0}(binIxSep))-Zsep(ix0))./...
%                 sqrt((Xsep(dsepIx{ix0}(binIxSep))-Xsep(ix0)).^2+(Ysep(dsepIx{ix0}(binIxSep))-Ysep(ix0)).^2)));
%         end
    end
    % Average the Wind Bins
    SxSep(ix0) = nanmean(tmpSxSep);
end
% Inpaint NaN
Sx = inpaint_nans(Sx,5);
SxSep = inpaint_nans(SxSep,5);
% Resize Sx
Sx = imresize(Sx,[m,n]);
% Resize SxSep
SxSep = imresize(SxSep,[m,n]);
% Calculate Slope Break Parameter
Sb = Sx - SxSep;
