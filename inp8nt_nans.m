function [A] = inp8nt_nans(A,mask,n)
% inp8nt_nans.m uses a local distance weighting to gap fill nan values

% input: A - 2D array
%        mask - 2D array of known nan values chosen to exclude from gapfill

% output: A - Gap Filled Array

if nargin < 2
    mask = ones(size(A));
elseif nargin < 3
    n = 10; % default 21x21 pixel subgrid box
end
[M,N] = size(A);
% Describe Known and Unknown Indicies
% ix = 1:numel(A);
nanix = find(isnan(A));
% ix(nanix) = [];
% [ixY,ixX] = ind2sub(size(A),ix);

maskix = find(mask);
nanix = setdiff(nanix,maskix);
% nanix(rmvix) = [];
[nanY,nanX] = ind2sub(size(A),nanix);
% n = 25; % number of nearest neighbors to average


% mykdtree=KDTreeSearcher([nanX(:) nanY(:)]); % searcher for NaN Values
% [IDX,D]=knnsearch(mykdtree,[ixX(:), ixY(:)],'k',10); % search the points in the LiDAR grid

% mykdtree=KDTreeSearcher([ixX(:), ixY(:)]); % searcher for NaN Values
% [IDX,D]=knnsearch(mykdtree,[nanX(:) nanY(:)],'k',n); % search the points in the LiDAR grid

% % Distance Weighting
% for kk = 1:numel(nanix)
%     dist = sqrt((nanX(kk)-ixX).^2+((nanY(kk)-ixY).^2));
%     [dists,distIx] = sort(dist);
% 
%      % calculate weights - bisquare kernel
%         w=15/16*(1-((dists(1:n)./n)+eps).^2).^2;
%         w=w(:);
%     A(nanix(kk)) = sum(w'.*A(ix(distIx(1:n))))./sum(w);
% end

% SubGrid Method
for kk = 1:numel(nanix)
    % Determine Y coordinates
    if nanY(kk) - n < 1
        by1 = 1; by2 = nanY+n;
    elseif nanY(kk) + n > M
        by1 = nanY(kk) - n; by2 = M;
    else
        by1 = nanY(kk) - n; by2 = nanY(kk)+n;
    end
    % Determine X Coordinates
    if nanX(kk) - n < 1
        bx1 = 1; bx2 = nanX+n;
    elseif nanX(kk) + n > N
        bx1 = nanX(kk) - n; bx2 = N;
    else
        bx1 = nanX(kk) - n; bx2 = nanX(kk)+n;
    end
    % Sort Coordinates and Remove neighboring NaNs
    bY = by1:by2; bX = bx1:bx2;
    [BX,BY] = meshgrid(bX,bY);
    BY = BY(:);BX = BX(:);
    box = A(bY,bX); box = box(:);
    rmvix = find(isnan(box));
    box(rmvix) = []; BY(rmvix) = []; BX(rmvix) = [];

    % Distance Weighting
    dist = sqrt((nanX(kk)-BX).^2+((nanY(kk)-BY).^2));

    % calculate weights - bisquare kernel
    w=15/16*(1-((dist./max(dist))+eps).^2).^2;
    w=w(:);
    tmp = sum(w.*box)./sum(w);
    if isnan(tmp)
        keyboard; % Need Expanding Window
    else
        A(nanix(kk)) = tmp;
    end
end
end