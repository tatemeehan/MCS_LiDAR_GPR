function [A] = gapfillMCS(A,deg)
% gapfillMCS.m infills nan values using inpaint_nans.m. This operation
% takes ages on the rasters which have ~43% NaN values due to edge padding.
% To remedy this, the raster is rotated to minimize NaN padding, then
% gapfilled, and lastly untransformed back to the original georefferenced
% orientation.

if nargin < 2
    deg = 36.5;% Default Rotation
end
tmpA = A;
%% Rotate Rasters
% Rotate C 
A = imrotate(A,-deg);
A(A<=0) = NaN;
% Trim NaNs
ix = find(~isnan(A));
[row,col] = ind2sub(size(A),ix);
% Bounding Box
minrow = min(row); maxrow = max(row);
mincol = min(col); maxcol = max(col);
A = A(minrow:maxrow,mincol:maxcol);
% NaN Padding Mask
nanMaskR = double(~isnan(A));
nanMaskR(nanMaskR ~=1) = NaN;
% Inpaint NaNs
A = inpaint_nans(A,5);
A = A.*nanMaskR;
%% Inverse Transform
A = imrotate(A,deg);
A(A<=0) = NaN;
% Trim NaNs
ix = find(~isnan(A));
[row,col] = ind2sub(size(A),ix);
% Bounding Box
minrow = min(row); maxrow = max(row);
mincol = min(col); maxcol = max(col);
A = A(minrow:maxrow,mincol:maxcol);
A = imresize(A,size(tmpA));
end