function [R,RMSE,ME] = movCorr2D(A,B,r)
%2D Moving Window Correlation
%   Computes Pearson Correlation Coefficient between two matricies using a
%   square window of user defined size.
% Paralell Processing Enabled
p = gcp('nocreate');
if isempty(p)
    % There is no parallel pool
    poolsize = 0;
else
    % There is a parallel pool of <p.NumWorkers> workers
    poolsize = p.NumWorkers;
end

if nargin<3
    r = 5;
end
[m,n] = size(A);
[mb,nb] = size(B);
if m~= mb || n ~= nb
    warning('Matricies are not the same size.. They will be resized to match.')
    B = imresize(B,[m,n],Method="bilinear");
end
R = zeros(m,n);RMSE = R;ME = RMSE;
r2 = floor(r./2);
if poolsize == 0
for ii = 1:m
    if ii <= r2
        mIx = 1:ii+r2;
    elseif m-ii <= r2
        mIx = ii-r2:m;
    else
        mIx = ii-r2:ii+r2;
    end
    for jj = 1:n
        if jj <= r2
            nIx = 1:jj+r2;
        elseif n-jj <= r2
            nIx = jj-r2:n;
        else
            nIx = jj-r2:jj+r2;
        end
        % Create 1D Indicies of Windowed Data
        IX = combvec(mIx,nIx);
        ix = sub2ind([m,n],IX(1,:),IX(2,:));
        % Nan Filter
        nanIx = (isnan(A(ix)) | isnan(B(ix)));
        ix(nanIx) = [];
        if isempty(ix) || numel(ix) < r2+1
            R(ii,jj) = NaN;
            RMSE(ii,jj) = NaN;
            ME(ii,jj) = NaN;
        else
            tmp =  corrcoef(A(ix),B(ix));
            R(ii,jj) = tmp(2,1);
            error = B(ix)-A(ix);
            RMSE(ii,jj) = sqrt(mean((error).^2));
            ME(ii,jj) = mean(error);
        end
    end
end
else
    % Paralell Processing
    parfor (ii = 1:m,poolsize)
        if ii <= r2
            mIx = 1:ii+r2;
        elseif m-ii <= r2
            mIx = ii-r2:m;
        else
            mIx = ii-r2:ii+r2;
        end
        for jj = 1:n
            if jj <= r2
                nIx = 1:jj+r2;
            elseif n-jj <= r2
                nIx = jj-r2:n;
            else
                nIx = jj-r2:jj+r2;
            end
            % Create 1D Indicies of Windowed Data
            IX = combvec(mIx,nIx);
            ix = sub2ind([m,n],IX(1,:),IX(2,:));
            % Nan Filter
            nanIx = (isnan(A(ix)) | isnan(B(ix)));
            ix(nanIx) = [];
            if isempty(ix) || numel(ix) < r2+1
                R(ii,jj) = NaN;
                RMSE(ii,jj) = NaN;
                ME(ii,jj) = NaN;
            else
                tmp =  corrcoef(A(ix),B(ix));
                R(ii,jj) = tmp(2,1);
                error = B(ix)-A(ix);
                RMSE(ii,jj) = sqrt(mean((error).^2));
                ME(ii,jj) = mean(error);
            end
        end
    end
end
end