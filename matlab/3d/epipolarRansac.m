function varargout = epipolarRansac(pts1, pts2, varargin)
% EPIPOLARRANSAC Applies ransac to estimate stereo fundamental matrix.
%
% Other m-files required: epipolarGeometry.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

defaultMaxiter              = 2000;
defaultErrorThreshold       = 0.05;
defaultAlsoInliersThreshold = 20;

validPts                  = @(x) isnumeric(x) && (size(x,2) == 2);
validMaxiter              = @(x) isscalar(x)  && (round(x) == x);
validErrorThreshold       = @(x) isscalar(x)  && (x > 0);
validAlsoInliersThreshold = @(x) isscalar(x)  && (round(x) == x);

p = inputParser;
addRequired(p,'pts1', validPts);
addRequired(p,'pts2', validPts);
addOptional(p,'maxiter', defaultMaxiter, validMaxiter);
addOptional(p,'errorThreshold', defaultErrorThreshold, validErrorThreshold);
addOptional(p,'alsoInliersThreshold', defaultAlsoInliersThreshold, validAlsoInliersThreshold);

parse(p,pts1,pts2,varargin{:});

npts = size(pts1, 1);

bestF     = zeros(3);
bestError = inf;
inliers   = [];

fprintf("Random consensus to find epipolar geometry with %d matches\n", npts);

for iter=1:p.Results.maxiter
    % Select random inlier indices (minimum is 9)
    maybeInliers = randperm(npts, p.Results.alsoInliersThreshold);

    % Match points in image 1
    maybePts1 = pts1(maybeInliers, :);
    
    % Match points in image 2
    maybePts2 = pts2(maybeInliers, :);

    % Estimate fundamental matrix from these matches
    maybeF = epipolarGeometry(maybePts1, maybePts2);

    alsoInliers = [];
    notMaybeInliers = setdiff(1:npts, maybeInliers);
    
    for i=1:numel(notMaybeInliers)
        pt1 = pts1(notMaybeInliers(i), :);
        pt2 = pts2(notMaybeInliers(i), :);

        error = abs([pt1, 1] * maybeF * [pt2'; 1]);
        
        if (error < p.Results.errorThreshold)
            alsoInliers = [alsoInliers, notMaybeInliers(i)];
        end
    end

    if (numel(alsoInliers) > p.Results.alsoInliersThreshold)
        maybeAlso     = [maybeInliers, alsoInliers];
        maybeAlsoPts1 = pts1(maybeAlso,:);
        maybeAlsoPts2 = pts2(maybeAlso,:);

        F = epipolarGeometry(maybeAlsoPts1, maybeAlsoPts2);

        n     = numel(maybeAlso);
        err   = diag([maybeAlsoPts1, ones(n,1)] * F * [maybeAlsoPts2'; ones(1,n)]);
        error = sqrt(mean(err .^ 2));

        if (error < bestError)
            inliers = maybeAlso;
            bestF = F;
            bestError = error;
            fprintf("Iteration %4d | No. inliers: %5d | RMSE: %.3f\n", iter, ...
                numel(inliers), bestError);
        end

    end
end

inlierPts1 = pts1(inliers,:);
inlierPts2 = pts2(inliers,:);

switch (nargout) 
    case 1
        varargout{1} = inliers;
    case 3
        varargout{1} = inliers;
        varargout{2} = inlierPts1;
        varargout{3} = inlierPts2;
    case 4
        varargout{1} = inliers;
        varargout{2} = inlierPts1;
        varargout{3} = inlierPts2;
        varargout{4} = bestF;
    otherwise
        warning("epipolarConstraint has no return variable");
end

fprintf("Complete.\n")
end