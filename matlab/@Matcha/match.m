function [inlierPts1, inlierPts2, F, pts1, pts2] = match(obj, im2, ratio, pnorm, maxiter, errorThresh, inlierThresh)
% MATCH compute nearest neighbors and filter with epipolar constraint
%
% Other m-files required: 3d/epipolarRansac.m, 3d/epipolarGeometry.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

assert(obj.COMPUTED_DESCRIPTIONS, 'Matcha: Must compute descriptions first.');

fprintf("Matcha: Looking for nearest neighbors. ")
[idx, distance1, distance2] = nearestNeighbors(obj.keys.features, im2.keys.features, pnorm);

pts1 = [obj.keys.x, obj.keys.y];
pts2 = [im2.keys.x(idx), im2.keys.y(idx)];

% Ratio distance to closest and second closest
ratioIdx = find(distance1 < ratio * distance2);

pts1 = pts1(ratioIdx,:);
pts2 = pts2(ratioIdx,:);

numMatches = size(pts1, 1);
fprintf("Found %d matches.\n", numMatches);

[~, inlierPts1, inlierPts2, F] = epipolarRansac(pts1, pts2, maxiter, errorThresh, inlierThresh);

obj.COMPUTED_MATCHES = true;
end

