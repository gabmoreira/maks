function [idx2, distance1, distance2] = nearestNeighbors(features1, features2, pnorm)
% NEARESTNEIGHBORS Computes nearest neighbors
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

idx1 = 1:size(features1, 1);
idx2 = zeros(size(idx1));

distance1 = zeros(size(idx1));
distance2 = zeros(size(idx1));

for i=1:size(features1, 1)
    % Select feature in frame 1
    feature = features1(i,:);
    
    % Distance from this feature to all others in 2
    d = vecnorm(features2 - feature, pnorm, 2);
    
    [distance1(i), idx2(i)] = min(d);
    
    d(idx2(i)) = inf;
    [distance2(i), ~] = min(d);
end
    
end

