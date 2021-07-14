function [X,res] = triangulate(x1, x2, P1, P2)
% TRIANGULATE Compute 3D points from 2D matches and projection matrices
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

assert(size(x1,1) == size(x2,1), 'Sets of points must have the same size.');

% Store triangulated points
X = zeros(size(x1,1), 3);

% Store residuals
res = zeros(size(x1,1), 1);

for i=1:size(x1,1)
    A = [x1(i,2) * P1(3,:) - P1(2,:); 
         P1(1,:) - x1(i,1) * P1(3,:); 
         x2(i,2) * P2(3,:) - P2(2,:);
         P2(1,:) - x2(i,1) * P2(3,:);];
     
    [~,S,V] = svd(A);
    
    res(i) = min(diag(S));
    X(i,:) = V(1:3,end)' ./ V(4,end);
end

end

