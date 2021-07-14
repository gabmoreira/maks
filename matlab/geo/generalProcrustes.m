function [model, rmse] = generalProcrustes(A, B)
% GENERALPROCRUSTES Computes translation and rotation that align A and B
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

assert(size(A,2) == 3, 'Argument A must be n x 3');
assert(size(B,2) == 3, 'Argument B must be n x 3');

centroidA = mean(A,1);
centroidB = mean(B,1);

ACentered = A - centroidA;
BCentered = B - centroidB;

model.t = reshape(centroidB - centroidA, 3, 1);

model.R = projectToSO3(ACentered' * BCentered);

rmse = sqrt(mean(vecnorm(ACentered * model.R - BCentered, 2, 2) .^ 2));

model.t = reshape(centroidB - centroidA * model.R, 3, 1);
end

